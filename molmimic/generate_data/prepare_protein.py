import os, sys
import subprocess
import shutil
import gzip
import glob
import re
import time
import itertools as it
from collections import defaultdict
from multiprocessing.pool import ThreadPool
import logging

import pandas as pd
import dask
import dask.dataframe as dd
from joblib import Parallel, delayed

from Bio import SeqIO

from toil.job import JobFunctionWrappingJob

from toil.realtimeLogger import RealtimeLogger

from botocore.exceptions import ClientError

from molmimic.parsers.Electrostatics import run_pdb2pqr
from molmimic.parsers.SCWRL import run_scwrl
from molmimic.parsers.MODELLER import run_ca2model
from molmimic.parsers.CNS import Minimize
from molmimic.parsers.mmtf_spark import PdbToMmtfFull
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.job_utils import cleanup_ids, map_job_rv, map_job
from molmimic.generate_data.util import get_file, get_first_chain, get_all_chains, number_of_lines, \
    iter_unique_superfams, SubprocessChain, get_jobstore_name, is_ca_model

#Auto-scaling on AWS with toil has trouble finding modules? Heres the workaround
PDB_TOOLS = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdb_tools")

def setup_dask(num_workers):
    dask.config.set(scheduler='multiprocessing')
    dask.config.set(pool=ThreadPool(num_workers))

def extract_domain(pdb_file, pdb, chain, sdi, rslices, domNo, sfam_id, rename_chain=None, striphet=True, cath=True, work_dir=None):
    """Extract a domain from a protein structure and cleans the output to make
    it in standard PDB format. No information in changed or added

    Prepare a protein structure for use in molmimic. This will
    0) Unzip if gzipped
    1) Cleanup PDB file and add TER lines in between gaps
    2) Remove HETATMS
    3) Split Chains into separate files
    4) Each chain is then protonated and minimized. See 'run_single_chain'
    for more info

    Parameters
    ----------
    pdb_file : str
        Path to PDB file
    chain : str or None
        Chain ID to split out, protonate, and mininize. If None (default),
        all chains will be split out, protonated, and mininized.

    Returns
    -------
    If no chain was specified a list of all chains is returned with the
    path to the prepared file, the pdb code, and chain. If a chain is
    specified, 3 values are returned: the path to the prepared file,
    the pdb code, and chain. See 'run_single_chain' for info.
    """
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    if work_dir is None:
        work_dir = os.getcwd()

    if cath:
        domain_file = os.path.join(work_dir, "{}.pdb".format(sdi))
    else:
        domain_file = os.path.join(work_dir, "{}_{}_sdi{}_d{}.pdb".format(
            pdb, chain, int(sdi), domNo))


    open_fn = gzip.open if pdb_file.endswith(".gz") else open

    if pdb_file.endswith(".gz"):
        input = domain_file+".full"
        with gzip.open(pdb_file, 'rt') as zipf, open(input, "w") as pdbf:
            pdbf.write(zipf.read())
    else:
        input = pdb_file

    with open(input) as f:
        if f.read() == "":
            raise RuntimeError("Error processing PDB: {}".format(input))

    commands = [
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", input],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain)],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_striphet.py")],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
    ]

    if rslices is not None:
        commands += [[sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+rslices]
    commands += [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]]

    if rename_chain is not None:
        commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"),
            "-{}".format("1" if isinstance(rename_chain, bool) and rename_chain else rename_chain)])

    RealtimeLogger.info("{}".format(commands))
    with open(domain_file, "w") as output:
        SubprocessChain(commands, output)

    with open(domain_file) as f:
        content = f.read().rstrip()
        if content == "":
            RealtimeLogger.info(open(domain_file+".full").read())
            raise RuntimeError("Error processing PDB: {}".format(domain_file))

    if pdb_file.endswith(".gz"):
        os.remove(domain_file+".full")

    return domain_file

def prepare_domain(pdb_file, chain, work_dir=None, pdb=None, domainNum=None, sdi=None, sfam_id=None, job=None, cleanup=True):
    """Prepare a single domain for use in molmimic. This method modifies a PDB
    file by adding hydrogens with PDB2PQR (ff=parse, ph=propka) and minimizing
    using rosetta (lbfgs_armijo_nonmonotone with tolerance 0.001). Finally,
    the domain PDB file is cleaned so that can be parsed by simple PDB parsers.

    Method called during 'run_protein'

    Note: AssertError raised if pdb_file is gzipped or contains 2+ chains

    Parameters
    ----------
    pdb_file : str
        Path to PDB file of chain to prepare. Must not be gzipped or contain
        more than one chain.

    Returns
    -------
    cleaned_minimized_file : str
        Path to prepared PDB file
    pdb : str
        PDB ID
    chain : str
        Chain ID
    domain : 2-tuple
        Structural domain ID and domain Number.
    """
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    if pdb_file.endswith(".gz"):
        raise RuntimeError("Cannot be a gzip archive, try 'run_protein' instead")

    if work_dir is None:
        work_dir = os.getcwd()

    ##FIME: Not sure why this none..
    num_chains = len(get_all_chains(pdb_file))
    if not len(get_all_chains(pdb_file)) == 1:
        raise RuntimeError("Must contain only one chain. There {} chains {}".format(num_chains, pdb_file))

    pdb_path = os.path.dirname(pdb_file)
    prefix = pdb_file.split(".", 1)[0]

    #Add hydrogens and/or correct sidechains
    scwrl_file = None
    full_model_file = None
    propka_file = prefix+".propka"
    pdb2pqr_parameters = ["--chain"] #["--ph-calc-method=propka", "--chain", "--drop-water"]

    try:
        pqr_file = run_pdb2pqr(pdb_file, whitespace=False, ff="parse", parameters=pdb2pqr_parameters, work_dir=work_dir, job=job)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e1:
        #Might have failed to missing too many heavy atoms
        #Try again, but first add correct side chains
        try:
            scwrl_file = run_scwrl(pdb_file, work_dir=work_dir, job=job)
            pqr_file = run_pdb2pqr(scwrl_file, whitespace=False, ff="parse", parameters=pdb2pqr_parameters, work_dir=work_dir, job=job)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e2:
            if is_ca_model(pdb_file):
                #Run modeller to predict full atom model
                RealtimeLogger.info("Building CA model")
                try:
                    full_model_file = run_ca2model(pdb_file, chain, work_dir=work_dir, job=job)
                    pqr_file = run_pdb2pqr(full_model_file, whitespace=False, ff="parse", parameters=pdb2pqr_parameters, work_dir=work_dir, job=job)
                except (SystemExit, KeyboardInterrupt):
                    raise
                except Exception as e2:
                    raise
            else:
                #It really failed, skip it and warn
                raise
                raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs.".format(pdb_file) +
                                   "Most likeley reason for failing is that the structure is missing too many heavy atoms. {} or {}".format(e1, e2))

    try:
        with open(pqr_file) as f:
            pass
    except IOError:
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs.  \
Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    #Minimize using CNS
    minimized_file, _ = Minimize(pqr_file, work_dir=work_dir, job=job)
    commands = [
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_stripheader.py"), minimized_file],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(chain)],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
    ]

    cleaned_file = prefix+".pdb"
    with open(cleaned_file, "w") as cleaned:
        SubprocessChain(commands, cleaned)

    attempts = 0
    while number_of_lines(cleaned_file) == 0:
    	if attempts >= 10:
    		raise RuntimeError("Invalid PDB file")
        time.sleep(0.2)
        attempts += 1

    if cleanup:
        for f in [pqr_file, scwrl_file, full_model_file, minimized_file]:
            if f is not None and os.path.isfile(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    return cleaned_file

def process_domain(job, sdi, pdbFileStoreID, force_chain=None, force_rslices=None, force=False, work_dir=None, cleanup=True, cath=True, memory="12G", preemptable=True):
    if work_dir is None:
        if job:
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    in_store = IOStore.get("aws:us-east-1:molmimic-pdb")
    pdb_info_file = copy_pdb_h5(job, pdbFileStoreID)

    if cath:
        all_sdoms = pd.read_hdf(str(pdb_info_file), "table")
        sdom = all_sdoms[all_sdoms.cath_domain==sdi]
        all_domains = all_sdoms[all_sdoms.cath_domain.str.startswith(sdi[:5])].cath_domain.drop_duplicates()
        if len(all_domains) == 1 and force_rslices is None:
            #Use full chain
            rslices = [":"]
        elif force_rslices is not None:
            rslices = force_rslices
        else:
            rslices = ["{}:{}".format(st, en) for st, en in \
                sdom.sort_values("nseg")[["srange_start", "srange_stop"]]\
                .drop_duplicates().itertuples(index=False)]

        sfam_id = sdom.iloc[0].cathcode
        pdb = sdom.iloc[0].pdb
        chain = force_chain if force_chain is not None else sdi[4]
        domNo = sdom.iloc[0].domain

        key = "{}/{}.pdb".format(sfam_id.replace(".", "/"), sdi)
        out_store = IOStore.get("aws:us-east-1:molmimic-cath-structures")

        RealtimeLogger.info("CATH_DOMAIN {} CHAIN {}".format(sdi, chain))

    else:


        all_sdoms = pd.read_hdf(str(pdb_info_file), "merged") #, where="sdi == {}".format(sdi))

        sdom = all_sdoms[all_sdoms["sdi"]==float(sdi)]
        if sdom.shape[0] == 0:
            RealtimeLogger.info("SDI {} does not exist".format(sdi))
            return None, None, None

        sfam_id = sdom.iloc[0].sfam_id
        pdb = sdom.iloc[0].pdbId
        chain = force_chain if force_chain is not None else sdom.iloc[0].chnLett
        domNo = sdom.iloc[0].domNo

        _whole_chain = all_sdoms[(all_sdoms["pdbId"]==pdb)&(all_sdoms["chnLett"]==chain)]
        all_domains = _whole_chain[_whole_chain["whole_chn"]!=1.0]["sdi"].drop_duplicates()
        whole_chain = _whole_chain[_whole_chain["whole_chn"]==1.0]["sdi"].drop_duplicates()

        if len(all_domains) == 1 and force_rslices is None:
            #Use full chain
            rslices = [":"]
        elif force_rslices is not None:
            rslices = force_rslices
        else:
            rslices = ["{}:{}".format(st, en) for st, en in sdom[["from", "to"]]\
                .drop_duplicates().itertuples(index=False)]


        out_store = IOStore.get("aws:us-east-1:molmimic-full-structures")

        get_key = lambda f, p, c, s, d: "{}/{}/{}_{}_sdi{}_d{}.pdb".format(int(f),
            p[1:3].lower(), p.upper(), c, s, d)

        key = get_key(int(sfam_id), pdb, chain, sdi, domNo)

    if not force and out_store.exists(key):
        pdb_path = os.path.join(work_dir, key)
        if not os.path.isdir(os.path.dirname(pdb_path)):
            os.makedirs(os.path.dirname(pdb_path))
        out_store.read_input_file(key, pdb_path)

        #Correct ca_alignments
        needs_update = False
        out_store.read_input_file(key, pdb_path)
        with open(pdb_path) as dom, open(pdb_path) as raw:
            atom1, atom2 = next(dom), next(raw)
            if atom1[17:20] != atom2[17:20]:
                needs_update = True

        if not needs_update and not out_store.exists(key+".prepared"):
            return pdb_path, key, sfam_id

    #pdb_path = os.path.join("by_superfamily", str(sfam_id), pdb[1:3].upper())
    #domain_file = os.path.join(pdb_path, "{}_{}_sdi{}_d{}.pdb".format(pdb, chain, sdi, domNo))

    pdb_file_base = os.path.join(pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    pdb_file = os.path.join(work_dir, "pdb{}.ent.gz".format(pdb.lower()))

    RealtimeLogger.info("RUNNING DOMAIN {} {} {} {} {}".format(pdb_file, pdb, chain, sdi, domNo))

    try:
        #Download PDB archive from JobStore
        RealtimeLogger.info("AWS GET: {}; Save to: {}".format(pdb_file_base, pdb_file))
        in_store.read_input_file(pdb_file_base, pdb_file)
    except ClientError:
        try:
            in_store.read_input_file("obsolete/"+pdb_file_base, pdb_file)
        except ClientError:
            from Bio.PDBList import PDBList
            import gzip
            obsolete = False
            pdb_file = PDBList().retrieve_pdb_file(pdb, pdir=work_dir, file_format="pdb")
            if not os.path.isfile(pdb_file):
                obsolete = True
                pdb_file = PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=work_dir, file_format="pdb")
                if not os.path.isfile(pdb_file):
                    raise IOError("{} not found".format(r))

            with open(pdb_file, 'rb') as f_in, gzip.open(pdb_file+'.gz', 'wb') as f_out:
                f_out.writelines(f_in)

            in_store.write_output_file(pdb_file, "{}{}".format("obsolete/" if obsolete else "", pdb_file_base))

            try:
                os.remove(pdb_file+".gz")
            except OSError:
                pass


        #Extract domain; cleaned but atomic coordinates not added or changed
        domain_file = extract_domain(pdb_file, pdb, chain, sdi, rslices, domNo, sfam_id, work_dir=work_dir)
        RealtimeLogger.info("Finished extracting domain: {}".format(domain_file))

        with open(domain_file) as f:
            pass

        domain_file_base = os.path.basename(domain_file)
        out_store.write_output_file(domain_file, key) #os.path.join(str(int(sfam_id)), pdb[1:3].lower(), domain_file_base))

        try:
            prepared_file = prepare_domain(domain_file, chain, work_dir=work_dir, job=job)
            #prepared_key = os.path.join(str(int(sfam_id)), pdb[1:3].lower(), os.path.basename(prepared_file))
            out_store.write_output_file(prepared_file, key+".prepared")
            RealtimeLogger.info("Wrote output file {}".format(key+".prepared"))
        except RuntimeError as e:
            RealtimeLogger.info(str(e))
            raise
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        RealtimeLogger.info("Cannot process {}.{}.d{} ({}), error: {}".format(pdb, chain, domNo, sdi, e))
        return

    if cleanup:
        for f in (pdb_file, domain_file, prepared_file):
            if os.path.isfile(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    return prepared_file, domain_file_base, sfam_id

def convert_pdb_to_mmtf(job, sfam_id, jobStoreIDs=None, clustered=True, preemptable=True):
    raise NotImplementedError()

    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    clustered = "clustered" if clustered else "full"

    pdb_path = os.path.join(work_dir, "pdb")
    if not os.path.isdir(pdb_path):
        os.makedirs(pdb_path)

    #Download all with same sfam
    if jobStoreIDs is None:
        in_store = IOStore.get("{}:molmimic-{}-structures".format(prefix), clustered)
        for f in in_store.list_input_directory(sfam_id):
            if f.endswith(".pdb"):
                in_store.read_input_file(f, os.path.join(work_dir, f))
    else:
        for jobStoreID in jobStoreIDs:
            job.fileStore.readGlobalFile(fileStoreID, userPath=pdb_path)

    PdbToMmtfFull(pdb_path, mmtf_path, work_dir=work_dir, job=job)

    out_store = IOStore.get("{}:molmimic-{}-mmtf".format(prefix, clustered))
    out_store.write_output_directory(mmtf_path, sfam_id)

def create_data_loader(job, sfam_id, preemptable=True):
    """Create H5 for Molmimic3dCNN to read

    Note: move this somewhere else
    """
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]

    pdb_path = os.path.join(work_dir, "pdb")
    if not os.path.isdir(pdb_path):
        os.makedirs(pdb_path)

    id_format = re.compile("^([A-Z0-9]{4})_([A-Za-z0-9]+)_sdi([0-9]+)_d([0-9]+)$")

    #Get all with keys same sfam, but do not download

    in_store = IOStore.get("{}:molmimic-clustered-structures".format(prefix))
    keys = [id_format.match(f).groups() for f in in_store.list_input_directory(sfam_id) \
        if f.endswith(".pdb") and id_format.match(f)]

    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)))
    clusters_file = os.path.join(pdb_path, "{}_nr.fasta".format(int(sfam_id)))


    try:
        pdb, chain, sdi, domain = list(zip(*[id_format.match(seq.id[:-2]).groups() \
            for s in SeqIO.parse(clusters_file, "fasta")]))
    except ValueError:
        RealtimeLogger.info("Unable to create data loading file for {}.".format(sfam_id))
        return

    domains = pd.DataFrame({"pdb":pdb, "chain":chain, "domNo":domain, "sdi":sdi})

    data_loader = os.path.join(pdb_path, "{}.h5".format(int(sfam_id)))
    domains.to_hdf(str(data_loader), "table", complevel=9, complib="bzip2")

def copy_pdb_h5(job, path_or_pdbFileStoreID):
    if os.path.isfile(path_or_pdbFileStoreID):
        return path_or_pdbFileStoreID
    else:
        work_dir = job.fileStore.getLocalTempDir()
        sdoms_file = os.path.join(work_dir, "PDB.h5")

        with job.fileStore.readGlobalFileStream(path_or_pdbFileStoreID) as fs_sdoms, open(sdoms_file, "w") as f_sdoms:
            for line in fs_sdoms:
                f_sdoms.write(line)

        return sdoms_file

def process_sfam(job, sfam_id, pdbFileStoreID, cath, clusterFileStoreID, further_parallize=False, cores=1):
    work_dir = job.fileStore.getLocalTempDir()

    if cath:
        clusters_file = get_file(job, "S100.txt", clusterFileStoreID, work_dir=work_dir)
        all_sdoms = pd.read_csv(clusters_file, delim_whitespace=True, header=None,
            usecols=[0, 1, 2, 3, 4], names=["cath_domain", "C", "A", "T", "H"])
        RealtimeLogger.info("ALL SDOMS {}".format(all_sdoms.head()))
        groups = all_sdoms.groupby(["C", "A", "T", "H"])
        RealtimeLogger.info("Running sfam {}".format(sfam_id))
        sdoms = groups.get_group(tuple(sfam_id))["cath_domain"].drop_duplicates().dropna()
    else:
        sdoms_file = copy_pdb_h5(job, pdbFileStoreID)
        sdoms = pd.read_hdf(str(sdoms_file), "merged")
        sdoms = sdoms[sdoms["sfam_id"]==float(sfam_id)]["sdi"].drop_duplicates().dropna()
    #sdoms = sdoms[:1]
    if further_parallize:
        if cores > 2:
            #Only makes sense for slurm or other bare-matal clsuters
            setup_dask(cores)
            d_sdoms = dd.from_pandas(sdoms, npartitions=cores)
            RealtimeLogger.info("Running sfam dask {}".format(sdoms))
            processed_domains = d_sdoms.apply(lambda row: process_domain(job, row.sdi,
                sdoms_file), axis=1).compute()
        else:
            map_job(job, process_domain, sdoms, pdbFileStoreID)
    else:
        for sdom in sdoms:
            process_domain(job, sdom, pdbFileStoreID)


def post_process_sfam(job, sfam_id, jobStoreIDs):
    cluster(job, sfam_id, jobStoreIDs)

    if False:
        convert_pdb_to_mmtf(job, sfam_id, jobStoreIDs)
        create_data_loader(job, sfam_id, jobStoreIDs)

def start_toil(job, cath=True):
    """Start the workflow to process PDB files"""
    work_dir = job.fileStore.getLocalTempDir()


    if cath:
        in_store = IOStore.get("aws:us-east-1:molmimic-cath")
        sdoms_file = os.path.join(work_dir, "cath-domain-description-file-small.h5")
        in_store.read_input_file("cath-domain-description-file-small.h5", sdoms_file)

        #Add pdb info into local job store
        pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

        clusters_file = os.path.join(work_dir, "cath-domain-list-S100.txt")
        in_store.read_input_file("cath-domain-list-S35.txt", clusters_file)

        with open(clusters_file) as f:
            sfams = list(set([tuple(map(int, l.split()[1:5])) for l in f if l and not l.startswith("#")]))

        clusterFileStoreID = job.fileStore.writeGlobalFile(clusters_file)

    else:
        in_store = IOStore.get("aws:us-east-1:molmimic-ibis")

        #Download PDB info
        sdoms_file = os.path.join(work_dir, "PDB.h5")
        in_store.read_input_file("PDB.h5", sdoms_file)

        #Add pdb info into local job store
        pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

        #Get all unique superfamilies
        sdoms = pd.read_hdf(str(sdoms_file), "merged")

        # skip_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "keep.csv")
        # if os.path.isfile(skip_file):
        #     skip = pd.read_csv(skip_file)
        #     sdoms = sdoms[sdoms["sdi"].isin(skip["sdi"])]
        #     RealtimeLogger.info("SKIPPING {} sdis; RUNIING {} sdis".format(skip.shape[0], sdoms.shape[0]))
        #
        sfams = sdoms["sfam_id"].drop_duplicates().dropna()

        clusterFileStoreID = None
    #sfams = sfams[:1]
    #sfams = ["653504"]

    # max_cores = job.fileStore.jobStore.config.maxCores if \
    #     job.fileStore.jobStore.config.maxCores > 2 else \
    #     job.fileStore.jobStore.config.defaultCores

    max_cores = job.fileStore.jobStore.config.defaultCores
    #Add jobs for each sdi
    job.addChildJobFn(map_job, process_sfam, sfams, pdbFileStoreID, cath, clusterFileStoreID,
        cores=max_cores)

    #Add jobs for to post process each sfam
    #job.addFollowOnJobFn(map_job, post_process_sfam, sfams, pdbFileStoreID,
    #    cores=max_cores)

    del sfams
    os.remove(sdoms_file)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    detector = logging.StreamHandler()
    logging.getLogger().addHandler(detector)

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
