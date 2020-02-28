import os, sys
import subprocess
import shutil
import gzip
import glob
import re
import time
import traceback
import itertools as it
from collections import defaultdict
from multiprocessing.pool import ThreadPool
import logging

import pandas as pd
from joblib import Parallel, delayed

from Bio import SeqIO

from toil.job import JobFunctionWrappingJob

from toil.realtimeLogger import RealtimeLogger

from botocore.exceptions import ClientError

from molmimic.parsers.Electrostatics import Pdb2pqr
from molmimic.parsers.SCWRL import SCWRL
from molmimic.parsers.MODELLER import MODELLER
#from molmimic.parsers.cns import CNSMinimize

from molmimic.util.iostore import IOStore
from molmimic.util.toil import map_job
from molmimic.util.hdf import get_file
from molmimic.util import SubprocessChain, safe_remove
from molmimic.util.pdb import get_first_chain, get_all_chains, PDB_TOOLS
from molmimic.util.cath import download_cath_domain

from molmimic.generate_data import data_stores

class PrepareProteinError(RuntimeError):
    def __init__(self, cath_domain, stage, message, errors=None, *args, **kwds):
        super().__init__(*args, **kwds)
        self.cath_domain = cath_domain
        self.stage = stage
        self.message = message
        self.errors = errors if isinstance(errors, list) else []

    def __str__(self):
        return "Error during {}: {}\nErrors:\n".format(self.stage, self.message,
            "\n".join(map(str, self.errors)))

    def save(self, store=None):
        if store is None:
            store = data_stores.prepared_cath_structures
        fail_file = "{}.{}".format(self.cath_domain, self.stage)
        with open(fail_file, "w") as f:
            print(self.message, file=f)
            print(self.errors, file=f)
        store.write_output_file(fail_file, "errors/"+os.path.basename(fail_file))
        safe_remove(fail_file)

def extract_domain(pdb_file, cath_domain, sfam_id, rename_chain=None,
  striphet=True, rslices=None, work_dir=None):
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

    #Name out output domain file
    domain_file = os.path.join(work_dir, "{}.pdb".format(cath_domain))

    files_to_remove = []

    if pdb_file.endswith(".gz"):
        #Open gzip file
        input = domain_file+".full"
        with gzip.open(pdb_file, 'rt') as zipf, open(input, "w") as pdbf:
            pdbf.write(zipf.read())

        #Remove ungzipped file at end
        files_to_remove.append(input)
    else:
        input = pdb_file

    #Make sure input is not empty
    with open(input) as f:
        if f.read() == "":
            raise RuntimeError("Error processing PDB: {}".format(input))

    chain = cath_domain[4]

    commands = [
        #Pick first model
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", input],
        #Select desired chain
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain)],
        #Remove altLocs
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
        #Remove HETATMS
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_striphet.py")],
    ]

    prep_steps = ["pdb_selmodel.py -1 {}".format(input), "pdb_selchain.py -{}".format(chain),
        "pdb_delocc.py", "pdb_striphet.py"]

    if rslices is not None:
        #Slice up chain with given ranges
        commands += [[sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+rslices]
        prep_steps += ["pdb_rslice.py {}".format(" ".join(rslices))]

    #Make it tidy
    commands += [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]]
    prep_steps += ["pdb_tidy.py"]

    if rename_chain is not None:
        #Rename chain to given name if set in arguments
        new_chain = "-{}".format("1" if isinstance(rename_chain, bool) and rename_chain \
            else rename_chain)
        commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"),
            new_chain])
        prep_steps += ["pdb_chain.py {}".format(new_chain)]

    #Run all commands
    with open(domain_file, "w") as output:
        SubprocessChain(commands, output)

    #Make sure output domain_file is not empty
    with open(domain_file) as f:
        content = f.read().rstrip()
        if content == "":
            raise RuntimeError("Error processing PDB: {}".format(domain_file))

    safe_remove(files_to_remove)

    return domain_file, prep_steps

def prepare_domain(pdb_file, chain, cath_domain, sfam_id=None,
  perform_cns_min=False, work_dir=None, job=None, cleanup=True):
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
        raise RuntimeError("Cannot be a gzip archive, try 'extract_domain' instead")

    if work_dir is None:
        work_dir = os.getcwd()

    num_chains = len(get_all_chains(pdb_file))
    if not len(get_all_chains(pdb_file)) == 1:
        raise RuntimeError("Must contain only one chain. There are {} chains in {}".format(
            num_chains, pdb_file))

    prefix = pdb_file.split(".", 1)[0]

    pdb2pqr = Pdb2pqr(work_dir=work_dir, job=job)
    scwrl = SCWRL(work_dir=work_dir, job=job).fix_rotamers
    modeller = MODELLER(work_dir=work_dir, job=job).remodel_structure

    errors = []
    files_to_remove = []
    for attempt, fixer in enumerate((None, scwrl, modeller)):
        try:
            #If failed 1st time, add correct sidechain rotamers (SCWRL)
            #If failed 2nd time, turn CA models into full atom models (MODELLER)
            fixed_pdb = fixer(pdb_file) if fixer is not None else pdb_file
            if fixer is not None:
                #Remove "fixed" pdb file at end
                files_to_remove.append(fixed_pdb)
        except (SystemExit, KeyError):
            raise
        except Exception as e:
            tb = traceback.format_exc()
            errors.append(tb)
            RealtimeLogger.info("Fixer {} failed: {}".format(fixer.__class__.__name__, tb))
            continue

        try:
            #Protonate PDB, Minimize structure, and assign partial charges
            protonated_pdb = pdb2pqr.create_tidy_pqr(
                fixed_pdb, keep_occ_and_b=True, ff="parse", chain=True,)
            break
        except Exception as error:
            #Failed, try again with different fixer
            tb = traceback.format_exc()
            errors.append(tb)
            RealtimeLogger.info("Pdb2pqr failed: {}".format(tb))

    else:
        #Completely failed, raise error
        raise PrepareProteinError(cath_domain, "prepare", "Unable to protonate" + \
            "{} using pdb2pqr. Please check pdb2pqr error logs.".format(pdb_file) + \
            "Most likeley reason for failing is that the structure is missing " + \
            "too many heavy atoms.", errors)

    if perform_cns_min:
        #Remove pdb2pqr file at end
        files_to_remove.append(protonated_pdb)

        #Perform CNS minimization
        minimize = CNSMinimize(work_dir=work_dir, job=job)
        protonated_pdb = minimize(protonated_pdb)

    #Remove pdb2pqr file (or CNS file) at end
    files_to_remove.append(protonated_pdb)

    #Clean
    commands = [
        #Remove header
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_stripheader.py"), protonated_pdb],
        #Rename chain to given chain
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(chain)],
        #Make it tidy
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
    ]

    cleaned_file = prefix+".pdb"
    with open(cleaned_file, "w") as cleaned:
        SubprocessChain(commands, cleaned)

    if cleanup:
        safe_remove(files_to_remove)

    prep_steps = [fixer.__class__.__name__.rsplit(".")[-1]] if fixer is not None \
        else []
    prep_steps += ["pdb2pqr", "pdb_stripheader.py", "pdb_chain.py -{}".format(chain),
        "pdb_tidy.py"]

    return cleaned_file, prep_steps

def _process_domain(job, cath_domain, cathcode, cathFileStoreID=None, force_chain=None,
  force_rslices=None, force=False, work_dir=None, get_from_pdb=False, cleanup=True,
  memory="12G", preemptable=True):
    """Process domains by extracting them, protonating them, and then minimizing.
    PrepareProteinError is raised on error"""
    if work_dir is None:
        if job and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    files_to_remove = []

    #Get cath domain file
    cath_key = "{}/{}.pdb".format(cathcode.replace(".", "/"), cath_domain)
    if data_stores.prepared_cath_structures.exists(cath_key):
        return

    try:
        #Download cath domain from s3 bucket or cath api
        domain_file = download_cath_domain(cath_domain, cathcode, work_dir=work_dir)
        files_to_remove.append(domain_file)
        assert os.path.isfile(domain_file), "Domain file not found: {}".format(domain_file)
    except (SystemExit, KeyboardInterrupt):
        raise
    except KeyError as e:
        if get_from_pdb:
            raise NotImplementedError("Download from PDB is not finishied")
            all_sdoms = pd.read_hdf(str(pdb_info_file), "table")
            sdom = all_sdoms[all_sdoms.cath_domain==sdi]
            all_domains = all_sdoms[all_sdoms.cath_domain.str.startswith(sdi[:5])].cath_domain.drop_duplicates()
            if len(all_domains) == 1 and force_rslices is None:
                #Use full chain
                rslices = None
            elif force_rslices is not None:
                rslices = force_rslices
            else:
                rslices = ["{}:{}".format(st, en) for st, en in \
                    sdom.sort_values("nseg")[["srange_start", "srange_stop"]]\
                    .drop_duplicates().itertuples(index=False)]

            pdb_file_base = os.path.join(pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
            pdb_file = os.path.join(work_dir, "pdb{}.ent.gz".format(pdb.lower()))

        else:
            #Cannot download
            tb = traceback.format_exc()
            raise PrepareProteinError(cath_domain, "download", tb)
    except Exception as e:
        tb = traceback.format_exc()
        raise PrepareProteinError(cath_domain, "unk_download", tb)

    if not data_stores.prepared_cath_structures.exists(cath_key+".raw"):
        #Extract domain; cleaned but atomic coordinates not added or changed
        try:
            domain_file, prep_steps = extract_domain(domain_file, cath_domain, cathcode,
                work_dir=work_dir)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            tb = traceback.format_exc()
            raise PrepareProteinError(cath_domain, "extract", tb)

        #Write raw domain file to store
        data_stores.prepared_cath_structures.write_output_file(domain_file, cath_key+".raw")

        #Write preperation steps
        prep_steps_file = os.path.join(work_dir, "{}.raw.prep".format(cath_domain))
        with open(prep_steps_file, "w") as fh:
            for step in prep_steps:
                print(step, file=fh)
        data_stores.prepared_cath_structures.write_output_file(prep_steps_file,
            cath_key+".raw.prep")

        files_to_remove.append(domain_file)
        RealtimeLogger.info("Finished extracting domain: {}".format(domain_file))

    chain = cath_domain[4]

    #Protonate and minimize, raises PrepareProteinError on error
    prepared_file, prep_steps = prepare_domain(domain_file, chain, cath_domain,
        work_dir=work_dir, job=job)

    if os.path.getsize(prepared_file) == 0:
        raise PrepareProteinError(cath_domain, "empty_file", "")

    #Write prepared domain file to store
    data_stores.prepared_cath_structures.write_output_file(prepared_file, cath_key)
    files_to_remove.append(prepared_file)
    RealtimeLogger.info("Finished preparing domain: {}".format(domain_file))

    #Write preperation steps
    prep_steps_file = os.path.join(work_dir, "{}.prep".format(cath_domain))
    with open(prep_steps_file, "w") as fh:
        for step in prep_steps:
            print(step, file=fh)
    data_stores.prepared_cath_structures.write_output_file(prep_steps_file,
        cath_key+".prep")

    if cleanup:
        safe_remove(files_to_remove)

def process_domain(job, cath_domain, cathcode, cathFileStoreID=None, force_chain=None,
  force_rslices=None, force=False, work_dir=None, get_from_pdb=False, cleanup=True,
  memory="12G", preemptable=True):
    try:
        try:
            _process_domain(job, cath_domain, cathcode, cathFileStoreID=cathFileStoreID,
                force_chain=force_chain, force_rslices=force_rslices, force=force,
                work_dir=work_dir, get_from_pdb=get_from_pdb, cleanup=cleanup,
                memory=memory, preemptable=preemptable)
        except (SystemExit, KeyboardInterrupt):
            raise
        except PrepareProteinError as e:
            e.save()
        except:
            tb = traceback.format_exc()
            raise PrepareProteinError(cath_domain, "unk_error", tb)
    except (SystemExit, KeyboardInterrupt):
        raise
    except PrepareProteinError as e:
        e.save()

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
    safe_remove(sdoms_file)

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
