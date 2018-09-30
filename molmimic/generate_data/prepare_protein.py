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

import pandas as pd
import dask
import dask.dataframe as dd
from joblib import Parallel, delayed

from Bio import SeqIO

from toil.job import JobFunctionWrappingJob

from molmimic.parsers.Electrostatics import run_pdb2pqr
from molmimic.parsers.CNS import Minimize
from molmimic.parsers.mmtf_spark import PdbToMmtfFull
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.job_utils import cleanup_ids, map_job_rv, map_job
from molmimic.generate_data.util import data_path_prefix, get_structures_path, \
    get_features_path, get_first_chain, get_all_chains, number_of_lines, \
    iter_unique_superfams, SubprocessChain, get_jobstore_name

def setup_dask(num_workers):
    dask.config.set(scheduler='multiprocessing')
    dask.config.set(pool=ThreadPool(num_workers))

def extract_domain(pdb_file, pdb, chain, sdi, rslices, domNo, sfam_id, rename_chain=None, striphet=True, work_dir=None):
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

    #pdb_path = os.path.join(work_dir, "by_superfamily", str(int(sfam_id)), pdb[1:3].upper())
    domain_file = os.path.join(work_file, "{}_{}_sdi{}_d{}.pdb.extracted".format(pdb, chain, int(sdi), domNo))

    # if not os.path.exists(pdb_path):
    #     os.makedirs(pdb_path)

    open_fn = gzip.open if pdb_file.endswith(".gz") else open

    if pdb_file.endswith(".gz"):
        input = domain_file+".full"
        with gzip.open(pdb_file, 'rt') as zipf, open(input, "w") as pdb:
            pdb.write(zipf.read())
    else:
        input = pdb_file

    commands = [
        [sys.executable, "-m", "pdb-tools.pdb_selmodel", "-1", input],
        [sys.executable, "-m", "pdb-tools.pdb_selchain", "-{}".format(chain)],
        [sys.executable, "-m", "pdb-tools.pdb_delocc"],
        [sys.executable, "-m", "pdb-tools.pdb_rslice"]+rslices,
        [sys.executable, "-m", "pdb-tools.pdb_striphet"],
        [sys.executable, "-m", "pdb-tools.pdb_tidy"]
    ]

    if rename_chain is not None:
        commands.append([sys.executable, "-m", "pdb-tools.pdb_chain.py",
            "-{}".format("1" if isinstance(rename_chain, bool) and rename_chain else rename_chain)])

    with open(domain_file, "w") as output:
        SubprocessChain(commands, output)

    if pdb_file.endswith(".gz"):
        os.remove(domain_file+".full")

    return domain_file

def prepare_domain(pdb_file, work_dir=None, pdb=None, chain=None, domainNum=None, sdi=None, sfam_id=None):
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

    #Add hydrogens
    propka_file = prefix+".propka"

    try:
        parameters = ["--ph-calc-method=propka", "--chain", "--drop-water"]
        pqr_file = run_pdb2pqr(pdb_file, whitespace=False, ff="parse", parameters=parameters, work_dir=work_dir)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs. \
Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    try:
        with open(pqr_file) as f:
            pass
    except IOError:
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs.  \
Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    #Minimize useing CNS
    minimized_file, score_file = Minimize(pqr_file, work_dir=work_dir)

    commands = [
        [sys.executable, "-m", "pdb-tools.pdb_stripheader", minimized_file],
        [sys.executable, "-m", "pdb-tools.pdb_chain", "-{}".format(chain)],
        [sys.executable, "-m", "pdb-tools.pdb_tidy"]
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

    return cleaned_file

def process_domain(job, sdi, pdbFileStoreID, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    pdb_info_file = job.fileStore.readGlobalFile(pdbFileStoreID)
    all_sdoms = pd.read_hdf(unicode(pdb_info_file), "merged")

    sdom = all_sdoms[all_sdoms["sdi"]==float(sdi)]
    if sdom.shape[0] == 0:
        job.log("SDI {} does not exist".format(sdi))
        return

    sfam_id = sdom.iloc[0].sfam_id
    pdb = sdom.iloc[0].pdbId
    chain = sdom.iloc[0].chnLett
    domNo = sdom.iloc[0].domNo

    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-pdb".format(prefix))
    out_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    pdb_path = os.path.join("by_superfamily", str(sfam_id), pdb[1:3].upper())
    domain_file = os.path.join(pdb_path, "{}_{}_sdi{}_d{}.pdb".format(pdb, chain, sdi, domNo))
    # if out_store.exists(domain_file):
    #     return

    job.log("RUNNING DOMAIN {} {} {} {} {}".format(pdb_file, pdb, chain, sdi, domNo))

    pdb_file_base = os.path.join(pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    pdb_file = os.path.join(work_dir, pdb_file_base)

    try:
        #Download PDB archive from JobStore
        in_store.read_input_file(pdb_file_base, pdb_file)

        #Multiple ranges for one sdi. Might be domain swaped?
        rslices = ["{}:{}".format(st, en) for st, en in sdom[["from", "to"]].drop_duplicates().itertuples(index=False)]

        #Extract domain; cleaned but atomic coordinates not added or changed
        domain_file = extract_domain(pdb_file, pdb, chain, sdi, rslices, domNo, sfam_id, work_dir=work_dir)
        domain_file_base = os.path.basename(domain_file)
        out_store.write_output_file(domain_file, os.path.join(str(int(sfam_id)), pdb[1:3].lower(), domain_file_base))

        try:
            prepared_file = prepare_domain(domain_file, work_dir=work_dir)
            out_store.write_output_file(prepared_file, os.path.join(str(int(sfam_id)), pdb[1:3].lower(), os.path.basename(prepared_file)))
            jobStoreID, job.fileStore.writeGlobalFile(prepared_file)
        except RuntimeError as e:
            job.log(str(e))
            pass
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        job.log("Cannot process {}.{}.d{} ({}), error: {}".format(pdb, chain, domNo, sdi, e))

    return prepared_file, domain_file_base, sfam_id

def cluster(job, sfam_id, jobStoreIDs, pdbFileStoreID, id=0.95, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-clustered-structures".format(prefix))

    sdoms_file = job.fileStore.readGlobalFile(pdbFileStoreID)

    sdoms = pd.read_hdf(unicode(sdoms_file), "merged")
    sdoms = sdoms[sdoms["sfam_id"]==sfam_id]
    sdoms = sdoms[["pdbId", "chnLett", "sdi", "domNo"]].drop_duplicates().dropna()

    #Save all domains to fasta
    domain_fasta = os.path.join(work_dir, "{}.fasta".format(int(sfam_id)))
    domain_ids = {}
    with open(domain_fasta, "w") as fasta:
        for jobStoreID, pdb_fname, _ in jobStoreIDs:
            with job.fileStore.readGlobalFileStream(jobStoreID) as f:
                try:
                    seq = subprocess.check_output([sys.executable, "-m", \
                        "pdb-tools.pdb_toseq.py"], stdin=f)
                    fasta.write(">{}\n{}\n".format(pdb_fname, "\n".join(seq.splitlines()[1:])))
                    domain_ids[pdb_fname] = jobStoreID
                except (KeyboardInterrupt, SystemExit):
                    raise
                except Exception as e:
                    job.log("Error getting fasta for : {} {}".format(sfam_id, fname))
                    pass

    # d_sdoms = dd.from_pandas(sdoms, npartitions=cores)
    # d_sdoms.apply(lambda row: subprocess.call([sys.executable,
    #     os.path.join(PDB_TOOLS, "pdb_toseq.py"), f])
    #
    # process_domain(job, dataset_name, row.sdi) \
    #     if not os.path.isfile(os.path.join(
    #         PDB_PATH, dataset_name, "by_superfamily", str(sfam_id),
    #         row.pdbId[1:3].upper(), "{}_{}_sdi{}_d{}.pdb".format(
    #             row.pdbId, row.chnLett, row.sdi, row.domNo
    #         )
    #     )) else None, axis=1).compute()


    #Save all domains to fasta
    # domain_fasta = os.path.join(pdb_path, "{}.fasta".format(int(sfam_id)))
    # with open(domain_fasta, "w") as fasta:
    #     for f in glob.glob(os.path.join(pdb_path, "*", "*.pdb")):
    #         subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_toseq.py"),
    #             f], stdout=fasta, stderr=subprocess.PIPE)

    #Order domains by resolution so the ones with the highest resolutions are centroids
    resolutions = pd.read_hdf(unicode(sdoms_file), "resolu")

    try:
        pdbs, ids, sequences = zip(*[(s.id.split("_", 1)[0].upper(), s.id, str(s.seq)) \
            for s in SeqIO.parse(domain_fasta, "fasta")])
    except ValueError:
        job.log("Unable to cluster {}. No PDBs passed the protonation/minization steps.".format(sfam_id))
        return

    domains = pd.DataFrame({"pdbId":pdbs, "domainId":ids, "sequence":sequences})
    domains = pd.merge(domains, resolutions, how="left", on="pdbId")
    xray = domains[domains["resolution"] >= 0.].sort_values("resolution")
    nmr = domains[domains["resolution"] < 0.]

    with open(domain_fasta, "w") as f:
        for row in it.chain(xray.itertuples(index=False), nmr.itertuples(index=False)):
            print >> f, ">{} [resolution={}]\n{}".format(row.domainId, row.resolution, row.sequence)

    sfam_key = "{0}/{0}.fasta".format(int(sfam_id))
    out_store.write_output_file(domain_fasta, sfam_key)

    clusters_file, uclust_file = run_usearch(["-cluster_fast",
        "{i}"+domain_fasta, "-id", str(id),
        "-centroids", "{{out}}{}_clusters.uc".format(int(sfam_id)),
        "-uc", "{{o}}{}_clusters.uc".format(int(sfam_id))])

    #Convert uclust to h5
    uclust = pd.read_table(unicode(uclust_file), comment="#", header=None, names=[
        "record_type",
        "cluster",
        "length",
        "pctId",
        "strand",
        "unk1",
        "unk2",
        "alignment",
        "label_query",
        "label_target"
    ])
    del uclust["unk1"]
    del uclust["unk2"]
    hdf_base = "{}_clusters.h5".format(int(sfam_id))
    hdf_file = os.path.join(work_dir, hdf_base)
    uclust.to_hdf(unicode(hdf_file), "table", complevel=9, complib="bzip2")
    out_store.write_output_file(hdf_file, "{}/{}".format(int(sfam_id), hdf_base))
    os.remove(uclust_file)

    #Upload clustered pdbs
    clustered_pdbs = []
    for seq in SeqIO.parse(clusters_file, "fasta"):
        pdb_base = seq.id
        jobStoreID = domain_ids[pdb_base]
        pdb_key = "{}/{}/{}".format(int(sfam_id), pdb_base[1:3], pdb_base)
        pdb_file = job.fileStore.readGlobalFile(jobStoreID)
        out_store.write_output_file(pdb_file, pdb_key)
        os.remove(pdb_file)
        clustered_pdbs.append((jobStoreID, pdb_base, sfam_id))

        #Remove jobStoreID for list so the actual file won't be removed
        del jobStoreIDs[jobStoreID.index(jobStoreID)]

    #Delete all reduant pdb files
    cleanup_ids(jobStoreIDs)

    return clustered_pdbs

def convert_pdb_to_mmtf(job, sfam_id, jobStoreIDs=None, clustered=True):
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
        pdb, chain, sdi, domain = zip(*[id_format.match(seq.id[:-2]).groups() \
            for s in SeqIO.parse(clusters_file, "fasta")])
    except ValueError:
        job.log("Unable to create data loading file for {}.".format(sfam_id))
        return

    domains = pd.DataFrame({"pdb":pdb, "chain":chain, "domNo":domain, "sdi":sdi})

    data_loader = os.path.join(pdb_path, "{}.h5".format(int(sfam_id)))
    domains.to_hdf(unicode(data_loader), "table", complevel=9, complib="bzip2")

def process_sfam(job, sfam_id, pdbFileStoreID, cores=1):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    sdoms_file = job.fileStore.readGlobalFile(pdbFileStoreID)

    sdoms = pd.read_hdf(unicode(sdoms_file), "merged")
    sdoms = sdoms[sdoms["sfam_id"]==sfam_id]["sdi"].drop_duplicates().dropna()

    if cores >= 20:
        setup_dask(cores)
        d_sdoms = dd.from_pandas(sdoms, npartitions=cores)
        processed_domains = d_sdoms.apply(lambda row: process_domain(job, row.sdi),
            axis=1).compute()
    else:
        processed_domains = job.addChildJobFn(map_job_rv, process_domain, sdoms,
            preemptable=True).rv()

    return processed_domains

def post_process_sfam(job, sfam_id, jobStoreIDs):
    cluster(job, sfam_id, jobStoreIDs)

    if False:
        convert_pdb_to_mmtf(job, sfam_id, jobStoreIDs)
        create_data_loader(job, sfam_id, jobStoreIDs)

def start_toil(job, name="prep_protein"):
    """Start the workflow to process PDB files"""
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    #Download PDB info
    sdoms_file = os.path.join(work_dir, "PDB.h5")
    in_store.read_input_file("PDB.h5", sdoms_file)

    #Add pdb info into local job store
    pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

    #Get all unique superfamilies
    sdoms = pd.read_hdf(unicode(sdoms_file), "merged")
    sfams = sdoms["sfam_id"].drop_duplicates().dropna()

    max_cores = job.fileStore.jobStore.config.maxCores if \
        job.fileStore.jobStore.config.maxCores > 2 else \
        job.fileStore.jobStore.config.defaultCores

    #Add jobs for each sdi
    job.addChildJobFn(map_job, process_sfam, sfams, pdbFileStoreID=pdbFileStoreID,
        cores=max_cores)

    #Add jobs for to post process each sfam
    job.addFollowOnJobFn(map_job, post_process_sfam, sfams, pdbFileStoreID=pdbFileStoreID,
        cores=max_cores)

    del sdoms
    os.remove(sdoms_file)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "never"

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)

    # meta = pd.DataFrame({c:[1] for c in cols})
    #
    # with ProgressBar():
    #     pdb_map = ddf.map_partitions(lambda _df: _df.apply(\
    #         lambda row: convert_row(job, row), axis=1), meta=meta)\
    #             .compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)

    # sdoms["sdi"].drop_duplicates().dropna().apply(
    #     lambda sdi: job.addChildJobFn(process_domain, dataset_name, sdi))
    # # job.log("STARTING {} SDIs".format(len(sdis)))
    # # for sdi in sdis: #.apply(lambda sdi:
    # #     #job.log(str(sdi))
    # #     job.addChildJobFn(process_domain, dataset_name, sdi)#)


    #
    # #Add jobs for to post process each sfam
    # sdoms["sfam_id"].drop_duplicates().dropna().apply(
    #     lambda sfam_id: job.addFollowOnJobFn(post_process_sfam, dataset_name, sfam_id))
    # job.log("STARTING {} sfams".format(len(sfams)))
    # for sfam_id in sfams: #.apply(lambda sfam_id: \
    #     #job.log(str(sfam_id))
    #     job.addFollowOnJobFn(post_process_sfam, dataset_name, sfam_id)#)


    #job.log("ADDED ALL FOLLOW JOBS {}".format(len(job._followOns)))
    # job.addFollowOnJobFn(cluster, dataset_name, sfam_id)
    # j2.addFollowOnJobFn(create_data_loader, dataset_name, sfam_id)
    # j2.addFollowOnJobFn(convert_pdb_to_mmtf, dataset_name, sfam_id)



# def toil_cdd(job, dataset_name, sfam_id):
#     pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)))
#
#     #Make CDD directory
#     if not os.path.isdir(pdb_path):
#         os.makedirs(pdb_path)
#
#     all_sdoms = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "merged")
#     sdoms_groups = all_sdoms.groupby("sfam_id")
#
#     try:
#         sdoms = sdoms_groups.get_group(sfam_id)
#     except KeyError:
#         job.log("No structural domains for {}".format(sfam_id))
#         return
#
#     #Make pdb group directories
#     sdoms["pdbId"].drop_duplicates().apply(lambda x: x[1:3]).drop_duplicates(). \
#         apply(lambda g: os.makedirs(os.path.join(pdb_path, g)) if not \
#             os.path.isdir(os.path.join(pdb_path, g)) else None)
#
#     #Some sdis have multiple entries for split structure, just use the one sdi
#     sdoms = sdoms[["pdbId", "chnLett", "sdi", "domNo"]].drop_duplicates()
#
#     #Process each domain
#     for pdbId, chnLett, sdi, domNo in sdoms.itertuples(index=False):
#         domain_pdb_path = os.path.join(pdb_path, pdbId[1:3].upper())
#         domain_file = os.path.join(domain_pdb_path, "{}_{}_sdi{}_d{}.pdb".format(pdbId, chnLett, sdi, domNo))
#         if not os.path.isfile(domain_file):
#             job.addChildJobFn(process_domain, dataset_name, pdbId, chnLett, sdi, domNo, sfam_id)
#             # if True:
#             #process_domain(job, dataset_name, pdbId, chnLett, sdi, domNo, sfam_id)

    # for cdd, sfam_id in iter_cdd(use_id=True, label="Ig"):
    #     toil_cdd(None, "default", cdd, sfam_id)
    #     cluster(None, "default", cdd)
    #     convert_pdb_to_mmtf(None, "default", cdd)


#
# def run_protein(pdb_file, chain=None, sdi=None, domainNum=None, cdd=None, process_chains=True, process_domains=True):
#     """Prepare a protein structure for use in molmimic. This will
#     0) Unzip if gzipped
#     1) Cleanup PDB file and add TER lines in between gaps
#     2) Remove HETATMS
#     3) Split Chains into separate files
#     4) Each chain is then protonated and minimized. See 'run_single_chain'
#     for more info
#
#     Parameters
#     ----------
#     pdb_file : str
#         Path to PDB file
#     chain : str or None
#         Chain ID to split out, protonate, and mininize. If None (default),
#         all chains will be split out, protonated, and mininized.
#
#     Returns
#     -------
#     If no chain was specified a list of all chains is returned with the
#     path to the prepared file, the pdb code, and chain. If a chain is
#     specified, 3 values are returned: the path to the prepared file,
#     the pdb code, and chain. See 'run_single_chain' for info.
#     """
#     print pdb_file
#     if not os.path.isfile(pdb_file):
#         raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))
#
#     base = os.path.basename(pdb_file)
#     if base.startswith("pdb") and base.endswith(".ent.gz"):
#         name_format = "^pdb([A-Za-z0-9]{4}).ent.gz"
#     else:
#         name_format = "^([A-Za-z0-9]{4}).pdb"
#
#     match = re.match(name_format, base)
#     if match and PDB_PATH is not None:
#         pdb = match.group(1)
#         if cdd is not None:
#             pdb_path = os.path.join(PDB_PATH, cdd, pdb[1:3].lower())
#         else:
#             pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
#         if not os.path.exists(pdb_path):
#             os.makedirs(pdb_path)
#     else:
#         print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
#         pdb_path = os.getcwd()
#         pdb = os.path.basename(os.path.splitext(pdb_file.replace(".gz", ""))[0])
#
#     #Unzip PDB file
#     if pdb_file.endswith(".gz"):
#         unzipped_pdb = ""
#         with gzip.open(pdb_file, 'rt') as f:
#             unzipped_pdb = f.read()
#
#     if chain is None:
#         print "run all chains"
#         #Split PDB into chains, 1 chain per file
#         if not pdb_file.endswith(".gz"):
#             subprocess.call([os.path.join(PDB_TOOLS, "pdb_splitchain.py"), pdb_file])
#         else:
#             splitchains = subprocess.Popen([os.path.join(PDB_TOOLS, "pdb_splitchain.py")])
#             splitchains.communicate(unzipped_pdb)
#
#         #Process all chains
#         if not process_chains:
#             for chain in get_all_chains(pdb_file):
#                 chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
#                 yield chain_file, pdb, chain, (None, None)
#         else:
#             for chain in get_all_chains(pdb_file):
#                 chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
#                 for domain_file, pdb_name, chain, (sdi, domainNum) in run_single_chain(chain_file):
#                     yield domain_file, pdb_name, chain, (sdi, domainNum)
#
#     else:
#         #Split desired chain in PDB into 1 file
#         chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
#         with open(chain_file, "w") as chainf:
#             if not pdb_file.endswith(".gz"):
#                 subprocess.call([os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain), pdb_file], stdout=chainf)
#             else:
#                 splitchain = subprocess.Popen([os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain)], stdin=subprocess.PIPE, stdout=chainf)
#                 splitchain.communicate(unzipped_pdb)
#
#         if not process_chains:
#             yield chain_file, pdb, chain, (None, None)
#         else:
#             for domain_file, pdb_name, chain, (sdi, domainNum) in run_single_chain(chain_file, domainNum=domainNum, sdi=sdi):
#                 yield domain_file, pdb_name, chain, (sdi, domainNum)
#
# def run_single_domain(chain_file, pdb, chain, chainNum=None, sdi=None, calculate_features=False):
