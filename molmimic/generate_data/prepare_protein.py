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

from util import data_path_prefix, get_structures_path, get_features_path, \
    get_first_chain, get_all_chains, number_of_lines, iter_unique_superfams

NUM_WORKERS = 20
dask.config.set(scheduler='multiprocessing', num_workers=NUM_WORKERS)
dask.config.set(pool=ThreadPool(NUM_WORKERS))

RAW_PDB_PATH = os.path.join(data_path_prefix, "pdb", "pdb")
PDB_PATH = os.path.join(data_path_prefix, "structures")
PDB_TOOLS = os.path.join(os.path.dirname(data_path_prefix), "pdb-tools")

def SubprocessChain(commands, output):
    if len(commands) > 2:
        prev_proc = subprocess.Popen(
            commands[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        for cmd in commands[1:-1]:
            proc = subprocess.Popen(
                cmd,
                stdin=prev_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            prev_proc = proc
        final_proc = subprocess.Popen(
            cmd[-1],
            stdin=prev_proc.stdout,
            stdout=output,
            stderr=subprocess.PIPE)
        return final_proc.communicate()
    elif len(commands) == 2:
        prev_proc = subprocess.Popen(
            commands[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        final_proc = subprocess.Popen(
            commands[1],
            stdin=prev_proc.stdout,
            stdout=output,
            stderr=subprocess.PIPE)
    elif len(commands) == 1:
        final_proc = subprocess.Popen(
            commands[0],
            stdout=output,
            stderr=subprocess.PIPE)
    else:
        raise RuntimeError


def extract_domain(dataset_name, pdb_file, pdb, chain, sdi, rslices, domNo, sfam_id, rename_chain=None, striphet=True):
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

    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)), pdb[1:3].upper())
    domain_file = os.path.join(pdb_path, "{}_{}_sdi{}_d{}.pdb.extracted".format(pdb, chain, int(sdi), domNo))

    if not os.path.exists(pdb_path):
        os.makedirs(pdb_path)

    open_fn = gzip.open if pdb_file.endswith(".gz") else open

    if pdb_file.endswith(".gz"):
        input = domain_file+".full"
        with gzip.open(pdb_file, 'rt') as zipf, open(input, "w") as pdb:
            pdb.write(zipf.read())
    else:
        input = pdb_file

    commands = [
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", input],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain)],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+rslices,
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_striphet.py")],
        [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
    ]

    if rename_chain is not None:
        commands.append([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"),
            "-{}".format("1" if isinstance(rename_chain, bool) and rename_chain else rename_chain)])

    with open(domain_file, "w") as output:
        SubprocessChain(commands, output)
        # splitmodel = subprocess.Popen(
        #     [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", input],
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE)
        # splitchain = subprocess.Popen(
        #     [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain)],
        #     stdin=splitmodel.stdout,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE)
        # delocc = subprocess.Popen(
        #     [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
        #     stdin=splitchain.stdout,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE)
        # splitdomain = subprocess.Popen(
        #     [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+rslices,
        #     stdin=delocc.stdout,
        #     stdout=subprocess.PIPE if rename_chain is not None else output,
        #     stderr=subprocess.PIPE)
        #
        # if rename_chain is not None:
        #     renamechain = subprocess.Popen(
        #         [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-1"],
        #         stdin=splitdomain.stdout,
        #         stdout=output,
        #         stderr=subprocess.PIPE)
        #     renamechain.communicate()
        # else:
        #     splitdomain.communicate()
        #
        # #Cleanup PDB file and add TER lines in between gaps
        # tidy = subprocess.Popen(
        #     [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), pdb_file],
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE)
        # #Remove HETATMS
        # striphet = subprocess.Popen(
        #     [sys.executable, os.path.join(PDB_TOOLS, "pdb_striphet.py")],
        #     stdin=tidy.stdout,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE)

    if pdb_file.endswith(".gz"):
        os.remove(domain_file+".full")

    return domain_file

def prepare_domain(pdb_file, pdb=None, chain=None, domainNum=None, sdi=None, sfam_id=None):
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

    ##FIME: Not sure why this none..
    num_chains = len(get_all_chains(pdb_file))
    if not len(get_all_chains(pdb_file)) == 1:
        raise RuntimeError("Must contain only one chain. There {} chains {}".format(num_chains, pdb_file))

    pdb_path = os.path.dirname(pdb_file)
    prefix = pdb_file.split(".", 1)[0]

    cleaned_file = prefix+".pdb.cleaned"
    # with open(cleaned_file, "w") as output:
    #     #Cleanup PDB file and add TER lines in between gaps
    #     tidy = subprocess.Popen(
    #         [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), pdb_file],
    #         stdout=subprocess.PIPE,
    #         stderr=subprocess.PIPE)
    #     #Remove HETATMS
    #     striphet = subprocess.Popen(
    #         [sys.executable, os.path.join(PDB_TOOLS, "pdb_striphet.py")],
    #         stdin=tidy.stdout,
    #         stdout=subprocess.PIPE,
    #         stderr=subprocess.PIPE)
    #     #Remove altLocs
    #     delloc = subprocess.Popen(
    #         [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
    #         stdin=striphet.stdout,
    #         stdout=output,
    #         stderr=subprocess.PIPE)
    #     delloc.communicate()

    #Add hydrogens -- FIXME biowulf and singilarity path hardcoded
    pqr_file = prefix+".pdb.pqr"
    propka_file = prefix+".propka"

    try:
        subprocess.call(["pdb2pqr",
            "--ff=parse",
            "--ph-calc-method=propka",
            "--chain",
            "--drop-water",
            pdb_file, pqr_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs. Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    try:
        with open(pqr_file) as f:
            pass
    except IOError:
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs. Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    #Minimize, assumes minimize from rosetta is in path
    score_file = prefix+".sc"
    minimized_file = prefix+".pdb_0001.pdb"

    try:
        subprocess.check_output(["minimize.static.linuxgccrelease",
            "-s", pqr_file,
            "-run:min_type", "lbfgs_armijo_nonmonotone",
            "-run:min_tolerance", "0.001",
            "-overwrite", "false", #Pandas apply calls first row twice so this is needed
            "-ignore_zero_occupancy", "false",
            "-out:file:scorefile", score_file,
            "-out:path:pdb", pdb_path,
            "-out:path:score", pdb_path],
            stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        raise RuntimeError("Unable to minimize file {}".format(pqr_file))

    attempts = 0
    while not os.path.isfile(minimized_file) or number_of_lines(minimized_file) == 0:
        if attempts >= 10:
            raise RuntimeError("Unable to minimize file {}".format(pqr_file))
        time.sleep(1)
        attempts += 1

    shutil.move(minimized_file, minimized_file+".rosetta")

    #Cleanup Rosetta PDB
    cleaned_minimized_file = prefix+".pdb.min"
    with open(cleaned_minimized_file, "w") as cleaned_minimized:
        subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_stripheader.py"), minimized_file+".rosetta"], stdout=cleaned_minimized)

    if len(get_all_chains(cleaned_minimized_file)) > 1:
        #Sometimes the chains are split in pdb2pqr, fixes it to one chain
        one_chain_clean_file = prefix+".min.pdb.one_chain"
        with open(one_chain_clean_file, "w") as one_chain_clean:
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(chain), cleaned_minimized_file], stdout=one_chain_clean)
        shutil.move(cleaned_minimized_file, cleaned_minimized_file+".multi_chain")

        cleaned_file = prefix+".pdb"
        tidyed_pdb_file = prefix+".min.pdb.tidy"
        with open(tidyed_pdb_file, "w") as cf:
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), one_chain_clean_file], stdout=cf)

    else:
        cleaned_file = prefix+".pdb"
        tidyed_pdb_file = prefix+".min.pdb.tidy"
        with open(tidyed_pdb_file, "w") as cf:
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), cleaned_minimized_file], stdout=cf)

    shutil.move(tidyed_pdb_file, cleaned_file)

    attempts = 0
    while number_of_lines(cleaned_file) == 0:
    	if attempts >= 10:
    		raise RuntimeError("Invalid PDB file")
        time.sleep(0.2)
        attempts += 1

    return cleaned_file

def process_domain(job, dataset_name, sdi, pdb=None, chain=None, domNo=None, sfam_id=None):
    job.log("SYS EXE = {}".format(sys.executable))
    all_sdoms = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "merged")

    sdom = all_sdoms[all_sdoms["sdi"]==float(sdi)]
    if sdom.shape[0] == 0:
        job.log("SDI {} does not exist".format(sdi))
        return

    sfam_id = sdom.iloc[0].sfam_id
    pdb = sdom.iloc[0].pdbId
    chain = sdom.iloc[0].chnLett
    domNo = sdom.iloc[0].domNo

    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(sfam_id), pdb[1:3].upper())
    domain_file = os.path.join(pdb_path, "{}_{}_sdi{}_d{}.pdb".format(pdb, chain, sdi, domNo))
    if os.path.isfile(domain_file):
        return

    pdb_file = os.path.join(RAW_PDB_PATH, pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    print pdb_file, pdb, chain, sdi, domNo

    if not os.path.isfile(pdb_file):
        job.log("Cannot process {}.{}.d{} ({}), pdb DNE".format(pdb, chain, domNo, sdi))
        return
    try:
        #Multiple ranges for one sdi. Might be domain swaped?
        rslices = ["{}:{}".format(st, en) for st, en in sdom[["from", "to"]].drop_duplicates().itertuples(index=False)]
        domain_file = extract_domain(dataset_name, pdb_file, pdb, chain, sdi, rslices, domNo, sfam_id)
        try:
            prepare_domain(domain_file)
            print "Wrote", domain_file
        except RuntimeError as e:
            job.log(str(e))
            pass
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        job.log("Cannot process {}.{}.d{} ({}), error: {}".format(pdb, chain, domNo, sdi, e))

def convert_pdb_to_mmtf(job, dataset_name, sfam_id):
    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)))
    chemcomp = os.path.join(data_path_prefix, "pdb", "chemcomp")
    if not os.path.isdir(chemcomp):
        os.makedirs(chemcomp)

    spark_env = os.environ.copy()
    spark_env["PDB_DIR"] = chemcomp
    spark_env["PDB_CACHE_DIR"] = chemcomp

    mmtf_path = os.path.join(data_path_prefix, "mmtf", dataset_name, "by_superfamily", str(int(sfam_id)))
    if not os.path.isdir(mmtf_path):
        os.makedirs(mmtf_path)

    #Convert PDBs to MMTF-Hadoop Sequence file directory
    try:
        subprocess.call(["spark-submit",
            "--class", "edu.sdsc.mmtf.spark.io.demos.PdbToMmtfFull",
            "{}/target/mmtf-spark-0.2.0-SNAPSHOT.jar".format(os.environ["MMTFSPARK_HOME"]),
            pdb_path, mmtf_path],
            env=spark_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        job.log("Error converting to MMTF: {}".format(e))
        pass

def cluster(job, dataset_name, sfam_id, id=0.95):
    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)))

    if not os.path.isdir(pdb_path):
        return

    #Save all domains to fasta
    domain_fasta = os.path.join(pdb_path, "{}.fasta".format(int(sfam_id)))
    with open(domain_fasta, "w") as fasta:
        for f in glob.glob(os.path.join(pdb_path, "*", "*.pdb")):
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_toseq.py"),
                f], stdout=fasta, stderr=subprocess.PIPE)

    #Order domains by resolution so the ones with the highest resolutions are centroids
    resolutions = pd.read_table(os.path.join(os.path.dirname(RAW_PDB_PATH), "resolu.idx"),
        header=None, names=["pdbId", "resolution"], skiprows=6, sep="\t;\t")

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

    #Cluster using uclust
    uclust_file = os.path.join(pdb_path, "{}_clusters.uc".format(int(sfam_id)))
    clusters_file = os.path.join(pdb_path, "{}_nr.fasta".format(int(sfam_id)))
    subprocess.call(["usearch", "-cluster_fast", domain_fasta, "-id", str(id),
        "-centroids", clusters_file, "-uc", uclust_file])

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
    uclust.to_hdf(unicode(os.path.join(pdb_path, "{}_clusters.h5".format(int(sfam_id)))), "table", complevel=9, complib="bzip2")
    os.remove(uclust_file)

    #Copy clustered to new directory
    cluster_path = os.path.join(pdb_path, "clustered")
    if not os.path.isdir(cluster_path):
        os.makedirs(cluster_path)

    for seq in SeqIO.parse(clusters_file, "fasta"):
        domainId = seq.id[:-2]
        pdb_group = domainId[1:3].upper()
        out_dir = os.path.join(cluster_path, pdb_group)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        fname = "{}.pdb".format(domainId)
        shutil.copy(os.path.join(pdb_path, pdb_group, fname), os.path.join(out_dir, fname))

def create_data_loader(job, dataset_name, sfam_id):
    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily", str(int(sfam_id)))
    clusters_file = os.path.join(pdb_path, "{}_nr.fasta".format(int(sfam_id)))
    id_format = re.compile("^([A-Z0-9]{4})_([A-Za-z0-9]+)_sdi([0-9]+)_d([0-9]+)$")

    try:
        pdb, chain, sdi, domain = zip(*[id_format.match(seq.id[:-2]).groups() \
            for s in SeqIO.parse(clusters_file, "fasta")])
    except ValueError:
        job.log("Unable to create data loading file for {}.".format(sfam_id))
        return

    domains = pd.DataFrame({"pdb":pdb, "chain":chain, "domNo":domain, "sdi":sdi})

    data_loader = os.path.join(pdb_path, "{}.h5".format(int(sfam_id)))
    domains.to_hdf(unicode(data_loader), "table", complevel=9, complib="bzip2")

def process_sfam(job, dataset_name, sfam_id, cores=NUM_WORKERS):
    sdoms = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "merged")
    sdoms = sdoms[sdoms["sfam_id"]==sfam_id]
    sdoms = sdoms[["pdbId", "chnLett", "sdi", "domNo"]].drop_duplicates().dropna()
    d_sdoms = dd.from_pandas(sdoms, npartitions=cores)
    d_sdoms.apply(lambda row: process_domain(job, dataset_name, row.sdi) \
        if not os.path.isfile(os.path.join(
            PDB_PATH, dataset_name, "by_superfamily", str(sfam_id),
            row.pdbId[1:3].upper(), "{}_{}_sdi{}_d{}.pdb".format(
                row.pdbId, row.chnLett, row.sdi, row.domNo
            )
        )) else None, axis=1).compute()

    # sdoms.apply(lambda row: job.addChildJobFn(process_domain, dataset_name, row.sdi) \
    #     if not os.path.isfile(os.path.join(
    #         PDB_PATH, dataset_name, "by_superfamily", str(sfam_id),
    #         row.pdbId[1:3].upper(), "{}_{}_sdi{}_d{}.pdb".format(
    #             row.pdbId, row.chnLett, row.sdi, row.domNo
    #         )
    #     )) else None, axis=1)

def post_process_sfam(job, dataset_name, sfam_id):
    cluster(job, dataset_name, sfam_id)
    #create_data_loader(job, dataset_name, sfam_id)


def start_toil(job, dataset_name, name="prep_protein", memory="120000M"):
    pdb_path = os.path.join(PDB_PATH, dataset_name, "by_superfamily")
    sdoms = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "merged")


    #Make sfam directories
    # sfams = dd.from_pandas(sdoms["sfam_id"].drop_duplicates().dropna(), npartitions=NUM_WORKERS)
    # sfams.apply(lambda g: os.makedirs(os.path.join(pdb_path, str(int(g)))) \
    #     if not os.path.isdir(os.path.join(pdb_path, str(int(g)))) else None).compute()
    #
    # #Make pdb group direcotries
    # sfam_pdbs = sdoms[["sfam_id","pdbId"]].drop_duplicates().dropna()
    # sfam_pdbs["pdbId"] = sfam_pdbs["pdbId"].apply(lambda s: s[1:3])
    # sfam_pdbs =  dd.from_pandas(sfam_pdbs.drop_duplicates(), npartitions=NUM_WORKERS)
    # sfam_pdbs.apply(lambda d: os.makedirs(
    #         os.path.join(pdb_path, str(int(d["sfam_id"])), d["pdbId"])
    #     ) if not os.path.isdir(
    #         os.path.join(pdb_path, str(int(d["sfam_id"])), d["pdbId"])) else None,
    #     axis=1).compute()
    # job.log("DONE DIRS")

    if dataset_name != "parallel":
        #Add jobs for each sdi
        sdoms["sfam_id"].drop_duplicates().dropna().apply(
            lambda sfam_id: job.addChildJobFn(process_sfam, dataset_name, sfam_id))

        print "ADDED ALL JOBS {}".format(len(job._children))

        #Add jobs for to post process each sfam
        sdoms["sfam_id"].drop_duplicates().dropna().apply(
            lambda sfam_id: job.addFollowOnJobFn(post_process_sfam, dataset_name, sfam_id))
    else:
        def addDomain(s, j, d, c):
            job.log("ADDING domain {} {} {}".format(s, j, d))
            return j.addChildJobFn(process_domain, d, s)
            #djob = JobFunctionWrappingJob(process_domain, f, s)
            #djob._addPredecessor(j)
            #c.append(djob)

        def addSfam(s, j, d, c):
            job.log("ADDING sfam {} {} {}".format(s, j, d))
            return j.addFollowOnJobFn(post_process_sfam, d, s)
            #sjob = JobFunctionWrappingJob(post_process_sfam, d, s)
            #sjob._addPredecessor(j)
            #c.append(sjob)

        job.log("STARTING ADD DOMAINS")
        Parallel(n_jobs=NUM_WORKERS)(delayed(addDomain)(s, job, dataset) for s in \
            sdoms["sdi"].drop_duplicates().dropna())
        job.log("ADDED DOMAINS: {}".format(len(job._children)))

        job.log("STARTING ADD SFAMS")
        Parallel(n_jobs=NUM_WORKERS)(delayed(addSfam)(s,job, dataset) for s in \
            sdoms["sfam_id"].drop_duplicates().dropna())
        job.log("STARTING ADD SFAMS")

        # job.log("STARTING ADD DOMAINS")
        # dd.from_pandas(sdoms["sdi"].drop_duplicates().dropna(),
        #     npartitions=NUM_WORKERS).apply(addDomain, args=(job, dataset_name)).compute().tolist()
        # job.log("ADDED DOMAINS: {}".format(len(job._children)))
        #
        # job.log("STARTING ADD SFAMS")
        # dd.from_pandas(sdoms["sfam_id"].drop_duplicates().dropna(),
        #     npartitions=NUM_WORKERS).apply(addSfam, args=(job, dataset_name)).compute().tolist()
        # job.log("ADDED SFAMS: {}".format(len(job._followOns)))


        # del sdis
        # print "ADDED ALL JOBS {}".format(len(job._children))
        #
        # #Add jobs for each sdi
        # sfams = dd.from_pandas(sdoms["sfam_id"].drop_duplicates().dropna(), npartitions=NUM_WORKERS)
        # job._followOns += sfams.apply(sfamJob, meta=('job', 'object')).compute().tolist()

    # def domainJob(sdi):
    #     djob = JobFunctionWrappingJob(process_domain, dataset_name, sdi)
    #     djob._addPredecessor(job)
    #     return djob
    #
    # def sfamJob(sfam):
    #     sjob = JobFunctionWrappingJob(post_process_sfam, dataset_name, sfam)
    #     sjob._addPredecessor(job)
    #     return sjob

    #Add jobs for each sdi
    # sdis = dd.from_pandas(sdoms["sdi"].drop_duplicates().dropna(), npartitions=NUM_WORKERS)
    # print "STARTING ADD DOMAINS"
    # job._children += sdis.apply(domainJob, meta=('job', 'object')).compute().tolist()
    # del sdis
    # print "ADDED ALL JOBS {}".format(len(job._children))
    #
    # #Add jobs for each sdi
    # sfams = dd.from_pandas(sdoms["sfam_id"].drop_duplicates().dropna(), npartitions=NUM_WORKERS)
    # job._followOns += sfams.apply(sfamJob, meta=('job', 'object')).compute().tolist()
    # del sfams

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


    job.log("ADDED ALL FOLLOW JOBS {}".format(len(job._followOns)))
    # job.addFollowOnJobFn(cluster, dataset_name, sfam_id)
    # j2.addFollowOnJobFn(create_data_loader, dataset_name, sfam_id)
    # j2.addFollowOnJobFn(convert_pdb_to_mmtf, dataset_name, sfam_id)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "never"
    options.leaderCores = 20
    dataset_name = options.jobStore.split(":")[-1]

    job = Job.wrapJobFn(start_toil, dataset_name)
    with Toil(options) as toil:
        if options.restart:
            toil.restart(job)
        else:
            toil.start(job)

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


# if __name__ == "__main__":
#     if len(sys.argv) == 3:
#     	#Submit IBIS interfaces file to get PDBs to prepare to batch system
#         submit_ibis(sys.argv[2], sys.argv[1])
#
#     elif len(sys.argv) == 4 and sys.argv[1] == "load-pdb-group":
#         PDB_PATH = get_structures_path(sys.argv[2])
#         load_pdb_group(sys.argv[2], sys.argv[3])
#
#     elif len(sys.argv) == 5 and "load" in sys.argv:
#     	#Load IBIS interfaces file to get PDBs to prepare
#         PDB_PATH = get_structures_path(sys.argv[2])
#         load_ibis(sys.argv[2], sys.argv[3], sys.argv[4])
#
#     elif len(sys.argv) == 4 and "load-families" in sys.argv:
#     	#Load IBIS interfaces file to get PDBs to prepare from in indivual protein family
#         PDB_PATH = get_structures_path(sys.argv[2])
#         submit_ibis_cdd(sys.argv[2], sys.argv[3])
#
#     elif len(sys.argv) in [4,5] and "protein" in sys.argv:
#     	#Prepare a single protein by splitting in chains, then calls run_single_chain for each chain
#     	PDB_PATH = get_structures_path(sys.argv[2])
#         chain = sys.argv[3] if len(sys.argv)==4 else None
#         for c in run_protein(sys.argv[3], chain):
#             pass
#
#     elif len(sys.argv) == 5 and "single-chain" in sys.argv:
#     	#Prepare a single chain
#     	PDB_PATH = get_structures_path(sys.argv[3], cdd=sys.argv[2])
#         for d in run_single_chain(sys.argv[3]):
#             pass
#
#     elif len(sys.argv) in [6, 8] and "single-domain" in sys.argv:
#     	#Split on domains after chain has been prepared
#     	PDB_PATH = get_structures_path(sys.argv[2])
#         for d in run_single_domain(*sys.argv[3:]):
#             pass
#
#     else:
#         print len(sys.argv), sys.argv
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
#     assert not chain_file.endswith(".gz"), "Cannot be a gzip archive, try 'run_protein' instead"
#     assert len(get_all_chains(chain_file)) == 1, "Must contain only one chain"
#
#     for domain_file, pdb_name, chain, (sdi, domainNum) in split_domains_in_pdb(chain_file, pdb, chain):
#         if calculate_features:
#             subprocess.call(["/data/draizene/3dcnn-torch-py2",
#                 "python", os.path.realpath(__file__), pdb, chain, "None", "False"])
#         yield domain_file, pdb_name, chain, (sdi, domainNum)
#
# def load_ibis(dataset_name, cdd, ibis_data, process_all_chains=False, add_to_job=None, wait_for_subjobs=False):
#     from molmimic.torch_model.torch_loader import IBISDataset
#
#     output_dir = get_structures_path(dataset_name)
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#
#     #chain_names_file = open(os.path.join(output_dir, "{}_raw_chains.txt".format(cdd.replace("/", ""))), "w")
#
#     try:
#         dataset = IBISDataset(ibis_data)
#     except:
#         chain_names_file.close()
#         return
#
#     data = dataset.data
#
#     name = os.path.splitext(os.path.basename(ibis_data))[0]
#
#     if isinstance(add_to_job, SwarmJob):
#         job = add_to_job
#     elif add_to_job==True:
#         job_name = os.environ.get("SLURM_JOB_NAME", "{}_prep_chains".format(dataset_name))
#     	job = SwarmJob("{}_prep_chains".format(job_name), modules="rosetta", walltime="8:00:00")
#     else:
#         job = None
#
#     i = 0
#
#     pdb_groups = data.groupby(lambda x: data["pdb"].loc[x][1:3])
#     for pdb_divided, pdb_group in pdb_groups:
#         pdb_path = os.path.join(output_dir, pdb_divided.lower())
#         feature_path = os.path.join(output_dir, "atom", pdb_divided.lower())
#         if not os.path.exists(pdb_path):
#             os.makedirs(pdb_path)
#         if not os.path.exists(feature_path):
#             os.makedirs(feature_path)
#
#         if not process_all_chains:
#             for _, pdb in pdb_group[["pdb", "chain"]].drop_duplicates().iterrows():
#                 pdb_file = "/pdb/pdb/{}/pdb{}.ent.gz".format(pdb_divided.lower(), pdb["pdb"].lower())
#                 try:
#                     for chain_file, pdb, chain, domain in run_protein(pdb_file, pdb["chain"], cdd=cdd, process_chains=False):
#                         if job is not None:
#                             job += "python {} single-chain {} {}\n".format(__file__, dataset_name, chain_file)
#                         #print >> chain_names_file, chain_file
#                 except RuntimeError:
#                     continue
#         else:
#             for pdb in pdb_group["pdb"].drop_duplicates():
#                 pdb_file = "/pdb/pdb/{}/pdb{}.ent.gz".format(pdb_divided.lower(), pdb.lower())
#                 try:
#                     for chain_file, pdb, chain, domain in run_protein(pdb_file, cdd=cdd, process_chains=False):
#                         if job is not None:
#                             job += "python {} single-chain {} {}\n".format(__file__, dataset_name, chain_file)
#                         #print >> chain_names_file, chain_file
#                 except RuntimeError:
#                     continue
#
#     #chain_names_file.close()
#
#     if isinstance(add_to_job, bool) and add_to_job:
#         print "=> Chain Job", job.run()
#
#     #return chain_names_file.name
#
# def all_chains_done(chain_names_file):
#     if isinstance(chain_names_file, (list, tuple)):
#         for f in chain_names_file:
#             if not all_chains_done(f):
#                 return False
#         return True
#
#     with open(chain_names_file) as f:
#         for line in f:
#             if not os.path.isfile(line.rstrip()+".done"):
#                 return False
#     return True
#
# def load_families(dataset_name, ibis_data, run_chain_swarm=True):
#     job_name = "{}_prep_chains".format(dataset_name)
#     if run_chain_swarm:
#         job_name = os.environ.get("SLURM_JOB_NAME", job_name)
#
#     CDD = pd.read_csv("MMDB/StructDomSfam.csv", usecols=["label"]).drop_duplicates().dropna()
#     CDD = sorted(CDD["label"].apply(lambda cdd: cdd.replace("/", "").replace("'", "\'")).tolist())
#
#     chain_names_files = []
#     job = SwarmJob(job_name, modules="rosetta", walltime="2:00:00")
#     for cdd in CDD:
#         f = os.path.join(ibis_data, "{}.raw".format(cdd))
#         chain_names = load_ibis(dataset_name, cdd, f, add_to_job=job)
#         chain_names_files.append(chain_names)
#
#     if run_chain_swarm:
#         job.run()
#
#         while not all_chains_done(chain_names_files):
#             #Wait for chains to finish
#             time.sleep(800)
#
# def load_pdb_group(dataset_name, group):
#     output_dir = get_structures_path(dataset_name)
#     with open(os.path.join("prepare_chains", group+".txt")) as f:
#         for line in f:
#             try:
#                 pdb, chain = line.rstrip().split("\t")
#             except ValueError:
#                 continue
#             chain_file = os.path.join(output_dir, group.lower(), "{}_{}.pdb".format(pdb.lower(), chain))
#             run_single_chain(chain_file)
#
# def submit_ibis_cdd(dataset_name, ibis_data, job_name="prep_proteins", dependency=None):
#     output_dir = get_structures_path(dataset_name)
#     CDD = pd.read_csv("MMDB/StructDomSfam.csv", usecols=["label"]).drop_duplicates().dropna()
#     CDD = sorted(CDD["label"].apply(lambda cdd: cdd.replace("/", "").replace("'", "\'")).tolist())
#
#     pdbs = defaultdict(set)
#
#     #full_job = SwarmJob(job_name+"_full", walltime="0:30:00")
#
#     for cdd in CDD:
#         cdd_f = os.path.join(ibis_data, "{}.raw".format(cdd.replace("/", "").replace("'", "\'")))
#         #full_job += "python {} load {} \"{}\" {}\n".format(__file__, dataset_name, cdd.replace("/", ""), cdd_f)
#
#         with open(cdd_f) as f:
#             for line in f:
#                 try:
#                     pdb, chain = line.split("\t")[:2]
#                 except IndexError:
#                     continue
#                 pdbs[pdb[1:3]].add((pdb, chain))
#
#
#     #jid = full_job.run(filter_unique=True, dependency=dependency)
#     #print jid
#
#     chain_job = SwarmJob(job_name+"_chain", modules="rosetta", walltime="2-00:00:00")
#     if not os.path.isdir("prepare_chains"):
#         os.makedirs("prepare_chains")
#     for pdb_group, pdb_chains in pdbs.iteritems():
#         with open(os.path.join("prepare_chains", pdb_group.lower()+".txt"), "w") as f:
#             for pdb, chain in pdb_chains:
#                 print >> f, "{}\t{}".format(pdb, chain)
#         chain_job += "python {} load-pdb-group {} {}\n".format(__file__, dataset_name, pdb_group)
#     jid1 = chain_job.run(filter_unique=True)#, dependency="afterany:"+jid)
#     return jid1
#
# def submit_ibis(dataset_name, ibis_data, job_name="prep_proteins", dependency=None, run_chain_swarm=False):
#     command = "load-families" if os.path.isdir(ibis_data) else "load"
#     job = SwarmJob(job_name, walltime="96:00:00", modules="rosetta", mem="1", individual=True)
#     job += "python {} {} {} {}\n".format(__file__, command, dataset_name, ibis_data)
#     job_id = job.submit_individual(dependency=dependency)
#     swarm_file = "{}_prep_chains.sh".format(dataset_name) if not run_chain_swarm else "{}.sh".format(dataset_name)
#     return job_id, swarm_file
#
# def submit_chain_swarm(swarm_file, num_tasks, job_name="prep_proteins", dependency=None):
#     job = SwarmJob(job_name, walltime="96:00:00", mem="1", modules="rosetta")
#     job += "awk 'NR==$SLURM_ARRAY_TASK_ID{print;exit}' {}".format(swarm_file)
#     job.run(filter_unique=True)
