import os, sys
sys.path.append("/data/draizene/pdb-tools")
sys.path.append("/data/draizene/molmimic")

import subprocess
import shutil
import gzip
import glob
import re
import time
from itertools import groupby

from molmimic.calculate_features import SwarmJob
from molmimic.split_domains import split_domains_in_pdb
from util import get_structures_path, get_features_path, get_first_chain, get_all_chains, number_of_lines

PDB_PATH = None

#os.environ.get("MOLMIMIC_PDB_PATH", os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdbs"))

def run_protein(pdb_file, chain=None, sdi=None, domainNum=None):
    """Prepare a protein structure for use in molmimic. This will
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
    print pdb_file
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    base = os.path.basename(pdb_file)
    if base.startswith("pdb") and base.endswith(".ent.gz"):
        name_format = "^pdb([A-Za-z0-9]{4}).ent.gz"
    else:
        name_format = "^([A-Za-z0-9]{4}).pdb"

    match = re.match(name_format, base)
    if match and PDB_PATH is not None:
        pdb = match.group(1)
        pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
    else:
        print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
        pdb_path = os.getcwd()
        pdb = os.path.basename(os.path.splitext(pdb_file.replace(".gz", ""))[0])

    #Unzip PDB file
    if pdb_file.endswith(".gz"):
        unzipped_pdb = ""
        with gzip.open(pdb_file, 'rt') as f:
            unzipped_pdb = f.read()

    if chain is None:
        print "run all chains"
        #Split PDB into chains, 1 chain per file
        if not pdb_file.endswith(".gz"):
            subprocess.call(["/data/draizene/pdb-tools/pdb_splitchain.py", pdb_file])
        else:
            splitchains = subprocess.Popen(["/data/draizene/pdb-tools/pdb_splitchain.py"])
            splitchains.communicate(unzipped_pdb)

        #Process all chains
        if SwarmJob.slurm_enabled():
            for chain in get_all_chains(pdb_file):
                chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
                yield chain_file, pdb, chain, (None, None)
        else:
            for chain in get_all_chains(pdb_file):
                chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
                for domain_file, pdb_name, chain, (sdi, domainNum) in run_single_chain(chain_file):
                    yield domain_file, pdb_name, chain, (sdi, domainNum)

    else:
        #Split desired chain in PDB into 1 file
        chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
        with open(chain_file, "w") as chainf:
            if not pdb_file.endswith(".gz"):
                subprocess.call(["/data/draizene/pdb-tools/pdb_selchain.py", "-{}".format(chain), pdb_file], stdout=chainf)
            else:
                splitchain = subprocess.Popen(["/data/draizene/pdb-tools/pdb_selchain.py", "-{}".format(chain)], stdin=subprocess.PIPE, stdout=chainf)
                splitchain.communicate(unzipped_pdb)

        if SwarmJob.slurm_enabled():
            yield chain_file, pdb, chain, (None, None)
        else:
            for domain_file, pdb_name, chain, (sdi, domainNum) in run_single_chain(chain_file, domainNum=domainNum, sdi=sdi):
                yield domain_file, pdb_name, chain, (sdi, domainNum)

def run_single_chain(pdb_file, domainNum=None, sdi=None, split_domains=False, calculate_features=False):
    """Prepare a single chain for use in molmimic (called during 'run_protein').
    This prepares the chains by:
    (1) Removing all altLoc's except the first
    (2) Add hydrogens using pdb2pqr (ff=parse, ph=propka)
    (3) Minimzed usinf rosetta (lbfgs_armijo_nonmonotone with tolerance 0.001)
    (4) Cleaned so that can be parsed by simple PDB parsers

    Note: raises an AssertError if pdb_file is gzipped or contains more than
    one chain

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
    print "Running", pdb_file
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    assert not pdb_file.endswith(".gz"), "Cannot be a gzip archive, try 'run_protein' instead"
    assert len(get_all_chains(pdb_file)) == 1, "Must contain only one chain"

    file_prefix = os.path.splitext(pdb_file)[0]
    base = os.path.basename(pdb_file)
    name_format = "^([A-Za-z0-9]{4})_([A-Za-z0-9]{1,3}).pdb$"
    match = re.match(name_format, base)
    if match and PDB_PATH is not None:
        pdb, chain = match.group(1), match.group(2)
        pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
    else:
        print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
        pdb_path = os.getcwd()
        pdb = os.path.basename(file_prefix)
        chain = get_first_chain(pdb_file)

    #Cleanup PDB file and add TER lines in between gaps
    tidyed_pdb_file = os.path.join(pdb_path, "{}_{}.pdb.tidy".format(pdb, chain))
    with open(tidyed_pdb_file, "w") as tidyf:
        subprocess.call(["/data/draizene/pdb-tools/pdb_tidy.py", pdb_file], stdout=tidyf)

    #Remove HETATMS
    no_het_file = os.path.join(pdb_path, "{}_{}.pdb.nohet".format(pdb, chain))
    with open(no_het_file, "w") as no_het_f:
        subprocess.call(["/data/draizene/pdb-tools/pdb_striphet.py", tidyed_pdb_file], stdout=no_het_f)

    #Remove altLocs
    delocc_file = os.path.join(pdb_path, "{}.delocc".format(pdb_file))
    with open(delocc_file, "w") as delocc:
        subprocess.call(["/data/draizene/pdb-tools/pdb_delocc.py", no_het_file], stdout=delocc)

    #Add hydrogens
    pqr_file = os.path.join(pdb_path, "{}.pqr.pdb".format(file_prefix))
    propka_file = os.path.join(pdb_path, "{}.propka".format(pqr_file))
    print "Run PDB2PQR", " ".join(["/data/draizene/3dcnn-torch-py2", "shell", "pdb2pqr", "--ff=parse", "--ph-calc-method=propka", "--chain", "--drop-water", delocc_file, pqr_file])

    subprocess.call(["/data/draizene/3dcnn-torch-py2", "shell", "pdb2pqr", "--ff=parse", "--ph-calc-method=propka", "--chain", "--drop-water", delocc_file, pqr_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print "Finished running PDB2PQR"
    try:
        #rename propKa file
        shutil.move(propka_file, "{}.propka".format(file_prefix))
    except IOError:
        pass

    try:
        with open(pqr_file) as f:
            pass
    except IOError:
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs. Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    #Minimize
    score_file = "{}.sc".format(file_prefix)
    minimized_file = "{}.pqr_0001.pdb".format(file_prefix)

    subprocess.call(["minimize",
        "-s", pqr_file,
        "-run:min_type", "lbfgs_armijo_nonmonotone",
        "-run:min_tolerance", "0.001",
        "-ignore_zero_occupancy", "false",
        "-out:file:scorefile", score_file,
        "-out:path:pdb", pdb_path,
        "-out:path:score", pdb_path])

    attempts = 0
    while not os.path.isfile(minimized_file) or number_of_lines(minimized_file) == 0:
    	if attempts >= 10:
    		raise RuntimeError("Unable to minimize file {}".format(pqr_file))
        time.sleep(0.2)
        attempts += 1

    #Cleanup Rosetta PDB
    cleaned_minimized_file = "{}.min.pdb".format(file_prefix)
    with open(cleaned_minimized_file, "w") as cleaned_minimized:
        subprocess.call(["/data/draizene/pdb-tools/pdb_stripheader.py", minimized_file], stdout=cleaned_minimized)

    if len(get_all_chains(cleaned_minimized_file)) > 1:
        #Sometimes the chains are split in pdb2pqr, fixes it to one chain
        one_chain_clean_file = "{}.min.one_chain.pdb".format(file_prefix)
        with open(one_chain_clean_file, "w") as one_chain_clean:
            subprocess.call(["/data/draizene/pdb-tools/pdb_chain.py", "-{}".format(chain), cleaned_minimized_file], stdout=one_chain_clean)

        cleaned_file = "{}.min.pdb".format(file_prefix)
        tidyed_pdb_file = "{}.min.pdb.tidy".format(file_prefix)
        with open(tidyed_pdb_file, "w") as cf:
            subprocess.call(["/data/draizene/pdb-tools/pdb_tidy.py", one_chain_clean_file], stdout=cf)

    else:
        cleaned_file = "{}.min.pdb".format(file_prefix)
        tidyed_pdb_file = "{}.min.pdb.tidy".format(file_prefix)
        with open(tidyed_pdb_file, "w") as cf:
            subprocess.call(["/data/draizene/pdb-tools/pdb_tidy.py", cleaned_minimized_file], stdout=cf)

    os.remove(cleaned_file)
    shutil.move(tidyed_pdb_file, cleaned_file)

    print "num lines cleaned", number_of_lines(cleaned_file)
    attempts = 0
    while number_of_lines(cleaned_file) == 0:
    	if attempts >= 10:
    		raise RuntimeError("Invalud PDB file")
        print "num lines cleaned", number_of_lines(cleaned_file)
        time.sleep(0.2)
        attempts += 1

    return_file = cleaned_file

    for domain_file, pdb_name, chain, (sdi, domainNum) in split_domains_in_pdb(return_file, pdb, chain):
        if domainNum is not None and sdi is not None:
            yield domain_file, pdb_name, chain, (sdi, domainNum)

    if domainNum is None and sdi is None:
        yield return_file, pdb, chain, (None, None)

def run_single_domain(chain_file, pdb, chain, chainNum=None, sdi=None, calculate_features=False):
    assert not chain_file.endswith(".gz"), "Cannot be a gzip archive, try 'run_protein' instead"
    assert len(get_all_chains(chain_file)) == 1, "Must contain only one chain"

    for domain_file, pdb_name, chain, (sdi, domainNum) in split_domains_in_pdb(chain_file, pdb, chain):
        if calculate_features:
            subprocess.call(["/data/draizene/3dcnn-torch-py2",
                "python", os.path.realpath(__file__), pdb, chain, "None", "False"])
        yield domain_file, pdb_name, chain, (sdi, domainNum)

def load_ibis(ibis_data, dataset_name, process_all_chains=False, add_to_job=None):
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data)
    data = dataset.data

    key = lambda p: p[1:3]
    pdbs = sorted(set(data["pdb"].values.tolist()), key=key)

    name = os.path.splitext(os.path.basename(ibis_data))[0]

    if isinstance(add_to_job, SwarmJob):
        job = add_to_job
    elif add_to_job==True:
        job_name = os.environ.get("SLURM_JOB_NAME", "{}_prep_chains".format(dataset_name))
    	job = SwarmJob("{}_prep_chains".format(job_name), modules="rosetta", walltime="8:00:00")
    else:
        job = None


    i = 0
    pdb_groups = data.groupby(lambda x: data["pdb"].loc[x][1:3])
    for pdb_divided, pdb_group in pdb_groups:
        pdb_path = os.path.join(get_structures_path(dataset_name), pdb_divided.lower())
        feature_path = os.path.join(get_features_path(dataset_name), "atom", pdb_divided.lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
        if not os.path.exists(feature_path):
            os.makedirs(feature_path)

        if not process_all_chains:
            for _, pdb in pdb_group[["pdb", "chain"]].drop_duplicates().iterrows():
                pdb_file = "/pdb/pdb/{}/pdb{}.ent.gz".format(pdb_divided.lower(), pdb["pdb"].lower())
                try:
                    for chain_file, pdb, chain, domain in run_protein(pdb_file, pdb["chain"]):
                        if job is not None:
                            job += "python {} single-chain {} {}\n".format(__file__, dataset_name, chain_file)
                except RuntimeError:
                    continue
        else:
            for pdb in pdb_group["pdb"].drop_duplicates():
                pdb_file = "/pdb/pdb/{}/pdb{}.ent.gz".format(pdb_divided.lower(), pdb.lower())
                try:
                    for chain_file, pdb, chain, domain in run_protein(pdb_file):
                        if job is not None:
                            job += "python {} single-chain {} {}\n".format(__file__, dataset_name, chain_file)
                except RuntimeError:
                    continue

    if isinstance(add_to_job, bool) and add_to_job:
    	job.run()

def load_families(ibis_data, dataset_name, run_chain_swarm=True):
    job_name = "{}_prep_chains".format(dataset_name)
    if run_chain_swarm:
        job_name = os.environ.get("SLURM_JOB_NAME", job_name)
    job = SwarmJob(job_name, modules="rosetta", walltime="2:00:00")
    for f in glob.glob(os.path.join(ibis_data, "*.tsv")):
        load_ibis(f, dataset_name, job=job)
    if run_chain_swarm:
        job.run(filter_unique=True, update_dependencies=True)

def submit_ibis(ibis_data, dataset, job_name="prep_proteins", dependency=None, run_chain_swarm=False):
    command = "load-families" if os.path.isdir(ibis_data) else "load"
    job = SwarmJob(job_name, walltime="96:00:00", modules="rosetta", individual=True)
    job += "python {} {} {} {}\n".format(__file__, command, dataset, ibis_data)
    job_id = job.submit_individual(dependency=dependency)
    swarm_file = "{}_prep_chains.sh".format(dataset) if not run_chain_swarm else "{}.sh".format(dataset)
    return job_id, swarm_file

def submit_chain_swarm(swarm_file, num_tasks, job_name="prep_proteins", dependency=None):
    job = SwarmJob(job_name, walltime="96:00:00", modules="rosetta")
    job += "awk 'NR==$SLURM_ARRAY_TASK_ID{print;exit}' {}".format(swarm_file)
    job.run(filter_unique=True)

if __name__ == "__main__":
    if len(sys.argv) == 3:
    	#Submit IBIS interfaces file to get PDBs to prepare to batch system
        submit_ibis(sys.argv[2], sys.argv[1])

    elif len(sys.argv) == 4 and "load" in sys.argv:
    	#Load IBIS interfaces file to get PDBs to prepare
        PDB_PATH = get_structures_path(sys.argv[2])
        load_ibis(sys.argv[3], sys.argv[2])

    elif len(sys.argv) == 4 and "load-families" in sys.argv:
    	#Load IBIS interfaces file to get PDBs to prepare from in indivual protein family
        PDB_PATH = get_structures_path(sys.argv[2])
        load_families(sys.argv[3], sys.argv[2], process_all_chains=False)

    elif len(sys.argv) in [4,5] and "protein" in sys.argv:
    	#Prepare a single protein by splitting in chains, then calls run_single_chain for each chain
    	PDB_PATH = get_structures_path(sys.argv[2])
        chain = sys.argv[3] if len(sys.argv)==4 else None
        for c in run_protein(sys.argv[3], chain):
            pass

    elif len(sys.argv) == 4 and "single-chain" in sys.argv:
    	#Prepare a single chain
    	PDB_PATH = get_structures_path(sys.argv[2])
        for d in run_single_chain(sys.argv[3]):
            pass

    elif len(sys.argv) in [6, 8] and "single-domain" in sys.argv:
    	#Split on domains after chain has been prepared
    	PDB_PATH = get_structures_path(sys.argv[2])
        for d in run_single_domain(*sys.argv[3:]):
            pass

    else:
        print len(sys.argv), sys.argv
