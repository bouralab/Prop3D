import os, sys
sys.path.append("/data/draizene/molmimic")

import requests

from molmimic import calculate_ibis_dataset
from molmimic import split_domains
from molmimic import prepare_protein
from molmimic import calculate_features
from molmimic import filter_ibis_dataset

def build_ibis_dataset(dataset_name):
    """Build dataset from IBIS to be used in molmimic. This script will submit
    jobs to a slurm controller to:
        1) Donwload raw ibis data if not present
        2) Parse the raw ibis data into qury domains
        3) Add domain information into ibis data split by domain
        4) Prepare all proteins in IBIS
        4) Calculatess all possible features for proteins in IBIS
        5) Removes interfaces that failed prepartion and feature calculation

    Note: Jobs are all submitted with the same name. Some jobs spawn many
    subjobs using swarm that will also have the same name. All succsive jobs
    from this function will wait for all previous jobs (and subjobs) that have
    same name using the 'dependency=singleton' option.
    """
    interface_dir = os.path.join(os.path.dirname(__file__), "..", "data", "interfaces")

    #1) Get ibis
    ibisdown_path = os.path.join(interface_dir, "ibisdown")
    print "1. Getting IBIS dataset"
    if not os.path.isdir(ibisdown_path):
        print "=> Downloading ibisdown from ftp://ftp.ncbi.nih.gov/pub/mmdb/ibis/ibisdown.tar.gz..."
        #Downlaod from NCBI
        import requests, zipfile, StringIO
        r = requests.get("ftp://ftp.ncbi.nih.gov/pub/mmdb/ibis/ibisdown.tar.gz", stream=True)
        z = zipfile.ZipFile(StringIO.StringIO(r.content))
        z.extractall(os.path.join(interface_dir))
    else:
        print "=> Exists in molmimic data directory"

    #2) Parse IBIS down, split into files based on query domain
    print "2. Submitting jobs to parse IBIS dataset by CDD Domains..."
    ibis_by_cdd = os.path.join(interface_dir, dataset_name)
    jid2 = calculate_ibis_dataset.submit_ibis(dataset_name, ibisdown_path, job_name="build_ibis")
    print "=>Job ID:", jid2

    #3) Add domain information
    print "3. Submitting jobs to add domain information..."
    jid3 = split_domains.submit_cdd(ibis_by_cdd, job_name="build_ibis", dependency="afterany:"+jid2)
    print "=>Job ID:", jid3

    #4) Prepare all proteins
    print "4. Submitting jobs to prepare all proteins in IBIS dataset (protonation with pdb2pqr and minimization with rosetta)..."
    jid4, chain_swarm_file = prepare_protein.submit_ibis(dataset_name, ibis_by_cdd, job_name="build_ibis", dependency="afterany:"+jid3)
    print "=>Job ID:", jid4

    #5) Calculate features
    print "5. Submitting jobs to calculate all features for all proteins in IBIS dataset..."
    jid5 = calculate_features.submit_ibis(ibis_by_cdd, dataset_name, job_name="build_ibis", dependency="afterany:"+jid4)
    print "=>Job ID:", jid5

    #6) Remove entries that have no structures and features
    print "6. Filter interfaces to keep those with prepared structures and features..."
    jid6 = filter_ibis_dataset.submit_filter(dataset_name, ibis_by_cdd, job_name="build_ibis", dependency="afterany:"+jid5)
    print "=>Job ID:", jid6

if __name__ == "__main__":
    assert len(sys.argv) in [2,3], "{} Usage: python {} dataset_name [ibisdown_path]".format(sys.argv, sys.argv[0])
    build_ibis_dataset(*sys.argv[1:])
