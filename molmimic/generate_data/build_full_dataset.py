import os, sys
sys.path.append("/data/draizene/molmimic")

import requests

from toil.common import Toil
from toil.job import Job

import download_data
import convert_mmdb_to_pdb
import get_structural_interactome
import calculate_bsa
import prepare_protein
import calculate_features
import filter_dataset

data_path_prefix = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

# from molmimic import calculate_ibis_dataset
# from molmimic import split_domains
# from molmimic import prepare_protein
# from molmimic import calculate_features
# from molmimic import filter_ibis_dataset

def build_ibis_dataset(dataset_name, force=False):
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

    #1) Get raw data
        #Download MMDB, IBIS observed, and IBIS inferred from Dropbox...
        #Download rsync snapshots of PDB (PDB format) SIFTS

    #2) Convert MMDB mapping to pdb
    jid1 = convert_mmdb_to_pdb.submit_jobs()

    #3) Calculate observed interactome
    jid2 = get_structural_interactome.submit_jobs_observed(dataset_name, dependency="afterok:"+jid1)

    #4) Calculate inferred interactome
    jid3 = get_structural_interactome.submit_jobs_inferred(dataset_name, dependency="afterok:"+jid2)

    #5) Calculate BSA for observed and inferred
    jid4 = calculate_bsa.submit_ibis_cdd(dataset_name, dependency="afterok:"+jid3)

    #6) Split PDB chains -- all PDBs from each family

    #7) Prepare all chains

    #8) Convert chains to MMTF

    #9) Calculate features

    #4) Prepare all proteins
    # print "4. Submitting jobs to prepare all proteins in IBIS dataset (protonation with pdb2pqr and minimization with rosetta)..."
    #
    # jid4, chain_swarm_file = prepare_protein.submit_ibis(dataset_name, ibis_by_cdd, run_chain_swarm=True, job_name="build_ibis")#, dependency="afterany:"+jid3)
    # print "=>Job ID:", jid4
    #
    # #5) Calculate features
    # print "5. Submitting jobs to calculate all features for all proteins in IBIS dataset..."
    # jid5 = calculate_features.submit_ibis(ibis_by_cdd, dataset_name, job_name="build_ibis", dependency="afterany:"+jid4)
    # print "=>Job ID:", jid5
    #
    # #6) Remove entries that have no structures and features
    # print "6. Filter interfaces to keep those with prepared structures and features..."
    # jid6 = filter_ibis_dataset.submit_filter(dataset_name, ibis_by_cdd, job_name="build_ibis", dependency="afterany:"+jid5)
    # print "=>Job ID:", jid6

def start_toil(dataset_name, options, use_data=False):
    if use_data:
        data = Job.wrapJobFn(download_data.start_toil).encapsulate()
        mmdb2pdb = data.addFollowOnJobFn(convert_mmdb_to_pdb.start_toil).encapsulate()
    else:
        mmdb2pdb = Job.wrapJobFn(convert_mmdb_to_pdb.start_toil).encapsulate()

    interactome = mmdb2pdb.addChildJobFn(get_structural_interactome.start_toil, dataset_name).encapsulate()
    bsa = interactome.addFollowOnJobFn(calculate_bsa.start_toil, dataset_name).encapsulate()

    prep_protein = mmdb2pdb.addChildJobFn(prepare_protein.start_toil, dataset_name).encapsulate()
    features = mmdb2pdb.addFollowOnJobFn(calculate_features.start_toil, dataset_name, name="features").encapsulate()

    filter = mmdb2pdb.addFollowOnJobFn(filter_dataset.start_toil, dataset_name, name="filter").encapsulate()

    with Toil(options) as toil:
        toil.start(mmdb2pdb if not use_data else data)

if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    dataset_name = options.jobStore.split(":")[-1]

    start_toil(dataset_name, options)
