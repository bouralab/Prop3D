import os, sys
sys.path.append("/data/draizene/molmimic")

import shutil
import glob

from molmimic.calculate_features import SwarmJob
from molmimic.util import get_features_path, get_structures_path

def filter_ibis(dataset_name, ibis_data, cleanup=False):
    from molmimic.torch_model.torch_loader import IBISDataset

    outpath = "{}.tsv".format(os.path.splitext(ibis_data)[0])
    with open(ibis_data) as f, open(outpath, "w") as new:
        for i, line in enumerate(f):
            if i==0:
                new.write(line)
            else:
                try:
                    pdb, chain, sdi, domNum, _ = line.split("\t")
                except ValueError:
                    continue
                pdb_id = "{}_{}_sdi{}_{}".format(pdb, chain, sdi, domNum)
                features_path = os.path.join(get_features_path(dataset_name), "atom", pdb[1:3].lower())
                structure_path = os.path.join(get_structures_path(dataset_name), pdb[1:3].lower())
                
                if os.path.isfile(os.path.join(features_path, "{}.npy".format(pdb_id))) and \
                  os.path.isfile(os.path.join(structure_path, "{}.pdb".format(pdb_id))):
                    new.write(line)

    if cleanup:
        os.remove(ibis_data)
        shutil.move(outpath, ibis_data)

def run_filter(dataset_name, ibis_data, job_name="filter_ibis", dependency=None):
    if os.path.isdir(ibis_data):
        ibis_data_files = glob.glob(os.path.join(ibis_data, "*.tsv"))
    else:
        ibis_data_files = [ibis_data]
    
    job = SwarmJob(job_name, mem="1", individual=len(ibis_data_files)==1)
    for ibis_data in ibis_data_files:
        job += "python {} filter {}\n".format(__file__, dataset_name, ibis_data)

    if len(ibis_data_files)==1:
        jid = job.submit_individual(dependency=dependency)
    else:
        jid = job.run(dependency=dependency)

    features_path = os.path.join(get_features_path(dataset_name), "done")
    endjob = SwarmJob(job_name, mem="1", individual=True)
    endjob += "touch {}\n".format(features_path)
    endjob.submit_individual(dependency="afterany:"+jid)

    while not os.path.isfile(features_path):
        time.sleep(800)

def submit_filter(dataset_name, ibis_data, job_name="filter_ibis", dependency=None):
    job = SwarmJob(job_name, walltime="2-00:00:00", mem="1", individual=True)
    job += "python {} {} {}\n".format(__file__, dataset_name, ibis_data)
    job_id = job.submit_individual(dependency=dependency)
    return job_id

if __name__ == "__main__":
    if len(sys.argv) == 3:
        run_filter(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4 and sys.argv[1] == "filter":
        filter_ibis(sys.argv[2], sys.argv[3])