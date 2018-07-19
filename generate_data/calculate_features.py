import sys
import os
import argparse
import time
from itertools import groupby
import glob

from util import get_structures_path, get_features_path, initialize_data

def calculate_features(dataset, pdb, chain, sdi=None, domain=None):
    from molmimic.biopdbtools import Structure
    print pdb, chain
    #Set paths correctly so Structure can find PDBs and Features
    initialize_data(dataset)
    s = Structure(pdb, chain, sdi=sdi, domain=domain, force_feature_calculation=True)
    s.get_flat_features()

def load_ibis(ibis_data, dataset_name, use_domains=0):
    from molmimic.torch_model.torch_loader import IBISDataset

    if os.path.isdir(ibis_data):
        ibis_data_files = glob.glob(os.path.join(ibis_data, "*.tsv"))
    else:
        ibis_data_files = [ibis_data]

    job_name = os.environ.get("SLURM_JOB_NAME", "{}_features".format(dataset_name))
    job = SwarmJob(job_name, walltime="06:00:00")

    features_path = get_features_path(dataset_name)
    all_features_file = os.path.join(features_path, "all_features.txt")
    with open(all_features_file, "w") as all_features:
        for ibis_data in ibis_data_files:
            dataset = IBISDataset(ibis_data, input_shape=(512,512,512))
            data = dataset.full_data #if course_grained else dataset.full_data
            i = 0
            pdb_groups = data.groupby(lambda x: data["pdb"].loc[x][1:3])
            for pdb_divided, pdb_group in pdb_groups:
                features_path = os.path.join(get_features_path(dataset_name), "atom")

                if not os.path.exists(features_path):
                    os.makedirs(features_path)
                for _, pdb in pdb_group.iterrows():
                    i += 1
                    if use_domains in [1, 2]:
                        print "Running {}/{}: {}.{}.{} ({})".format(i, data.shape[0], pdb["pdb"], pdb["chain"], pdb["domnum"], pdb["sdi"])
                        job += "/data/draizene/3dcnn-torch-py2 python {} features {} {} {} {} {}\n".format(os.path.realpath(__file__), dataset_name, pdb["pdb"], pdb["chain"], pdb["sdi"], pdb["domnum"])
                        print >> all_features, os.path.join(features_path, pdb["pdb"][1:3].lower(), "{}_{}_sdi{}_d{}.npy".format(pdb["pdb"], pdb["chain"], pdb["sdi"], pdb["domnum"]))
                    if use_domains in [0, 2]:
                        print "Running {}/{}: {}.{}".format(i, data.shape[0], pdb["pdb"], pdb["chain"])
                        job += "/data/draizene/3dcnn-torch-py2 python {} features {} {} \n".format(os.path.realpath(__file__), dataset_name, pdb["pdb"], pdb["chain"])
                        print >> all_features, os.path.join(features_path, pdb["pdb"][1:3].lower(), "{}_{}.npy".format(pdb["pdb"], pdb["chain"]))
            print
    jid = job.run()

    endjob = SwarmJob(job_name, individual=True)
    endjob += "touch {}.done\n".format(all_features_file)
    endjob.submit_individual(dependency="afterany:"+jid)

    while not os.path.isfile(all_features_file+".done"):
        time.sleep(800)

def submit_ibis(ibis_data, dataset_name, use_domains=2, job_name="build_ibis", dependency=None):
    domains = ["", " --only-domains ", " --domains "]
    job = SwarmJob(job_name, walltime="96:00:00", mem="1", individual=True)
    job += "python {} run{}{} {}\n".format(__file__, domains[use_domains], dataset_name, ibis_data)
    return job.submit_individual(dependency=dependency)

def toil_cdd(job, dataset_name, cdd, pdb_group):
    pdb_dir = os.path.join(get_structures_path(dataset_name), cdd, pdb_group)
    fname_re = re.compile("(\w{4})_(\w{1,2})_d(\d+)_sdi(\d+).min.pdb$")
    for f in glob.glob(os.path.join(pdb_dir, "*.min.pdb")):
        m = fname_re.match(os.path.basename(f))
        if not match:
            continue
        pdb, chain, domNo, sdi = m.groups()
        calculate_features(f, pdb, chain, domNo, sdi)

def start_toil(job, dataset_name):
    pdb_dir = get_structures_path(dataset_name)
    for cdd in iter_cdd():
        cdd = cdd.replace("/", "")
        for pdb_group in next(os.walk(os.path.join(pdb_dir, cdd)))[1]:
            job.addChildJobFn(toil_cdd, dataset_name, cdd, pdb_group)

if __name__ == "__main__":
    if len(sys.argv) in [3, 4]:
        if "--domains" in sys.argv:
            use_domains = 2
        elif "--only-domains" in sys.argv:
            use_domains = 1
        else:
            use_domains = 0

        submit_ibis(sys.argv[-1], sys.argv[-2], use_domains)
    elif len(sys.argv) in [4, 5] and sys.argv[1] == "run":
        if "--domains" in sys.argv:
            use_domains = 2
        elif "--only-domains" in sys.argv:
            use_domains = 1

        load_ibis(sys.argv[-1], sys.argv[-2], use_domains)
    elif len(sys.argv) == 4 and sys.argv[1]=="features":
        import re
        f = os.path.splitext(os.path.basename(sys.argv[-1]))
        m = re.match("^(.{4})_(.{1,2})_sdi(\d+)_d(\d+)\.pdb$", f)
        if m:
            calculate_features(sys.argv[2], m.groups(1), m.groups(2), sdi=m.groups(3), domain=m.groups(4))
        else:
            raise RuntimeError("PDB file must have the format: {PDB}_{Chain}_sdi{sdi}.d{domNo}.pdb")

    elif len(sys.argv) in [5, 7] and sys.argv[1]=="features":
        calculate_features(*sys.argv[2:])
    else:
        print len(sys.argv), sys.argv
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain, resi, id, use_course_grained)")
