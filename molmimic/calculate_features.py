import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/usr/share/pdb2pqr")

import os
import argparse
import time
import math
from datetime import datetime
from itertools import groupby
import subprocess

class SwarmJob(object):
    start_date = datetime.today().strftime('%Y-%m-%d')
    max_jobs = 1000
    max_cpus = 56

    def __init__(self, name, cpus=1, mem="60", walltime="72:00:00", gpus=None, user_parameters=None, modules=None, individual=False, merge_output=False, threads_per_job=None):
        parameters = []

        self.name = name
        self.cpus = cpus
        self.mem = mem
        self.modules = modules
        self.merge_output = merge_output
        self.cmd_file = "{}.sh".format(self.name)
        self.job_id = None
        self.walltime = walltime
        self.cmd = "#!/bin/sh\n"

        if self.modules is not None and individual:
            for module in self.modules.split(","):
                self.cmd += "module load {}\n".format(module)

        self.cpu_parameters = "--cpus-per-task={}".format(cpus)

        self.parameters = [
            "--job-name={}".format(name),
            "--mem={}g".format(mem),
            "--time={}".format(walltime)
        ]
        self.parameters += [self.cpu_parameters]

        if isinstance(gpus, int):
            self.gpus = gpus
            self.gpu_parameters = [
                "--partition=gpu",
                "--gres=gpu:k80:{}".format(self.gpus)
                ]
            self.parameters += self.gpu_parameters
        elif gpus is not None and gpus==True:
            self.gpus = math.ceil(self.cpus/14.)
            self.gpu_parameters = [
                "--partition=gpu",
                "--gres=gpu:k80:{}".format(gpus if isinstance(gpus, int) and gpus > 0 else 4)
                ]
            self.parameters += self.gpu_parameters
        else:
            self.gpus = 0
            self.gpu_parameters = []

        self.threads_per_job = threads_per_job

        if user_parameters is not None and isinstance(user_parameters, (list, tuple)):
            self.parameters += user_parameters

    def __iadd__(self, new):
        self.cmd += new
        return self

    @staticmethod
    def slurm_enabled():
        try:
            subprocess.check_output(["sacct"])
            return True
        except:
            return False

    @staticmethod
    def number_of_jobs_running():
        num_lines = len(subprocess.check_output(["sacct", "--state", "R", "--starttime", SwarmJob.start_date]).split("\n"))
        return num_lines-2 if num_lines > 2 else 0

    @staticmethod
    def number_of_jobs_pending():
        num_lines = len(subprocess.check_output(["/usr/local/slurm/bin/sacct", "--state", "PENDING", "--starttime", SwarmJob.start_date]).split("\n"))
        return num_lines-2 if num_lines > 2 else 0

    @staticmethod
    def can_add_job():
        return SwarmJob.number_of_jobs_running()+SwarmJob.number_of_jobs_pending() < SwarmJob.max_jobs

    def add_individual_parameters(self):
        for param in self.parameters:
            self.cmd += "#SBATCH {}\n".format(param)

    def write(self):
        with open(self.cmd_file, "w") as f:
            f.write(self.cmd)

    def submit_individual(self, write=True, hold_jid=None):
        if write:
            self.write()

        while not SwarmJob.can_add_job():
            time.sleep(0.3)

        cmd = ["/usr/local/slurm/bin/sbatch"]

        if hold_jid:
            cmd.append("--dependency=afterany:{}".format(",".join(hold_jid)))

        cmd.append(self.cmd_file)

        for _ in xrange(5):
            try:
                self.job_id = subprocess.check_output(cmd)
                break
            except subprocess.CalledProcessError:
                time.sleep(0.3)

        time.sleep(0.3)
        return self.job_id

    def run(self, write=True):
        if write:
            self.write()

        cmd = ["swarm", "--file", self.cmd_file,
               "--logdir", "{}_logs".format(self.name),
               "-g", str(self.mem),
               "--job-name", self.name
               ]
        if self.walltime is not None:
            cmd += ["--time", self.walltime]

        if self.cpus > 1:
            cmd += ["--sbatch", "'{}'".format(self.cpu_parameters)]

        if self.gpus >= 1:
            cmd += self.gpu_parameters

        if self.threads_per_job:
            cmd += ["-t", str(self.threads_per_job)]

        if self.modules is not None:
            cmd += ["--module", self.modules]

        if self.merge_output:
            cmd += ["--merge-output"]

        for _ in xrange(5):
            try:
                self.job_id = subprocess.check_output(cmd)
                break
            except subprocess.CalledProcessError:
                time.sleep(0.3)

        time.sleep(0.3)
        return self.job_id

def calculate_features(pdb, chain, resi, course_grained):
    from molmimic.biopdbtools import Structure
    resi = None if resi == "None" else resi
    course_grained = course_grained == "True"
    course_grained = bool(course_grained)
    Structure.features_from_string(pdb, chain, resi=resi, course_grained=course_grained, force_feature_calculation=True, grid=False)

def load_ibis(ibis_data, course_grained=False, cellular_organisms=False, unclustered=False):
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data, input_shape=(512,512,512), cellular_organisms=cellular_organisms)
    data = dataset.data #if course_grained else dataset.full_data
    job = SwarmJob("ibis_features", walltime="06:00:00")
    values = [tuple(x) for x in data[["pdb", "chain"]].values]

    key = lambda p: p[0][1:3]
    pdbs = sorted(values, key=key)

    i = 0
    for pdb_divided, pdb_group in groupby(pdbs, key=key):
        #id = dataset.full_data.loc[(dataset.full_data["pdb"]==row["pdb"])&(dataset.full_data["chain"]==row["chain"])].iloc[0]["gi"]
        #if os.path.isfile("/data/draizene/molmimic/features3/atom/{}/{}_{}.npy".format(row["pdb"], row["chain"])):
        #    continue
        feature_path = "/data/draizene/molmimic/features_Ig/atom/{}".format(pdb_divided.lower())
        if not os.path.exists(feature_path):
            os.makedirs(feature_path)
        for pdb, chain in pdb_group:
            i += 1
            print "Running {}/{}: {}.{}".format(i, data.shape[0], pdb, chain)
            job += "/data/draizene/3dcnn-torch-py2 python {} {} {} None {} \n".format(os.path.realpath(__file__), pdb, chain, course_grained)
    print job.run()


if __name__ == "__main__":
    if len(sys.argv) in [2, 3]:
        course_grained = "--course-grained" in sys.argv
        cellular_organisms = "--cellular-organisms" in sys.argv
        load_ibis(sys.argv[-1], course_grained, cellular_organisms)
    elif len(sys.argv) == 5:
        calculate_features(*sys.argv[1:])
    else:
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain, resi, id, use_course_grained)")
