import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/usr/share/pdb2pqr")

import os
import argparse
import time
from datetime import datetime
import subprocess

class SwarmJob(object):
    start_date = datetime.today().strftime('%Y-%m-%d')
    max_jobs = 1000
    max_cpus = 56

    def __init__(self, name, cpus=1, mem="60", walltime="72:00:00", gpus=False, user_parameters=None):
        parameters = []

        self.name = name
        self.cpus = cpus
        self.mem = mem
        self.cmd_file = "{}.sh".format(self.name)
        self.job_id = None
        self.cmd = "#!/bin/sh\n"

        self.parameters = [
            "--job-name={}".format(name),
            "--cpus-per-task={}".format(cpus),
            "--mem={}g".format(mem),
            "--time={}".format(walltime)
        ]

        if gpus:
            self.parameters += [
                "--partition=gpu",
                "--gres=gpu:k80:{}".format(gpus if isinstance(gpus, int) and gpus > 0 else 4)
                ]

        if user_parameters is not None and isinstance(user_parameters, (list, tuple)):
            self.parameters += user_parameters


    def __iadd__(self, new):
        self.cmd += new
        return self

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
               "--logdir", "logs",
               "-g", str(self.mem),
               ]

        for _ in xrange(5):
            try:
                self.job_id = subprocess.check_output(cmd)
                break
            except subprocess.CalledProcessError:
                time.sleep(0.3)

        time.sleep(0.3)
        return self.job_id

def calculate_features(pdb, chain, resi, id, course_grained):
    from molmimic.biopdbtools import Structure
    resi = None if resi == "None" else resi
    course_grained = course_grained == "True"
    course_grained = bool(course_grained)
    Structure.features_from_string(pdb, chain, resi=resi, id=id, course_grained=course_grained, force_feature_calculation=True, grid=False)

def load_ibis(ibis_data, course_grained=False):
    from molmimic.torch_model.torch_loader import IBISDataset
    print "Loading"
    dataset = IBISDataset(ibis_data, input_shape=(512,512,512))
    print "Loaded"
    data = dataset.data #if course_grained else dataset.full_data
    parsing = True
    job = SwarmJob("ibis_features")
    for i, row in data.iterrows():
        #if row["pdb"]=="3OXQ" and row["chain"] == "A":
        #    parsing = True
        #    continue
        if not parsing:
            continue

        id = dataset.full_data.loc[(dataset.full_data["pdb"]==row["pdb"])&(dataset.full_data["chain"]==row["chain"])].iloc[0]["gi"]
        print "Running {}: {}.{}".format(id, row["pdb"], row["chain"])
        resi = None
        # else:
        #     print "Running {} ({}.{}): {}".format(row["unique_obs_int"], row["pdb"], row["chain"], ",".join(["{}{}".format(i,n) for i, n in zip(row["resi"].split(","), row["resn"].split(","))]))
        #     id = row["unique_obs_int"]
        #     resi = row["resi"]

        job += "/data/draizene/3dcnn-torch-py2 python {} {} {} {} {} {}\n".format(os.path.realpath(__file__), row["pdb"], row["chain"], resi, id, course_grained)
    job.run()

if __name__ == "__main__":
    if len(sys.argv) in [2, 3]:
        if len(sys.argv) == 3:
            if sys.argv[1] == "--course-grained":
                load_ibis(sys.argv[2], True)
            elif sys.argv[2] == "--course-grained":
                load_ibis(sys.argv[1], True)
        else:
            load_ibis(sys.argv[1])
    elif len(sys.argv) == 6:
        calculate_features(*sys.argv[1:])
    else:
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain, resi, id, use_course_grained)")
