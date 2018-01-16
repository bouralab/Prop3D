import sys
sys.path.append("/data/draizene/3DUnetCNN")
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

    def __init__(self, name, cpus=1, mem="60g", walltime="72:00:00", gpus=False, user_parameters=None):
        parameters = []

        self.name = name
        self.cpus = cpus
        self.cmd_file = "{}.sh".format(self.name)
        self.job_id = None
        self.cmd = "#!/bin/sh\n"

        self.parameters = [
            "--job-name={}".format(name),
            "--cpus-per-task={}".format(cpus),
            "--mem={}".format(mem),
            "--time={}".format(walltime)
        ]

        if gpus:
            self.parametersa += [
                "--partition=gpu",
                "--gres=gpu:k80:{}".format(gpus if isinstance(gpus, int) and gpus > 0 else 4)
                ]

        if user_parameters is not None and isinstance(user_parameters, (list, tuple)):
            self.parameters += user_parameters

        for param in self.parameters:
            self.cmd += "#SBATCH {}\n".format(param)

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

    def write(self):
        with open(self.cmd_file, "w") as f:
            f.write(self.cmd)

    def submit(self, write=True, hold_jid=None):
        if write:
            self.write()

        while not SwarmJob.can_add_job():
            time.sleep(0.1)

        cmd = ["/usr/local/slurm/bin/sbatch"]

        if hold_jid:
            cmd.append("--dependency=afterany:{}".format(",".join(hold_jid)))

        cmd.append(self.cmd_file)
        
        self.job_id = subprocess.check_output(cmd)
        time.sleep(0.1)
        return self.job_id

def calculate_features(pdb, chain, resi, id):
    from molmimic.biopdbtools import Structure
    print Structure.features_from_string(pdb, chain, resi, id=id, force_feature_calculation=True)

def load_ibis(ibis_data):
    from molmimic.keras_model.pdb_generator import IBISGenerator
    data = IBISGenerator(ibis_data, input_shape=(96,96,96,59))
    for i, row in data.data.iterrows():
        print "Running {} ({}.{}): {}".format(row["unique_obs_int"], row["pdb"], row["chain"], ",".join(["{}{}".format(i,n) for i, n in zip(row["resi"].split(","), row["resn"].split(","))]))
        job = SwarmJob(row["unique_obs_int"])
        job += "/data/draizene/3dcnn python {} {} {} {} {}\n".format(os.path.realpath(__file__), row["pdb"], row["chain"], row["resi"], row["unique_obs_int"])
        job.submit()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        load_ibis(sys.argv[1])
    elif len(sys.argv) == 5:
        calculate_features(*sys.argv[1:])
    else:
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain, resi, id)")
