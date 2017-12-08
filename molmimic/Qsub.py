import sys
import time
import subprocess

class Qsub:
    max_jobs = 500

    def __init__(self, name, threads=8):
        self.name = name
        self.threads = threads
        self.cmd_file = "{}.sh".format(self.name)
        self.job_id = None
        self.cmd = """#!/bin/sh
# Job Name
#$ -N {name}

# Only email if aborted
#$ -m a

# Execute the script from the Current Working Directo
#$ -cwd

# Tell the job your memory requirements
#$ -l h_rt=72:00:00,mem_free=80G,h_vmem=80G,reserve_mem=80G
{cores}

source /panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/6CysDB_env/bin/activate
""".format(
    name=self.name,
    cores="#$ -pe multicore {}".format(self.threads) if threads > 1 else "")

    def __iadd__(self, new):
        self.cmd += new
        return self

    @staticmethod
    def number_of_jobs_running():
        num_lines = len(subprocess.check_output(["qstat"]).split("\n"))
        return num_lines-2 if num_lines > 2 else 0

    @staticmethod
    def can_add_job():
        return Qsub.number_of_jobs_running() < Qsub.max_jobs

    def write(self):
        with open(self.cmd_file, "w") as f:
            f.write(self.cmd)

    def submit(self, write=True, hold_jid=None):
        if write:
            self.write()

        while not Qsub.can_add_job():
            time.sleep(0.1)

        cmd = ["qsub"]
        if hold_jid:
            cmd += ["-hold_jid", ",".join(hold_jid)]
        cmd.append(self.cmd_file)
        o = subprocess.check_output(cmd)
        print o
        time.sleep(0.1)
        self.job_id = o.split()[2]
        return self.job_id
