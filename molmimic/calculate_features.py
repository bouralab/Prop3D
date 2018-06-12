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
import glob
import shutil

from molmimic.util import get_features_path, initialize_data

try:
    subprocess.call(["which", "sacct"])
    SLURM_AVAIL = True
except CalledProcessError:
    SLURM_AVAIL = False

try:
    subprocess.call(["which", "swarm"])
    SWARM_AVAIL = True
except CalledProcessError:
    SWARM_AVAIL = False

class SwarmJob(object):
    start_date = datetime.today().strftime('%Y-%m-%d')
    max_jobs = 1000
    max_cpus = 56

    def __init__(self, name, cpus=1, mem="60", walltime="72:00:00", gpus=None, user_parameters=None, modules=None, individual=False, merge_output=False, threads_per_job=None, cmd_file=None):
        parameters = []

        self.name = name
        self.cpus = cpus
        self.mem = mem
        self.modules = modules
        self.merge_output = merge_output
        self.cmd_file_name = cmd_file or "{}.sh".format(self.name)
        self.cmd_file = open(self.cmd_file_name, "w" if cmd_file is None else "a+")
        self.job_id = None
        self.walltime = walltime
        self += "#!/bin/sh\n"

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
        self.user_parameters = user_parameters

        if user_parameters is not None and isinstance(user_parameters, (list, tuple)):
            self.parameters += user_parameters

        if individual and cmd_file is None:
            self.add_individual_parameters()

    def __iadd__(self, new):
        self.cmd_file.write(new)
        return self

    @staticmethod
    def slurm_enabled():
        return SWARM_AVAIL

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

    @staticmethod
    def update_dependencies(job_name, old_job_id=None, new_job_id=None):
        old_job_id = old_job_id or os.environ.get("SLURM_JOB_ID", None)
        job_info = subprocess.check_output(["squeue", "-u", os.environ.get("SLURM_JOB_ACCOUNT", os.environ.get("USER", None)), "-o", "%A %j %V %E"])
        job_info = sorted([line.rstrip("\n").split() for line in job_info.split("\n")[1:]], key=lambda x: x[0] if isinstance(x, (list, tuple)) and len(x)>0 else x)
        set_wait_for_new_job = False
        for _job_name, job_id, starttime, dependency in job_info:
            if _job_name == job_name and _job_id != old_job_id:
                if dependency == "singelton" and not set_wait_for_new_job:
                    subprocess.call(["scontrol", "update", "JobId={}".format(_job_id), "Dependency=afterany:{}".format(new_job_id)])
                    set_wait_for_new_job = False
                elif dependency == "singelton":
                    subprocess.call(["scontrol", "update", "JobId={}".format(_job_id), "Dependency=singleton"])
                elif dependency.startswith("after") and new_job_id is not None:
                    dep_type, dependencies = dependency.split(":")
                    dependencies = ",".join([new_job_id if d==str(job_id) else d for d in dependencies.split(",")])
                    subprocess.call(["scontrol", "update", "JobId={}".format(_job_id), "Dependency={}:{}".format(dep_type, dependencies)])

    def add_individual_parameters(self):
        for param in self.parameters:
            self += "#SBATCH {}\n".format(param)
        if self.modules is not None:
            for module in self.modules.split(","):
                self += "module load {}\n".format(module)

    def write(self):
        with open(self.cmd_file, "w") as f:
            f.write(self.cmd)

    def submit_individual(self, dependency=None, hold_jid=None, update_dependencies=False):
        self.cmd_file.close()

        # while not SwarmJob.can_add_job():
        #     time.sleep(0.3)

        cmd = ["sbatch"]

        if dependency is not None:
            cmd.append("--dependency={}".format(dependency))
        elif hold_jid is not None:
            cmd.append("--dependency=afterany:{}".format(",".join(hold_jid)))

        cmd.append(self.cmd_file_name)

        for _ in xrange(5):
            try:
                self.job_id = subprocess.check_output(cmd).rstrip()
                break
            except subprocess.CalledProcessError:
                time.sleep(0.3)

        if update_dependencies:
            SwarmJob.update_dependencies(self.name, new_job_id=self.job_id)

        time.sleep(0.3)
        return self.job_id

    def run(self, filter_unique=False, split_jobs=False, dependency=None, update_dependencies=False):
        self.cmd_file.close()

        if filter_unique:
            with open("{}.uniq".format(self.cmd_file_name), "w") as uniq:
                subprocess.call(["uniq", self.cmd_file_name], stdout=uniq)
            self.cmd_file_name = "{}.uniq".format(self.cmd_file_name)

        if split_jobs:
            num_subjobs = int(subprocess.check_output(["wc", "-l", self.cmd_file_name]).split()[0])-1
            if num_subjobs/100. > 100.:
                bundle_size = num_subjobs/(2*1000)
                split_by = 2
                while bundle_size > 60:
                    split_by += 1
                    bundle_size = num_subjobs/(split_by*1000)
                print "Splitting Swarm into", split_by, "jobs"
                subprocess.call("tail -n +2 {0} > {0}.nohead".format(self.cmd_file_name), shell=True)
                subprocess.call(["split", "-d", "-l{}".format(int(math.ceil(float(num_subjobs/split_by)))), self.cmd_file_name+".nohead", os.path.splitext(os.path.basename(self.cmd_file_name))[0]+".split"])
                os.remove(self.cmd_file_name+".nohead")

                jobs = []
                for job_num in xrange(split_by):
                    split_file = "{}.split0{}".format(os.path.splitext(os.path.basename(self.cmd_file_name))[0], job_num)
                    subprocess.call("head -n 1 {0} | cat - {0} > {1}.fix".format(self.cmd_file_name, split_file), shell=True)
                    os.remove(split_file)
                    shutil.move(split_file+".fix", split_file)

                    job = SwarmJob(
                        self.name,
                        cpus=self.cpus,
                        mem=self.mem,
                        walltime=self.walltime,
                        gpus=self.gpus,
                        user_parameters=self.user_parameters,
                        modules=self.modules,
                        merge_output=self.merge_output,
                        threads_per_job=self.threads_per_job,
                        cmd_file=split_file)
                    jid = job.run(dependency=dependency)
                    jobs.append(jid)
                return ",".join(jobs)

        if SWARM_AVAIL:
            cmd = ["swarm", "--file", self.cmd_file_name,
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

            if dependency is not None:
                cmd += ["--dependency={}".format(dependency)]
        else:
            num_subjobs = int(subprocess.check_output(["wc", "-l", self.cmd_file_name]).split()[0])-1
            bundle_size = 1
            bundle_size_inc = math.log10(num_subjobs)
            num_jobs = num_subjobs
            while num_jobs < 10000:
                bundle_size += bundle_size_inc
                num_jobs = (num_subjobs/bundle_size)+1

            self.user_parameters["--array"] = "0-{}".format(num_jobs)
            
            job = SwarmJob(
                self.name, 
                cpus=self.cpus,
                mem=self.mem,
                walltime=self.walltime,
                gpus=self.gpus,
                modules=self.modules,
                merge_output=self.merge_output,
                threads_per_job=self.threads_per_job,
                individual = True,
                user_parameters=self.user_parameters)
            
            for i in xrange(bundle_size):
                job += "eval `awk 'NR==${{SLURM_ARRAY_TASK_ID}}*{}+{} {{print;exit}}' {}`".format(bundle_size, i+2, self.cmd_file_name)
            job.run()

        for _ in xrange(5):
            try:
                self.job_id = subprocess.check_output(cmd).rstrip()
                break
            except subprocess.CalledProcessError:
                time.sleep(0.3)

        if update_dependencies:
            SwarmJob.update_dependencies(self.name, new_job_id=self.job_id)

        time.sleep(0.3)
        return self.job_id

def calculate_features(dataset, pdb, chain, sdi=None, domain=None):
    from molmimic.biopdbtools import Structure
    print pdb, chain
    #Set paths correctly so Structure can find PDBs and Features
    initialize_data(dataset)
    Structure.features_from_string(pdb, chain, sdi=sdi, domain=domain, force_feature_calculation=True, grid=False)

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
