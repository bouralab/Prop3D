import os
import shutil

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

def run_cx(pdb_path, work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    if apiDockerCall is not None and job is not None:
        if not os.path.abspath(os.path.dirname(pdb_path)) == os.path.abspath(work_dir):
            shutil.copy(pdb_path, work_dir)
        parameters = [os.path.basename(pdb_path)]
        print "BEFORE", os.listdir(work_dir)
        with open(pdb_path) as f:
            print "IT EXISTS", next(f)
        cx_f = apiDockerCall(job,
                             image='edraizen/cx:latest',
                             working_dir="/data",
                             volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                             parameters=parameters)
        print "DONE RUNINGF CX", cx_f
        # cx_f = open(os.path.join(work_dir, os.path.basename(pdb_path)+".cx"))
        # delCx = os.path.join(work_dir, os.path.basename(pdb_path)+".cx")
    else:
        with open(pdb_path) as f:
            cx_f = subprocess.check_output("cx", stdin=f)
    cx_f = iter(cx_f.splitlines())

    #Read in b-factor from PDB file. CX sometimes introduces invalid characters
    #so the Bio.PDB parser cannot be used
    result = {}
    for l in cx_f:
        if l[:6].strip() in ["ATOM", "HETATM"]:
            try:
                result[int(l[6:11].strip())] = float(l[60:66].strip())
            except ValueError as e:
                pass

    # if delCx is not None:
    #     cx_f.close()
    #     os.remove(delCx)

    return result
