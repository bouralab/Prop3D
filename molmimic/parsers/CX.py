import subprocess
from joblib import Memory

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None

memory = Memory(verbose=0)

@memory.cache
def run_cx(struct, modified=None):
    full_pdb_path = struct.path
    pdb_path = os.path.basename(full_pdb_path)

    if apiDockerCall is not None and struct.job is not None:
        work_dir = job.fileStore.getLocalTempDir()
        parameters = [os.path.join("/data", pdb_path), os.path.join("/data", pdb_path+".cx")]
        apiDockerCall(struct.job,
                      image='edraizen/cx:latest',
                      working_dir=work_dir,
                      parameters=parameters)
        cx_f = open(os.path.join(work_dir, pdb_path+".cx"))
        delCx = os.path.join(work_dir, pdb_path+".cx")
    else:
        with open(full_pdb_path) as f:
            cx_f = subprocess.check_output("cx", stdin=f)
        cx_f = iter(cx_f.splitlines())
        delCx = None

    #Read in b-factor from PDB file. CX sometimes introduces invalid characters
    #so the Bio.PDB parser cannot be used
    result = {}
    for l in lines:
        if l[:6].strip() in ["ATOM", "HETATM"]:
            try:
                result[int(l[6:11].strip())] = float(l[60:66].strip())
            except ValueError as e:
                print("    Error, maybe the next line contains it?")

    if delCx is not None:
        cx_f.close()
        os.remove(delCx)

    return result
