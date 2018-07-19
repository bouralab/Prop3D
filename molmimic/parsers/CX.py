import subprocess
from tempfile import mkdtemp

from joblib import Memory

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

@memory.cache
def run_cx(struct, modified=False):
    if modified:
        pdbfd, tmp_pdb_path = tempfile.mkstemp()
        with os.fdopen(pdbfd, 'w') as tmp:
            struct.save_pdb(tmp, True)
    else:
        tmp_pdb_path = struct.path

    with open(tmp_pdb_path) as f:
        cx_f = subprocess.check_output("cx", stdin=f)

    if modified:
        os.remove(tmp_pdb_path)

    #Read in b-factor from PDB file. CX sometimes introduces invalid characters
    #so the Bio.PDB parser cannot be used
    result = {}
    lines = iter(cx_f.splitlines())
    for l in lines:
        if l[:6].strip() in ["ATOM", "HETATM"]:
            try:
                result[int(l[6:11].strip())] = float(l[60:66].strip())
            except ValueError as e:
                print("    Error, maybe the next line contains it?")

    return result
