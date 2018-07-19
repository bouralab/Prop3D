import subprocess

from joblib import Memory
import pybel
from Bio.PDB import DSSP

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

@memory.cache
def run_dssp(struct, modified=False):
    if modified:
        pdbfd, tmp_pdb_path = tempfile.mkstemp()
        with os.fdopen(pdbfd, 'w') as tmp:
            struct.save_pdb(tmp)
    else:
        tmp_pdb_path = struct.path

    try:
        dssp = DSSP(self.structure[0], tmp_pdb_path, dssp='dssp')
    except KeyboardInterrupt:
        raise
    except NameError:
        dssp = None
    except Exception as e:
        raise InvalidPDB("Cannot run dssp for {}".format(self.id))

    if modified:
        os.remove(tmp_pdb_path)

    return dssp
