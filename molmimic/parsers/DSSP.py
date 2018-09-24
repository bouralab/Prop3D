import os
from joblib import Memory
from Bio.PDB import DSSP

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None

memory = Memory(verbose=0)

@memory.cache
def run_dssp(struct, modified=None):
    full_pdb_path = struct.path
    pdb_path = os.path.basename(full_pdb_path)
    work_dir = os.path.dirname(full_pdb_path)

    if apiDockerCall is not None and struct.job is not None:
        work_dir = job.fileStore.getLocalTempDir()
        parameters = ['-i', os.path.join("/data", pdb_path), '-o', os.path.join("/data", pdb_path+".dssp")]
        apiDockerCall(struct.job,
                      image='edraizen/dssp:latest',
                      working_dir=work_dir,
                      parameters=parameters)
        dssp = DSSP(struct.structure[0], os.path.join(work_dir, pdb_path+".dssp"), dssp='dssp')
    else:
        try:
            dssp = DSSP(struct.structure[0], full_pdb_path, dssp='dssp')
        except KeyboardInterrupt:
            raise
        except NameError:
            dssp = None
        except Exception as e:
            raise InvalidPDB("Cannot run dssp for {}".format(self.id))

    return dssp
