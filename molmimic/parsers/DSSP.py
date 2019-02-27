import os
from joblib import Memory
from Bio.PDB import DSSP

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None

memory = Memory(verbose=0)

@memory.cache
def run_dssp(biopdb, pdb_path, work_dir=None, job=None):
    work_dir = work_dir or os.getcwd()

    if apiDockerCall is not None and job is not None:
        if not os.path.abspath(os.path.dirname(pdb_path)) == os.path.abspath(work_dir):
            shutil.copy(pdb_file, work_dir)
        parameters = ['-i', os.path.join("/data", os.path.basename(pdb_path)),
            '-o', os.path.join("/data", os.path.basename(pdb_path)+".dssp")]
        apiDockerCall(job,
                      image='edraizen/dssp:latest',
                      working_dir=work_dir,
                      parameters=parameters)
        dssp = DSSP(biopdb[0], os.path.join(work_dir,
            os.path.basename(pdb_path)+".dssp"), file_type='DSSP')
    else:
        try:
            dssp = DSSP(biopdb, pdb_path, dssp='dssp')
        except KeyboardInterrupt:
            raise
        except NameError:
            dssp = None
        except Exception as e:
            raise InvalidPDB("Cannot run dssp for {}".format(pdb_path))

    return dssp
