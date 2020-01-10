import os
from joblib import Memory
from Bio.PDB import DSSP as bioDSSP

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None

from docker.errors import ContainerError

from molmimic.generate_data.util import delocc_pdb

memory = Memory(verbose=0)

class DSSP(Container):
    IMAGE = 'edraizen/dssp:latest'
    LOCAL = ["dssp"]
    PARAMETERS = ["-i", ("in_file", "path:in"), ("out_file", "path:out")]

    def get_dssp(self, bioPDB, pdb_path, clean=True):
        dssp = self(in_file=pdb_path, out_file=pdb_path+".dssp")
        result = bioDSSP(bioPDB, dssp, file_type='DSSP')

        if clean:
            try:
                os.remove(dssp)
            except OSError:
                pass

# @memory.cache
# def run_dssp(biopdb, pdb_path, work_dir=None, job=None):
#     work_dir = work_dir or os.getcwd()
#     remove = False
#     if apiDockerCall is not None and job is not None:
#         if not os.path.abspath(os.path.dirname(pdb_path)) == os.path.abspath(work_dir):
#             shutil.copy(pdb_file, work_dir)
#             remove = True
#         pdb_path = os.path.join(work_dir, os.path.basename(pdb_path))
#         parameters = ['-i', os.path.join("/data", os.path.basename(pdb_path)),
#             '-o', os.path.join("/data", os.path.basename(pdb_path)+".dssp")]
#         try:
#             apiDockerCall(job,
#                           image='edraizen/dssp:latest',
#                           working_dir=work_dir,
#                           parameters=parameters)
#         except ContainerError as e:
#             if "DSSP could not be created due to an error" in str(e):
#                 #Try again without altlocs
#                 delocc_file = delocc_pdb(pdb_path)
#                 parameters = ['-i', os.path.join("/data", os.path.basename(delocc_file)),
#                     '-o', os.path.join("/data", os.path.basename(delocc_file)+".dssp")]
#                 apiDockerCall(job,
#                               image='edraizen/dssp:latest',
#                               working_dir=work_dir,
#                               parameters=parameters)
#         dssp = DSSP(biopdb[0], os.path.join(work_dir,
#             os.path.basename(pdb_path)+".dssp"), file_type='DSSP')
#     else:
#         try:
#             dssp = DSSP(biopdb, pdb_path, dssp='dssp')
#         except KeyboardInterrupt:
#             raise
#         except NameError:
#             dssp = None
#         except Exception as e:
#             raise InvalidPDB("Cannot run dssp for {}".format(pdb_path))
#
#     if remove:
#         try:
#             os.remove(pdb_path)
#         except OSError:
#             pass
#
#     return dssp
