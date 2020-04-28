import os
from joblib import Memory
from Bio.PDB import DSSP as bioDSSP

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None

from molmimic.parsers.container import Container
from molmimic.util.pdb import delocc_pdb

memory = Memory(verbose=0)

class DSSP(Container):
    IMAGE = 'docker://edraizen/dssp:latest'
    LOCAL = ["dssp"]
    PARAMETERS = ["-i", ("in_file", "path:in"), "-o", ("out_file", "path:out")]
    RETURN_FILES=True


    def get_dssp(self, bioPDB, pdb_path, out_file=None, clean=True):
        if out_file is None:
            out_file = os.path.join(self.work_dir, os.path.basename(pdb_path)+".dssp")

        try:
            dssp = self.__call__(in_file=pdb_path, out_file=out_file)
        except RuntimeError as e:
            if "DSSP could not be created due to an error" in str(e):
                delocc_file = delocc_pdb(pdb_path)
                self.files_to_remove.append(delocc_file)
                out_file = delocc_file+".dssp"
                try:
                    dssp = self.__call__(in_file=delocc_file, out_file=out_file)
                except RuntimeError as e:
                    raise
            raise

        result = bioDSSP(bioPDB[0], dssp, file_type='DSSP')

        if clean:
            try:
                os.remove(out_file)
            except OSError:
                pass

        return result

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
