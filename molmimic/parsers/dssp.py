import os
from Bio.PDB import DSSP as bioDSSP

from molmimic.parsers.container import Container
from molmimic.util.pdb import delocc_pdb

class DSSP(Container):
    IMAGE = 'docker://edraizen/dssp:latest'
    LOCAL = ["dssp"]
    PARAMETERS = [
        ("in_file", "path:in", "i"),
        ("out_file", "path:out", "o")]
    RETURN_FILES=True
    ARG_START="-"


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
