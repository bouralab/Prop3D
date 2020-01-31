import os
import shutil

from molmimic.parsers.container import Container

class CX(Container):
    IMAGE = 'docker://edraizen/cx:latest'
    LOCAL = ["cx"]
    PARAMETERS = [("in_file", "path:in:stdin")]

    def get_concavity(self, pdb_file):
        cx_f = self(in_file=pdb_file)
        return self.parse_cx(cx_f.splitlines())

    @staticmethod
    def parse_cx(cx_file):
        """cx_file: FIle-like object or list of lines
        """
        #Read in b-factor from PDB file. CX sometimes introduces invalid characters
        #so the Bio.PDB parser cannot be used
        result = {}
        for l in cx_file:
            if l[:6].strip() in ["ATOM", "HETATM"]:
                try:
                    result[int(l[6:11].strip())] = float(l[60:66].strip())
                except ValueError as e:
                    pass

        return result
