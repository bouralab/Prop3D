import os
import re
from molmimic.parsers.superpose import Superpose

from toil.realtimeLogger import RealtimeLogger

def fix_ter(pdb_file):
    if not os.path.isfile(pdb_file+".noter") or os.stat(pdb_file+".noter").st_size == 0:
        print(f"Making file {pdb_file}.noter")
        remove_ter_lines(pdb_file, pdb_file+".noter")
    return pdb_file+".noter"

class mTMAlign(Superpose):
    IMAGE = "docker://edraizen/mtmalign"
    ARG_START="-"
    RULES = {"make_file_list": "make_file_list"}
    DETACH=True
    PARAMETERS = [
        ("input_list", "make_file_list", ["-i", "{}"]),
        (":out_file", "str", "o"),
    ]
    EXTRA_CONTAINER_KWDS = {"shm_size":"64g"}
    LOCAL = ["/home/bournelab/mTM-align/mTM-align"]
    FORCE_LOCAL = True

    def make_file_list(self, key, file_list, format_in_path=True):
        file_list = [fix_ter(pdb) for pdb in file_list]
        return super().make_file_list(key, file_list, format_in_path=format_in_path)

    def all_vs_all(self, pdb_list, out_file="results.pdb"):
        out = self(input_list=pdb_list)
        print(out_file, out)
        return out

    def one_vs_all(self, experiment, pdb_list, table_out_file=None, in_all_vs_all=False, **kwds):
        raise RuntimeError("mTMAlign can only be run with all vs all")

    def one_vs_one(self, moving_pdb_file, fixed_pdb_file, table_out_file=None, **kwds):
        raise RuntimeError("mTMAlign can only be run with all vs all")
