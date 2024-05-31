import os, sys
import subprocess
import shutil

from joblib import Memory

from Prop3D.parsers.container import Container

class SCWRL(Container):
    IMAGE = 'docker://edraizen/scwrl4:latest'
    ENTRYPOINT = "/opt/scwrl4/Scwrl4"
    LOCAL = ["scwrl4"]
    PARAMETERS = [
        ("in_file", "path:in", "i"),
        ("out_file", "path:out", "o"),
        (":frame_file", "path:in", "f"),
        (":seq_file", "path:in", "s"),
        (":param_file", "path:in", "p"),
        (":in_crystal", "store_true", "#"),
        (":remove_hydrogens", "store_true", "h"),
        (":remove_h_n_term", "store_true", "t"),
        ]
    RETURN_FILES = True
    ARG_START="-"

    def fix_rotamers(self, in_file, out_file=None, **kwds):
        if out_file is None:
            out_file = str(Path(self.work_dir) / f"{Path(pdb_file).stem}.scwrl.pdb")
        return self(in_file, out_file, **kwds)
