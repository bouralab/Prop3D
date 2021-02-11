import os

import numpy as np
from Bio.PDB import PDBParser
from molmimic.util import safe_remove
from molmimic.parsers.container import Container

class Multivalue(Container):
    IMAGE = 'docker://edraizen/multivalue:latest'
    LOCAL = ["multivalue"]
    PARAMETERS = [
        ("coordinates_file", "path:in", ""),
        ("dx_file", "path:in", ""),
        ("out_file", "path:out", "")
        ]
    RETURN_FILES = True

    def compute_electrostatics_at_coordinates(self, coordinates, dx_file, out_file=None):
        """Modified from MaSIFs computeAPBS"""

        file_base = os.path.splitext(dx_file)[0]

        if out_file is None:
            out_file = f"{file_base}_out.csv"

        if isinstance(coordinates, (list, tuple, np.ndarray)):
            coordinates_file = os.path.join(self.work_dir, f"{file_base}.csv")
            with open(coordinates_file, "w") as f:
                for coord in coordinates:
                    print(f"{coord[0]},{coord[1]},{coord[2]}", file=f)
        elif isinstance(coordinates, str) and os.path.isfile(coordinates):
            coordinates_file = coordinates
        else:
            raise RuntimeError(f"Invlaid vertices -- must be list, array, or file, not {type(coordinates)}")

        print(dx_file)

        out_file = self(coordinates_file=coordinates_file, dx_file=dx_file,
            out_file=out_file)

        if isinstance(coordinates, (list,tuple)):
            safe_remove(coordinates_file)

        return out_file

    def get_electrostatics_at_coordinates(self, coordinates, dx_file, out_file=None):
        charge_file = self.compute_electrostatics_at_coordinates(coordinates,
            dx_file, out_file=out_file)

        return self.read_charges(charge_file)

    @staticmethod
    def read_charges(charge_file):
        """Modified from MaSIFs computeAPBS"""
        # charges = np.array([0.0] * len(vertices))
        # with open(charge_file) as f:
        #     for ix, line in enumerate(f):
        #         charges[ix] = float(line.split(",")[3])
        with open(charge_file) as f:
            charges = np.array([float(line.split(",")[3]) for line in f])

        return charges
