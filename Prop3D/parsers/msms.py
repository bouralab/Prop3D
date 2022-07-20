import os
import random

import numpy as np
from Bio.PDB import PDBParser
from Prop3D.util import safe_remove
from Prop3D.parsers.container import Container

#Values from MaSIF
# radii for atoms in explicit case.
_radii = {
    "N": "1.540000", "N": "1.540000", "O": "1.400000", "C": "1.740000",
    "H": "1.200000", "S": "1.800000", "P": "1.800000", "Z": "1.39",
    "X": "0.770000"}  ## Radii of CB or CA in disembodied case.

# This  polar hydrogen's names correspond to that of the program Reduce.
_polarHydrogens = {
    "ALA": ["H"], "GLY": ["H"], "SER": ["H", "HG"], "THR": ["H", "HG1"],
    "LEU": ["H"], "ILE": ["H"], "VAL": ["H"], "ASN": ["H", "HD21", "HD22"],
    "GLN": ["H", "HE21", "HE22"], "ARG": ["H", "HH11", "HH12", "HH21", "HH22", "HE"],
    "HIS": ["H", "HD1", "HE2"], "TRP": ["H", "HE1"], "PHE": ["H"],
    "TYR": ["H", "HH"], "GLU": ["H"], "ASP": ["H"],
    "LYS": ["H", "HZ1", "HZ2", "HZ3"], "PRO": [], "CYS": ["H"], "MET": ["H"]}

class MSMS(Container):
    IMAGE = 'docker://edraizen/msms:latest'
    ENTRYPOINT = "/usr/local/bin/msms"
    LOCAL = ["msms"]
    PARAMETERS = [
        ("in_file", "path:in", ["-if", "{}"]),
        ("out_file", "path:in", ["-of", "{}"]),
        (":asa_file", "path:in", ["-af", "{}"]),
        (":probe:1.5", "str"),
        (":no_header", "store_true"),
        (":density:3.0", "str"),
        (":hdensity:3.0", "str"),
        (":no_area", "store_true"),
        (":surface", "str"),
        (":socket", "str"),
        (":sinetd", "store_true"),
        (":noh", "store_true"),
        (":no_rest_on_pbr", "store_true"),
        (":no_rest", "store_true"),
        (":free_vertices", "store_true"),
        (":all_components", "store_true"),
        (":one_cavity", "str"),
        ]
    ARG_START = "-"

    def compute_surface_from_pdb(self, pdb_file, **kwds):
        xyzrnfilename = self.output_pdb_as_xyzrn(pdb_file,
            radii=kwds.pop("radii", None),
            polarHydrogens=kwds.pop("polarHydrogens", None))
        file_base = os.path.splitext(xyzrnfilename)[0]
        self(in_file=xyzrnfilename, out_file=file_base, asa_file=file_base, **kwds)
        return file_base

    def get_surface_and_area_from_pdb(self, pdb_file, return_msms_file=False, **kwds):
        file_base = self.compute_surface_from_pdb(pdb_file, **kwds)

        vert_file = f"{file_base}.vert"
        assert os.path.isfile(vert_file), "Vert file does not exist"
        area_file = f"{file_base}.area"
        assert os.path.isfile(area_file), "Area file does not exist"

        return self.read_msms(file_base, return_msms_file=return_msms_file)

    def output_pdb_as_xyzrn(self, pdbfilename, xyzrnfilename=None, radii=None, polarHydrogens=None):
        """Read a pdb file and output it is in xyzrn for use in MSMS

        Pablo Gainza - LPDI STI EPFL 2019
        This file is part of MaSIF.
        Released under an Apache License 2.0

        Parameters
        ----------
            pdbfilename: input pdb filename
            xyzrnfilename: output in xyzrn format.
        """
        radii = radii if isinstance(radii, dict) else _radii
        polarHydrogens = polarHydrogens if isinstance(polarHydrogens, dict) else _polarHydrogens

        parser = PDBParser()
        struct = parser.get_structure(pdbfilename, pdbfilename)

        if xyzrnfilename is None:
            randnum = random.randint(1,10000000)
            xyzrnfilename = os.path.join(self.work_dir, f"msms_{str(randnum)}.xyzrn")

        with open(xyzrnfilename, "w") as outfile:
            for atom in struct.get_atoms():
                name = atom.get_name()
                residue = atom.get_parent()
                # Ignore hetatms.
                if residue.get_id()[0] != " ":
                    continue
                resname = residue.get_resname()
                reskey = residue.get_id()[1]
                chain = residue.get_parent().get_id()
                atomtype = name[0]

                color = "Green"
                coords = None
                if atomtype in radii and resname in polarHydrogens:
                    if atomtype == "O":
                        color = "Red"
                    if atomtype == "N":
                        color = "Blue"
                    if atomtype == "H":
                        if name in polarHydrogens[resname]:
                            color = "Blue"  # Polar hydrogens
                    coords = "{:.06f} {:.06f} {:.06f}".format(
                        atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
                    )
                    insertion = "x"
                    if residue.get_id()[2] != " ":
                        insertion = residue.get_id()[2]
                    full_id = "{}_{:d}_{}_{}_{}_{}".format(
                        chain, residue.get_id()[1], insertion, resname, name, color
                    )
                if coords is not None:
                    outfile.write(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")

        return xyzrnfilename

    @staticmethod
    def read_msms(file_root, return_msms_file=False):
        vertices, faces, normals, names = MSMS.read_surface(file_root)
        areas = MSMS.read_areas(file_root)

        if return_msms_file:
            return vertices, faces, normals, names, areas, file_root

        return vertices, faces, normals, names, areas


    @staticmethod
    def read_surface(file_root):
        """
        Read an msms output file containing the surface that was output by MSMS.
        MSMS outputs two files: {file_root}.vert and {file_root}.face

        Pablo Gainza - LPDI STI EPFL 2019
        Released under an Apache License 2.0
        """
        vertfile = open(file_root + ".vert")
        meshdata = (vertfile.read().rstrip()).split("\n")
        vertfile.close()

        # Read number of vertices.
        count = {}
        header = meshdata[2].split()
        count["vertices"] = int(header[0])
        ## Data Structures
        vertices = np.zeros((count["vertices"], 3))
        normalv = np.zeros((count["vertices"], 3))
        atom_id = [""] * count["vertices"]
        res_id = [""] * count["vertices"]
        for i in range(3, len(meshdata)):
            fields = meshdata[i].split()
            vi = i - 3
            vertices[vi][0] = float(fields[0])
            vertices[vi][1] = float(fields[1])
            vertices[vi][2] = float(fields[2])
            normalv[vi][0] = float(fields[3])
            normalv[vi][1] = float(fields[4])
            normalv[vi][2] = float(fields[5])
            atom_id[vi] = fields[7]
            res_id[vi] = fields[9]
            count["vertices"] -= 1

        # Read faces.
        facefile = open(file_root + ".face")
        meshdata = (facefile.read().rstrip()).split("\n")
        facefile.close()

        # Read number of vertices.
        header = meshdata[2].split()
        count["faces"] = int(header[0])
        faces = np.zeros((count["faces"], 3), dtype=int)
        normalf = np.zeros((count["faces"], 3))

        for i in range(3, len(meshdata)):
            fi = i - 3
            fields = meshdata[i].split()
            faces[fi][0] = int(fields[0]) - 1
            faces[fi][1] = int(fields[1]) - 1
            faces[fi][2] = int(fields[2]) - 1
            count["faces"] -= 1

        assert count["vertices"] == 0
        assert count["faces"] == 0

        return vertices, faces, normalv, res_id

    @staticmethod
    def read_areas(file_base):
        areas = {}
        with open(file_base+".area") as ses_file:
            next(ses_file) # ignore header line
            for line in ses_file:
                fields = line.split()
                areas[fields[3]] = fields[1]
        return areas
