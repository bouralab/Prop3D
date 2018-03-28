import os
import gzip
import itertools as it
import re
import subprocess
from cStringIO import StringIO
from collections import Iterable, Counter, defaultdict
import tempfile

import numpy as np
import pandas as pd
from scipy import spatial
from sklearn.decomposition import PCA
from Bio import PDB
from Bio import SeqIO
from Bio.PDB.NeighborSearch import NeighborSearch

try:
    import freesasa
except ImportError:
    freesasa = None

try:
    import pybel
except ImportError:
    pybabel = None

from get_pdb import get_pdb
from molmimic.util import silence_stdout, silence_stderr
from molmimic.map_residues import map_residues

try:
    from numba import jit
except ImportError:
    jit = lambda x: x

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

hydrophobicity_scales = {
    "kd": {'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
           'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
           'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
           'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 },
    "biological": {
           "A": -0.11,"C":  0.13,"D": -3.49,"E": -2.68,"F": 0.32,
           "G": -0.74,"H": -2.06,"I":  0.60,"K": -2.71,"L": 0.55,
           "M":  0.10,"N": -2.05,"P": -2.23,"Q": -2.36,"R": -2.58,
           "S": -0.84,"T": -0.52,"V":  0.31,"W": -0.30,"Y": -0.68},
    "octanal":{
           "A": -0.50, "C":  0.02, "D": -3.64, "E": -3.63,
           "F":  1.71, "G": -1.15, "H": -0.11, "I":  1.12,
           "K": -2.80, "L":  1.25, "M":  0.67, "N": -0.85,
           "P": -0.14, "Q": -0.77, "R": -1.81, "S": -0.46,
           "T": -0.25, "V":  0.46, "W":  2.09, "Y":  0.71,}
    }

#vdw_radii = {'Ru': 1.2, 'Re': 4.3, 'Ra': 2.57, 'Rb': 2.65, 'Rn': 2.3, 'Rh': 1.22, 'Be': 0.63, 'Ba': 2.41, 'Bi': 1.73, 'Bk': 1.64, 'Br': 1.85, 'D': 1.2, 'H': 1.2, 'P': 1.9, 'Os': 1.58, 'Es': 1.62, 'Hg': 1.55, 'Ge': 1.48, 'Gd': 1.69, 'Ga': 1.87, 'Pr': 1.62, 'Pt': 1.72, 'Pu': 1.67, 'C': 1.775, 'Pb': 2.02, 'Pa': 1.6, 'Pd': 1.63, 'Cd': 1.58, 'Po': 1.21, 'Pm': 1.76, 'Ho': 1.61, 'Hf': 1.4, 'K': 2.75, 'He': 1.4, 'Md': 1.6, 'Mg': 1.73, 'Mo': 1.75, 'Mn': 1.19, 'O': 1.45, 'S': 1.8, 'W': 1.26, 'Zn': 1.39, 'Eu': 1.96, 'Zr': 1.42, 'Er': 1.59, 'Ni': 1.63, 'No': 1.59, 'Na': 2.27, 'Nb': 1.33, 'Nd': 1.79, 'Ne': 1.54, 'Np': 1.71, 'Fr': 3.24, 'Fe': 1.26, 'Fm': 1.61, 'B': 1.75, 'F': 1.47, 'Sr': 2.02, 'N': 1.5, 'Kr': 2.02, 'Si': 2.1, 'Sn': 2.17, 'Sm': 1.74, 'V': 1.06, 'Sc': 1.32, 'Sb': 1.12, 'Se': 1.9, 'Co': 1.13, 'Cm': 1.65, 'Cl': 1.75, 'Ca': 1.95, 'Cf': 1.63, 'Ce': 1.86, 'Xe': 2.16, 'Lu': 1.53, 'Cs': 3.01, 'Cr': 1.13, 'Cu': 1.4, 'La': 1.83, 'Li': 1.82, 'Tl': 1.96, 'Tm': 1.57, 'Lr': 1.58, 'Th': 1.84, 'Ti': 1.95, 'Te': 1.26, 'Tb': 1.66, 'Tc': 2.0, 'Ta': 1.22, 'Yb': 1.54, 'Dy': 1.63, 'I': 1.98, 'U': 1.75, 'Y': 1.61, 'Ac': 2.12, 'Ag': 1.72, 'Ir': 1.22, 'Am': 1.66, 'Al': 1.5, 'As': 0.83, 'Ar': 1.88, 'Au': 1.66, 'At': 1.12, 'In': 1.93}
vdw_radii = {
    "H" : 1.2,
    "Li" : 1.82,
    "Na" : 2.27,
    "K" : 2.75,
    "C" : 1.7,
    "N" : 1.55,
    "O" : 1.52,
    "F" : 1.47,
    "P" : 1.80,
    "S" : 1.80,
    "Cl" : 1.75,
    "Br" : 1.85,
    "Se" : 1.90,
    "Zn" : 1.39,
    "Cu" : 1.4,
    "Ni" : 1.63,
}

vdw_aa_radii = {
    'ALA': 2.52,
    'ARG': 3.60,
    'ASN': 3.0,
    'ASP': 2.94,
    'CYS': 2.92,
    'GLN': 3.25,
    'GLU': 3.21,
    'GLY': 2.25,
    'HIS': 3.308,
    'ILE': 3.22,
    'LEU': 3.20,
    'LYS': 3.42,
    'MET': 3.37,
    'PHE': 3.47,
    'PRO': 2.93,
    'SER': 2.67,
    'THR': 2.90,
    'TRP': 3.73,
    'TYR': 3.55,
    'VAL': 3.01
 }

surface_areas = {atom:4.*np.pi*(radius**2) for atom, radius in vdw_radii.items()}

maxASA = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0, "K": 201.0, "L": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}

class InvalidPDB(RuntimeError):
    pass

class Structure(object):
    def __init__(self, path, pdb, chain, id=None, course_grained=False, volume=264, voxel_size=1.0, rotate=True, input_format="pdb", force_feature_calculation=False, unclustered=False):
        self.path = path

        if not os.path.isfile(self.path):
            raise InvalidPDB("Cannot find file {}".format(self.path))

        if self.path.endswith(".gz"):
            raise InvalidPDB("Gzipped archives not allowed. Please use constructor or util.get_pdb")

        self.input_format = input_format
        if self.input_format in ["pdb", "pqr"]:
            parser = PDB.PDBParser()
        elif self.input_format == "mmcif":
            parser = PDB.FastMMCIFParser()
        elif self.input_format == "mmtf":
            parser = PDB.MMTFParser()
        else:
            raise RuntimeError("Invalid PDB parser (pdb, mmcif, mmtf)")

        try:
            self.structure = parser.get_structure(pdb, self.path)    
        except KeyError:
            #Invalid mmcif file
            raise InvalidPDB("Invalid PDB file: {} (path={})".format(pdb, self.path))

        self.chain = chain.split("_", 1)[0]

        try:
            all_chains = list(self.structure[0].get_chains())
        except (KeyError, StopIteration):
            raise InvalidPDB("Error get chains for {} {}".format(pdb, self.path))

        if len(all_chains) > 1:
            raise InvalidPDB("Only accepts PDBs with 1 chain")

        only_chain = all_chains[0]
        
        self.modified_pdb_file = False

        if self.chain != only_chain.get_id():
            #Reset chain if Bio.PDB changes name after extraction
            self.chain = only_chain.get_id()

        self.altloc = getattr(self.structure, "altloc", " ")

        self.pdb = pdb
        self.dssp = None
        self.pqr = None
        self.cx = None
        self.sasa = None
        self.sasa_struct = None
        self.mean_coord = np.zeros(3)
        self.mean_coord_updated = False
        self.starting_residue = None
        self.starting_index_seq = None
        self.starting_index_struc = None
        self.volume = volume
        self.id = "{}_{}{}".format(self.pdb, self.chain, "_{}".format(id) if id else "")
        self.course_grained = course_grained
        self.nFeatures = Structure.number_of_features(course_grained=course_grained)

        self.voxel_size = voxel_size

        if not rotate:
            self.shift_coords_to_volume_center()
            self.set_voxel_size(voxel_size)
        else:
            next(self.rotate())

        clustered_path = "features2" if unclustered else "features3"
        if force_feature_calculation:
            precalc_features_path = os.environ.get("MOLMIMIC_FEATURES", os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", clustered_path, "residue" if self.course_grained else "atom"))
            precalc_features_path = os.path.join(precalc_features_path, "{}.npy".format(self.id))
        else:
            precalc_features_path = os.environ.get("MOLMIMIC_FEATURES", os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", clustered_path, "residue" if self.course_grained else "atom"))
            precalc_features_path = os.path.join(precalc_features_path, "{}.npy".format(self.id))

        if course_grained:
            shape = max([r.get_id()[1] for r in self.structure.get_residues()])
        else:
            shape = max([a.serial_number for a in self.structure.get_atoms()])

        self.precalc_features = None
        if force_feature_calculation:
            print "force features"
            self.precalc_features = np.memmap(precalc_features_path, dtype=np.float, mode="w+", shape=(shape, self.nFeatures))
            self.force_feature_calculation = True
        else:
            try:
                self.precalc_features = np.memmap(precalc_features_path, dtype=np.float, mode="r", shape=(shape, self.nFeatures))
                self.force_feature_calculation = False
            except IOError as e:
                print e
                self.precalc_features = None
                self.force_feature_calculation = True

    @staticmethod
    def number_of_features(only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, course_grained=False):
        if course_grained:
            if non_geom_features:
                return 29
            elif only_aa:
                return 21
            else:
                return 39
        else:
            if non_geom_features:
                return 13
            elif only_atom:
                return 5
            elif only_aa:
                return 21
            elif use_deepsite_features:
                return 8
            else:
                return 67

    @classmethod
    def from_pdb(cls, pdb, chain, id=None, input_shape=(96,96,96), batch_index=None, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, course_grained=False, grid=True, return_data=True, return_truth=True, rotate=True, force_feature_calculation=False, expand_atom=False, include_full_protein=False, undersample=False, unclustered=False):
        path, pdb, chain, input_format = get_pdb(pdb, chain)
        return cls(
            path,
            pdb,
            chain,
            id=id,
            volume=np.max(input_shape),
            course_grained=course_grained,
            force_feature_calculation=force_feature_calculation,
            input_format=input_format,
            unclustered=unclustered,
            rotate=rotate)

    @staticmethod
    def features_from_string(pdb, chain, resi=None, id=None, input_shape=(96,96,96), batch_index=None, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, course_grained=False, grid=True, return_data=True, return_truth=True, rotate=True, force_feature_calculation=False, expand_atom=False, include_full_protein=False, undersample=False, unclustered=False):
        """Get features

        Paramters
        ---------
        pdb : str
            PDB ID
        chain : str
            Chain to extract
        resi : str
            Residue numbers separated by whitespace corresponding to positionts in the full protein sequence
        input_shape : 3-tuple
        rotations : int
        """
        path, pdb, chain, input_format = get_pdb(pdb, chain)
        s = Structure(
            path,
            pdb,
            chain,
            id=id,
            volume=np.max(input_shape),
            course_grained=course_grained,
            force_feature_calculation=force_feature_calculation,
            input_format=input_format,
            unclustered=unclustered,
            rotate=rotate)

        #s.orient_to_pai()
        if resi:
            binding_site = s.align_seq_to_struc(resi, return_residue=True)
        else:
            binding_site = None

        if grid:
            try:
                features = s.get_features(
                    input_shape=input_shape,
                    residue_list=binding_site,
                    batch_index=batch_index,
                    only_aa=only_aa,
                    only_atom=only_atom,
                    non_geom_features=non_geom_features,
                    use_deepsite_features=use_deepsite_features,
                    expand_atom=expand_atom,
                    include_full_protein=include_full_protein,
                    undersample=undersample)
            except Exception as e:
                print(e)
                raise
        else:
            features = [s.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features) \
                for atom in s.structure.get_atoms()]

        if s.precalc_features is not None and force_feature_calculation:
            s.precalc_features.flush()

        return features

    def get_atoms(self, include_hetatms=False, exclude_atoms=None):
        for a in self.filter_atoms(self.structure.get_atoms(), include_hetatms=include_hetatms, exclude_atoms=exclude_atoms):
            yield a

    def filter_atoms(self, atoms, include_hetatms=False, exclude_atoms=None):
        for a in atoms:
            hetflag, resseq, icode = a.get_parent().get_id()
            if not include_hetatms and hetflag is not ' ':
                continue
            if exclude_atoms is not None and a.get_name().strip() in exclude_atoms:
                continue
            yield a

    def save_pdb(self, path=None, file_like=False):
        lines = not path and not file_like
        if path is None:
            path = StringIO()

        writer = PDB.PDBIO()
        writer.set_structure(self.structure)
        writer.save(path)

        if file_like:
            path.seek(0)

        if lines:
            path = path.read()

        return path

    def _get_dssp(self):
        if self.dssp is None:
            if self.modified_pdb_file:
                pdbfd, tmp_pdb_path = tempfile.mkstemp()
                with os.fdopen(pdbfd, 'w') as tmp:
                    self.save_pdb(tmp)
            else:
                tmp_pdb_path = self.path

            try:
                self.dssp = PDB.DSSP(self.structure[0], tmp_pdb_path, dssp='mkdssp')
            except NameError:
                self.dssp = None

            if self.modified_pdb_file:
                os.remove(tmp_pdb_path)

        return self.dssp

    def _get_sasa(self):
        if freesasa is None:
            print("SASA not installed! SASA will be 0")
            return None, None
        if self.sasa is None:
            with silence_stdout(), silence_stderr():
                if self.modified_pdb_file:
                    #Note: need to remove hydrogens
                    #self.sasa_struct = freesasa.structureFromBioPDB(self.structure)
                    structure = Structure()
                    classifier = Classifier()

                    atoms = self.structure.get_atoms()

                    for a in atoms:
                        if a.element == "H":
                            #Ignore Hydrogens
                            continue

                        r = a.get_parent()
                        hetflag, resseq, icode = r.get_id()

                        c = r.get_parent()
                        v = a.get_vector()

                        structure.addAtom(a.get_fullname(), r.get_resname(), resseq, c.get_id(),
                                          v[0], v[1], v[2])

                    structure.setRadiiWithClassifier(classifier)
                else:
                    #Automatically removes hydrogens
                    self.sasa_struct = freesasa.Structure(self.path)
                self.sasa = freesasa.calc(self.sasa_struct)

        return self.sasa, self.sasa_struct

    def _get_pqr(self):
        """Run PDB2PQR to get charge for each atom

        TODO: figure out why unknown residues fro hetatms are absent
        """
        if self.pqr is None:
            if not self.modified_pdb_file and self.input_format == "pqr":
                #Don't recalculate, just use occupancy field!
                self.pqr = {a.get_id():a.occupancy for a in self.get_atoms()}
                return self.pqr

            if not self.modified_pdb_file and self.path.endswith(".min.pdb") and os.path.isfile("{}.pqr.pdb".format(self.path[:-8])):
                remove_pqr = False
                tmp_pqr_path = "{}.pqr.pdb".format(self.path[:-8])
            else:
                remove_pqr = True

                if self.modified_pdb_file:
                    pdbfd, tmp_pdb_path = tempfile.mkstemp()
                    with os.fdopen(pdbfd, 'w') as tmp:
                        self.save_pdb(tmp, True)
                else:
                    tmp_pdb_path = self.path

                _, tmp_pqr_path = tempfile.mkstemp()

                with silence_stdout(), silence_stderr():
                    subprocess.call(["/usr/share/pdb2pqr/pdb2pqr.py", "--ff=amber", "--whitespace", tmp_pdb_path, tmp_pqr_path])

                if self.modified_pdb_file:
                    os.remove(tmp_pdb_path)

            self.pqr = {}
            with open(tmp_pqr_path) as pqr:
                for line in pqr:
                    for line in pqr:
                        if not line.startswith("ATOM  ") or line.startswith("HETATM"): continue

                        fields = line.rstrip().split()
                        if len(fields) == 11:
                            recordName, serial, atomName, residueName, chainID, residueNumber, X, Y, Z, charge, radius = fields
                        elif len(fields) == 10:
                            recordName, serial, atomName, residueName, residueNumber, X, Y, Z, charge, radius = fields
                        else:
                            print fields
                            raise RuntimeError("Invalid PQR File")

                        resseq = int("".join([i for i in residueNumber if i.isdigit()]))

                        icode = "".join([i for i in residueNumber if not i.isdigit()])
                        if icode == "":
                            icode = " "

                        if recordName == "HETATM":  # hetero atom flag
                            if residueName in ["HOH", "WAT"]:
                                hetero_flag = "W"
                            else:
                                hetero_flag = "H_{}".format(residueName)
                        else:
                            hetero_flag = " "
                        residue_id = (hetero_flag, resseq, icode)

                        key = (residue_id, (atomName.strip(), ' '))
                        self.pqr[key] = float(charge)

            if remove_pqr:
                os.remove(tmp_pqr_path)

        return self.pqr

    def _get_cx(self):
        if self.cx is None:
            if self.modified_pdb_file:
                pdbfd, tmp_pdb_path = tempfile.mkstemp()
                with os.fdopen(pdbfd, 'w') as tmp:
                    self.save_pdb(tmp, True)
            else:
                tmp_pdb_path = self.path

            with open(tmp_pdb_path) as f:
                cx_f = subprocess.check_output("cx", stdin=f)

            if self.modified_pdb_file:
                os.remove(tmp_pdb_path)

            #Read in b-factor from PDB file. CX sometimes introduces invalid characters
            #so the Bio.PDB parser cannot be used
            self.cx = {}
            lines = iter(cx_f.splitlines())
            for l in lines:
                if l[:6].strip() in ["ATOM", "HETATM"]:
                    try:
                        self.cx[int(l[6:11].strip())] = float(l[60:66].strip())
                    except ValueError as e:
                        print("    Error, maybe the next line contains it?")
            # self.cx = {int(l[6:11].strip()):float(l[60:66].strip()) \
            #     for l in cx_f.splitlines() if l[:6].strip() in ["ATOM", "HETATM"]}

        return self.cx

    def get_autodock_features(self, atom):
        """Modified from pybel write PDBQT c++"""
        if not hasattr(self, "autodock_features"):
            if self.modified_pdb_file:
                mol = pybel.readstring("pdb", self.save_pdb(file_like=True).read())
            else:
                mol = next(pybel.readfile("pdb", self.path))

            mol.addh()

            pdbqt = mol.write("pdbqt")
            self.autodock_features = {} #defaultdict(lambda: ("  ", False)
            for atom_index in range(mol.OBMol.NumAtoms()):
                a = mol.OBMol.GetAtom(atom_index + 1)

                if a.IsCarbon() and a.IsAromatic():
                    element = 'A'
                elif a.IsOxygen():
                    element = 'OA'
                elif a.IsNitrogen() and a.IsHbondAcceptor():
                    element = 'NA'
                elif a.IsSulfur() and a.IsHbondAcceptor():
                    element ='SA'
                else:
                    element = "".join([c for c in a.GetType() if c.isalnum()])

                self.autodock_features[a.GetIdx()] = (element, a.IsHbondDonor())

        try:
            return self.autodock_features[atom.serial_number]
        except KeyError:
            return "  ", False

    def get_mean_coord(self):
        if not self.mean_coord_updated:
            self.mean_coord = np.mean(self.get_coords(), axis=0)
            self.mean_coord_updated = True
        return self.mean_coord

    def shift_coords_to_origin(self):
        mean_coord = self.get_mean_coord()
        coords = self.get_coords()-mean_coord
        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.mean(coords, axis=0)

    def shift_coords_to_volume_center(self):
        mean_coord = self.get_mean_coord()
        coords = self.get_coords()-mean_coord
        coords += np.array([self.volume/2]*3)
        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.mean(coords, axis=0)

    def get_coords(self, include_hetatms=False, exclude_atoms=None):
        return np.array([a.get_coord() for a in self.get_atoms(
            include_hetatms=include_hetatms, exclude_atoms=exclude_atoms)])

    def orient_to_pai(self, flip_axis=(0.2, 0.2, 0.2)):
        coords = PCA(n_components = 3).fit_transform(self.get_coords())
        coords = flip_around_axis(coords, axis=flip_axis)
        self.update_coords(coords)

    def rotate(self, num=1):
        """Rotate structure in randomly in place"""
        for r in range(num):
            M, theta, phi, z = rotation_matrix(random=True)
            self.shift_coords_to_origin()
            coords = np.dot(self.get_coords(), M)
            self.update_coords(coords)
            self.shift_coords_to_volume_center()
            self.set_voxel_size(self.voxel_size)
            yield r, theta, phi, z

    def update_coords(self, coords):
        for atom, coord in zip(self.structure.get_atoms(), coords):
            atom.set_coord(coord)
        self.mean_coord = None
        self.mean_coord_updated = False

    def create_full_volume(self, input_shape=(96, 96, 96)):
        truth_grid = np.zeros(list(input_shape)+[1])
        for atom in self.get_atoms():
            grid = self.get_grid_coord(atom, vsize=input_shape[0])
            truth_grid[grid[0], grid[1], grid[2], 0] = 1
        return truth_grid

    def get_residue_from_resseq(self, resseq, model=0, chain=None):
        chain = chain or self.chain
        try:
            return self.structure[model][chain][resseq]
        except KeyError as e:
            for r in self.structure[model][chain]:
                if r.get_id()[1] == resseq:
                    return r
            else:
                return None

    def align_seq_to_struc(self, *seq_num, **kwds):
        return_residue=kwds.get("return_residue", False)

        if return_residue:
            mapped_residues = [self.get_residue_from_resseq(pdbnum) \
                for current_resi, resn, pdbnum, ncbi in map_residues(self.pdb, self.chain, seq_num)]

            if mapped_residues.count(None) > len(mapped_residues)/2:
                raise InvalidPDB("Binding Site ({:4}%) missing from structure ({}.{})".format(mapped_residues.count(None)/float(len(mapped_residues)), self.pdb, self.chain))

            mapped_residues = [r for r in mapped_residues if r is not None]
        else:
            mapped_residues = [pdb for current_resi, resn, pdb, ncbi in map_residues(self.pdb, self.chain, seq_num)]

        return mapped_residues

    @staticmethod
    def get_feature_names(only_aa=False, only_atom=False):
        if only_atom:
            return ["C", "N", "O", "S", "Unk_element"]
        if only_aa:
            return PDB.Polypeptide.aa3
        feature_names = ["C", "CT", "CA", "N", "N2", "N3", "NA", "O", "O2", "OH", "S", "SH", "Unk_atom", "C", "N", "O", "S", "Unk_element", "vdw_volume", "charge", "neg_charge", "pos_charge", "neutral_charge", "cx", "is_concave", "is_convex", "is_concave_and_convex", "hydrophobicity", "is_hydrophbic", "is_hydrophilic", "atom_asa", "residue_asa", "residue_buried", "residue_exposed"]
        feature_names += PDB.Polypeptide.aa3
        feature_names += ["Unk_residue", "is_helix", "is_sheet", "Unk_SS"]
        return feature_names

    def get_features_per_atom(residue_list):
        """Get features for eah atom, but not organized in grid"""
        features = [self.get_features_for_atom(a) for r in residue_list for a in r]
        return features

    def get_features(self, input_shape=(96, 96, 96), residue_list=None, batch_index=None, return_data=True, return_truth=True, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, expand_atom=False, include_full_protein=False, undersample=False, voxel_size=1.0):
        if self.course_grained:
            return self.map_residues_to_voxel_space(
                binding_site_residues=residue_list,
                include_full_protein=include_full_protein,
                only_aa=only_aa,
                non_geom_features=non_geom_features,
                undersample=undersample
            )
        return self.map_atoms_to_voxel_space(
            expand_atom=expand_atom,
            binding_site_residues=residue_list,
            include_full_protein=include_full_protein,
            only_aa=only_aa,
            only_atom=only_atom,
            use_deepsite_features=use_deepsite_features,
            non_geom_features=non_geom_features,
            undersample=undersample)

    def get_features_in_grid(self, input_shape=(96, 96, 96), sparse=True, residue_list=None, batch_index=None, return_data=True, return_truth=True, only_aa=False, expand_atom=True, include_full_protein=False, voxel_size=0.5):
        if return_truth and residue_list is None:
            raise RuntimeError("Can only output truth if given residue_list")

        if not sparse:
            raise RuntimeError("Dense gris not supported; Too much memory")

        exclude_atoms = None #["N", "C", "O"] #Ignore backbone

        if residue_list is not None:
            if not include_full_protein:
                atoms = sorted((a for r in residue_list for a in r), key=lambda a: a.get_serial_number())
                atoms = list(self.filter_atoms(atoms, exclude_atoms=exclude_atoms))
                nAtoms = len(atoms)
                binding_site_atoms = [a.get_serial_number() for a in atoms]
            else:
                atoms = list(self.get_atoms(include_hetatms=True, exclude_atoms=exclude_atoms))
                nAtoms = len(atoms)
                binding_site_atoms = [a.get_serial_number() for r in residue_list for a in r]
        else:
            atoms = list(self.get_atoms(include_hetatms=True, exclude_atoms=exclude_atoms))
            nAtoms = len(atoms)
            truth_atoms = None
            binding_site_atoms = []

        if only_atom:
            nFeatures = 5
        elif only_aa:
            nFeatures = 21
        else:
            nFeatures = Structure.nFeatures

        if expand_atom:
            elems = Counter()
            for atom in atoms:
                elems[atom.element.title()] += 1

            num_voxels = sum([atom_spheres[elem].shape[0]*count for elem, count in elems.items()])
        else:
            num_voxels = len(atoms)

        indices = np.zeros((num_voxels, 3))
        data = np.zeros((num_voxels, nFeatures))
        truth = np.zeros((num_voxels, 1))

        last_index = 0

        for i, atom in enumerate(atoms):
            grid = self.get_grid_coord(atom, vsize=input_shape[0])
            if grid[0] >= input_shape[0] or grid[1] >= input_shape[1] or grid[2] >= input_shape[2]:
                continue

            features = self.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom)
            truth_val = int(atom.get_serial_number() in binding_site_atoms)

            if expand_atom:
                sphere_voxels = self.expand_atom(atom)
                next_index = last_index+sphere_voxels.shape[0]
                indices[last_index:next_index] = sphere_voxels
                data[last_index:next_index] = features
                truth[last_index:next_index] = truth_val
                last_index = next_index
            else:
                indices[i] = grid
                data[i] = features
                truth[i] = truth

        if return_data and return_truth:
            return indices, data, truth
        elif return_data:
            return indices, data
        else:
            return indices, truth

    def get_neighbors(self, atom, radius=2.):
        if self.kdtree is None:
            self.kdtree = NeighborSearch(self.get_atoms())

        neighbors = self.kdtree.search(atom.get_coord(), radius)
        return neighbors

    def map_atoms_to_voxel_space(self, expand_atom=False, binding_site_residues=None, include_full_protein=False, only_aa=False, only_atom=False, use_deepsite_features=False, non_geom_features=False, undersample=False, only_surface=True):
        """Map atoms to sparse voxel space.

        Parameters
        ----------
        expand_atom : boolean
            If true, atoms will be converted into spheres with their Van der walls
            radii. The features for the atom are copied into all voxels that contain
            the atom and features from overlapping atoms are summed. If false, an atom
            will only occupy one voxel, where overlapping features for overlapping atoms
            are summed.
        binding_site_residues : list of Bio.PDB.Residue objects or None
            If a binding is known, add the list of Bio.PDB.Residue objects, usually
            obtained by Structure.align_seq_to_struc()
        include_full_protein : boolean
            If true, all atoms from the protein are used. Else, just the atoms from the
            defined binding site. Only makes sense if binding_site_residues is not None
        Returns
        -------
        indices : np.array((nVoxels,3))
        data : np.array((nVoxels,nFeatures))
        """
        if binding_site_residues is not None:
            if not include_full_protein:
                atoms = sorted((a for r in binding_site_residues for a in r), key=lambda a: a.get_serial_number())
                atoms = list(self.filter_atoms(atoms))
                binding_site_atoms = [a.get_serial_number() for a in atoms]
                non_binding_site_atoms = []
            else:
                atoms = list(self.get_atoms(include_hetatms=True))
                nAtoms = len(atoms)
                binding_site_atoms = [a.get_serial_number() for r in binding_site_residues for a in r]

                if undersample:
                    non_binding_site_residues = []
                    for r in self.structure.get_residues():
                        if r in binding_site_residues: continue

                        if only_surface and bool(self.precalc_features[r.get_list()[0].serial_number-1][33]):
                            continue

                        non_binding_site_residues.append(r.get_id()[1])

                    try:
                        non_binding_site_residues = np.random.choice(non_binding_site_residues, len(binding_site_residues))
                        non_binding_site_atoms = []
                        for r in non_binding_site_residues:
                            try:
                                r = self.structure[0][self.chain][r]
                            except KeyError:
                                continue
                            for a in r:
                                non_binding_site_atoms.append(a.get_serial_number())
                    except ValueError as e:
                        print(e)
                        #Might give over balance
                        non_binding_site_atoms = []
                else:
                    non_binding_site_atoms = []

                    # non_binding_site_atoms = [a.get_serial_number() for a in atoms if a.get_serial_number() not in binding_site_atoms]
                    # try:
                    #     non_binding_site_atoms = np.random.choice(non_binding_site_atoms, len(binding_site_atoms))
                    # except ValueError:
                    #     #Might give over balance
                    #     non_binding_site_atoms = []

                atoms = list(atoms)

        else:
            atoms = list(self.get_atoms(include_hetatms=True))
            nAtoms = len(atoms)
            binding_site_atoms = []
            non_binding_site_atoms = []


        nFeatures = Structure.number_of_features(
            only_aa=only_aa,
            only_atom=only_atom,
            non_geom_features=non_geom_features,
            use_deepsite_features=use_deepsite_features,
            course_grained=False)

        data_voxels = defaultdict(lambda: np.zeros(nFeatures))
        truth_voxels = {}

        skipped = 0
        for atom in atoms:
            truth = atom.get_serial_number() in binding_site_atoms

            if not truth and undersample and atom.get_serial_number() not in non_binding_site_atoms:
                skipped += 1
                continue

            if only_surface:
                features, is_buried = self.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features, warn_if_buried=True)
                if not truth and is_buried:
                    skipped += 1
                    continue
            else:
                features = self.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features)

            truth = np.array([float(truth)])

            for atom_grid in self.get_grid_coords_for_atom_by_kdtree(atom):
                atom_grid = tuple(atom_grid.tolist())
                try:
                    data_voxels[atom_grid] = np.maximum(features, data_voxels[atom_grid])
                except ValueError:
                    print(nFeatures, data_voxels[atom_grid].shape, features.shape)
                    raise
                truth_voxels[atom_grid] = np.array([0.,1.]) if truth else np.array([1.,0.])

        try:
            if binding_site_residues is None:
                return np.array(list(data_voxels)), np.array(list(data_voxels.values()))
            else:
                truth = np.array([truth_voxels[grid] for grid in data_voxels])
                return np.array(list(data_voxels)), np.array(list(data_voxels.values())), truth
        except Exception as e:
            print(e)
            raise

    def map_residues_to_voxel_space(self, binding_site_residues=None, include_full_protein=False, non_geom_features=True, only_aa=False, use_deepsite_features=False, undersample=False):
        if binding_site_residues is not None:
            if not include_full_protein:
                residues = binding_site_residues
                binding_site_residues = [r.get_id()[1] for r in residues]
            else:
                residues = self.structure.get_residues()
                binding_site_residues = [r.get_id()[1] for r in binding_site_residues]

                if undersample:
                    non_binding_site_residues = [r.get_id()[1] for r in self.structure.get_residues() if r not in binding_site_residues]
                    try:
                        non_binding_site_residues = np.random.choice(non_binding_site_residues, len(binding_site_residues))
                    except ValueError as e:
                        print(e)
                        #Might give over balance
                        non_binding_site_residues = []
        else:
            residues = self.structure.get_residues()
            binding_site_residues = []

        nFeatures = Structure.number_of_features(
            only_aa=only_aa,
            non_geom_features=non_geom_features,
            course_grained=True,
            use_deepsite_features=use_deepsite_features)

        data_voxels = defaultdict(lambda: np.zeros(nFeatures))
        truth_voxels = {}

        residues = list(residues)

        for residue in residues:
            #assert not residue.get_id()[1] in binding_site_residues
            truth = residue.get_id()[1] in binding_site_residues
            if not truth and undersample and residue.get_id()[1] not in non_binding_site_residues:
                continue

            truth = np.array([int(truth)])

            try:
                features = self.get_features_for_residue(residue, only_aa=only_aa, non_geom_features=non_geom_features, use_deepsite_features=False)
            except Exception as e:
                print(e)
                raise
            for residue_grid in self.get_grid_coords_for_residue_by_kdtree(residue):
                residue_grid = tuple(residue_grid.tolist())
                try:
                    data_voxels[residue_grid] = np.maximum(features, data_voxels[residue_grid])
                except ValueError:
                    print(nFeatures, data_voxels[residue_grid].shape, features.shape)
                    raise
                truth_voxels[residue_grid] = truth

        if binding_site_residues is None:
            return np.array(list(data_voxels)), np.array(list(data_voxels.values()))
        else:
            truth = np.array([truth_voxels[grid] for grid in data_voxels.keys()])
            return np.array(list(data_voxels)), np.array(list(data_voxels.values())), truth

    def voxel_set_insection_and_difference(self, atom1, atom2):
        A = self.atom_spheres[atom1.get_serial_number()]
        B = self.atom_spheres[atom2.get_serial_number()]

        nrows, ncols = A.shape
        dtype={'names':['f{}'.format(i) for i in range(ncols)],
               'formats':ncols * [A.dtype]}

        intersection = np.intersect1d(A.view(dtype), B.view(dtype))
        intersection = intersection.view(A.dtype).reshape(-1, ncols)

        onlyA = np.setdiff1d(A.view(dtype), B.view(dtype))
        onlyA = onlyA.view(A.dtype).reshape(-1, ncols)

        onlyB = np.setdiff1d(B.view(dtype), A.view(dtype))
        onlyB = onlyA.view(A.dtype).reshape(-1, ncols)

        return intersection, onlyA, onlyB

    def get_features_for_atom(self, atom, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, warn_if_buried=False, preload=True):
        """Calculate FEATUREs"""
        if isinstance(atom, PDB.Atom.DisorderedAtom):
            #All altlocs have been removed so onlt one remains
            atom = atom.disordered_get_list()[0]

        if preload and not self.force_feature_calculation and self.precalc_features is not None:
            try:
                features = self.precalc_features[atom.serial_number-1]
                is_buried = bool(features[35]) #Residue asa #[self.precalc_features[a.serial_number-1][31] for a in atom.get_parent()]
                # if asa > 0.0:
                #     asa /= surface_areas.get(atom.element.title(), 1.0)
                #     is_buried = asa <= 0.2
                # else:
                #     is_buried = False

                if use_deepsite_features:
                    feats = np.concatenate((
                        features[61:67],
                        features[20:22]))
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                if only_atom:
                    feats = features[13:18]
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                elif only_aa:
                    feats = features[37:58]
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                elif non_geom_features:
                    feats = np.concatenate((
                        features[13:18],
                        features[19:23],
                        features[27:30], #1]))
                        np.array([float(is_buried)])))
                    if warn_if_buried:
                        return feats, is_buried
                    else:
                        return feats
                else:
                    if warn_if_buried:
                        return features, is_buried
                    else:
                        return features
            except ValueError as e:
                # print e
                # pass
                raise

        if use_deepsite_features:
            if warn_if_buried:
                return self.get_deepsite_features(atom), self.get_accessible_surface_area(atom)[-2]
            else:
                return self.get_deepsite_features(atom)
        elif only_atom:
            if warn_if_buried:
                return self.get_element_type(atom), self.get_accessible_surface_area(atom)[-2]
            else:
                return self.get_element_type(atom)
        elif only_aa:
            return self.get_residue(atom)
        elif non_geom_features:
            features = np.zeros(13)
            features[0:5] = self.get_element_type(atom)
            features[5:9] = self.get_charge(atom)
            features[9:13] = self.get_hydrophobicity(atom)
            if warn_if_buried:
                return features, self.get_accessible_surface_area(atom)[-2]
            else:
                return features
        else:
            print "getting full features"
            features = np.zeros(self.nFeatures)

            #atom_type
            features[0:13]  = self.get_atom_type(atom)
            features[13:18] = self.get_element_type(atom)
            features[18:19] = self.get_vdw(atom)
            features[19:23] = self.get_charge(atom)
            features[23:27] = self.get_concavity(atom)
            features[27:31] = self.get_hydrophobicity(atom)
            features[31:37] = self.get_accessible_surface_area(atom)
            features[37:58] = self.get_residue(atom)
            features[58:61] = self.get_ss(atom)
            features[61:67] = self.get_deepsite_features(atom, calc_charge=False)

            if self.precalc_features is not None:
                self.precalc_features[atom.serial_number-1] = features

            if warn_if_buried:
                return features, bool(features[35])
            else:
                return features

    def get_features_for_residue(self, residue, only_aa=False, non_geom_features=False, use_deepsite_features=False, preload=True):
        """Calculate FEATUREs"""
        if preload and self.precalc_features is not None:
            print("Using precalc features")
            try:
                features = self.precalc_features[residue.get_id()[1]-1]
                if non_geom_features:
                    return np.concatenate((
                        features[15:36],
                        features[0:4],
                        features[8:12],
                        ))
                elif only_aa:
                    return features[15:36]
                else:
                    return features[:self.nFeatures]
            except ValueError:
                pass

        if non_geom_features:
            features = np.zeros(29)
            features[0:21] = self.get_residue(residue)
            features[21:25] = self.get_charge(residue)
            features[25:29] = self.get_hydrophobicity(residue)
            return features
        elif only_aa:
            return self.get_residue(residue)
        else:
            features = np.zeros(39)

            features[0:4] = self.get_charge(residue)
            features[4:8] = self.get_concavity(residue)
            features[8:12] = self.get_hydrophobicity(residue)
            features[12:15] = self.get_accessible_surface_area(residue)
            features[15:36] = self.get_residue(residue)
            features[36:39]  = self.get_ss(residue)

            if self.precalc_features is not None:
                self.precalc_features[residue.get_id()[1]-1] = features

            return features

    def get_atom_type(self, atom):
        """
        ATOM_TYPE_IS_C
        ATOM_TYPE_IS_CT
        ATOM_TYPE_IS_CA
        ATOM_TYPE_IS_N
        ATOM_TYPE_IS_N2
        ATOM_TYPE_IS_N3
        ATOM_TYPE_IS_NA
        ATOM_TYPE_IS_O
        ATOM_TYPE_IS_O2
        ATOM_TYPE_IS_OH
        ATOM_TYPE_IS_S
        ATOM_TYPE_IS_SH
        ATOM_TYPE_IS_OTHER"""
        atom_types = ["C", "CT", "CA", "N", "N2", "N3", "NA", "O", "O2", "OH", "S", "SH"]
        atom_type = np.zeros(13)
        try:
            index = atom_types.index(atom.get_name().strip())
            atom_type[index] = 1.
        except ValueError:
            atom_type[12] = 1.
        return atom_type

    def get_element_type(self, atom):
        """ELEMENT_IS_ANY
        ELEMENT_IS_C
        ELEMENT_IS_N
        ELEMENT_IS_O
        ELEMENT_IS_S
        ELEMENT_IS_OTHER"""
        elems = "CNOS"
        elem_type = np.zeros(5)
        try:
            index = elems.index(atom.element)
            elem_type[index] = 1.
        except ValueError:
            elem_type[4] = 1.
        return elem_type

    def get_vdw(self, atom):
        return np.array([vdw_radii.get(atom.element.title(), 2.0)])

    def get_charge(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            charge_value = np.sum([self.get_charge(a)[0] for a in residue])
            charge = np.zeros(4)
            charge[0] = charge_value
            charge[1] = float(charge_value < 0)
            charge[2] = float(charge_value > 0)
            charge[3] = float(charge_value == 0)
            return charge
        else:
            raise RuntimeError("Input must be Atom or Residue: {}".format(type(atom_or_residue)))

        pqr = self._get_pqr()
        atom_id = atom.get_full_id()[3:5]

        if atom_id[1][1] != " ":
            #pdb2pqr removes alternate conformations and only uses the first
            atom_id = (atom_id[0], (atom_id[1][0], " "))
        charge_value = pqr.get(atom_id, np.NaN)

        charge = np.zeros(4)
        charge[0] = charge_value
        print "charge", charge_value, float(charge_value < 0)
        charge[1] = float(charge_value < 0)
        charge[2] = float(charge_value > 0)
        charge[3] = float(charge_value == 0)
        return charge

    def get_concavity(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            concavity_value = np.mean([self.get_concavity(a)[0] for a in residue])
            #concave, convex, or both
            concavity = np.zeros(4)
            concavity[0] = concavity_value
            concavity[1] = float(concavity_value <= 2)
            concavity[2] = float(concavity_value > 5)
            concavity[3] = float(2 < concavity_value <= 5)
            return concavity
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        cx = self._get_cx()
        concavity_value = cx.get(atom.serial_number, np.NaN)
        #concave, convex, or both
        concavity = np.zeros(4)
        concavity[0] = concavity_value
        concavity[1] = float(concavity_value <= 2)
        concavity[2] = float(concavity_value > 5)
        concavity[3] = float(2 < concavity_value <= 5)
        return concavity

    def get_hydrophobicity(self, atom_or_residue, scale="kd"):
        assert scale in hydrophobicity_scales.keys()
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        hydrophobicity = np.zeros(4)
        try:
            resname = PDB.Polypeptide.three_to_one(residue.get_resname())

            hydrophobicity_value = hydrophobicity_scales[scale][resname]
            hydrophobicity[0] = hydrophobicity_value
            hydrophobicity[1] = float(hydrophobicity_value < 0)
            hydrophobicity[2] = float(hydrophobicity_value > 0)
            hydrophobicity[3] = float(hydrophobicity_value == 0)
        except KeyError:
            hydrophobicity[0] = np.NaN
            hydrophobicity[1] = np.NaN
            hydrophobicity[2] = np.NaN
            hydrophobicity[3] = np.NaN

        return hydrophobicity

    def get_accessible_surface_area(self, atom_or_residue):
        """Returns the ASA value from freesasa (if inout is Atom) and the DSSP
        value (if input is Atom or Residue)

        Returns
        -------
        If input is residue a 3-vector is returned, otherwise a 4-vector is returned
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if is_atom:
            total_area = surface_areas.get(atom.element.title(), 1.0)
            try:
                sasa, sasa_struct = self._get_sasa()
                selection = "sele, chain {} and resi {} and name {}".format(self.chain, atom.get_parent().get_id()[1], atom.get_id()[0])
                with silence_stdout(), silence_stderr():
                    selections = freesasa.selectArea([selection], sasa_struct, sasa)
                    atom_area = selections["sele"]
                fraction = atom_area/total_area
            except (KeyError, AssertionError, AttributeError, TypeError):
                raise
                atom_area = np.NaN
                fraction = np.NaN

        dssp = self._get_dssp()
        try:
            residue_area = dssp[residue.get_full_id()[2:]][3]
        except (KeyError, AssertionError, AttributeError, TypeError):
            try:
                #Remove HETATMs
                residue_area = dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][3]
                if residue_area == "NA":
                    residue_area = np.NaN
            except (KeyError, AssertionError, AttributeError, TypeError):
                residue_area = np.NaN

        if is_atom:
            asa = np.zeros(6)
            asa[0] = atom_area
            asa[1] = float(fraction <= 0.2) #buried
            asa[2] = float(fraction > 0.2) #exposed
            asa[3] = residue_area
            asa[4] = float(residue_area < 0.2)
            asa[5] = float(residue_area >= 0.2)
        else:
            asa = np.zeros(3)
            asa[0] = residue_area
            asa[1] = float(residue_area < 0.2)
            asa[2] = float(residue_area >= 0.2)

        return asa

    def get_residue(self, atom_or_residue):
        """
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        residues = [0.]*(len(PDB.Polypeptide.aa1)+1) #one hot residue

        try:
            residues[PDB.Polypeptide.three_to_index(residue.get_resname())] = 1.
        except (ValueError, KeyError) as e:
            residues[-1] = 1.
        return residues

    def get_ss(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        dssp = self._get_dssp()
        try:
            atom_ss = dssp[residue.get_full_id()[2:]][2]
        except (KeyError, AssertionError, AttributeError, TypeError):
            try:
                #Remove HETATMs
                atom_ss = dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][2]
            except (KeyError, AssertionError, AttributeError, TypeError):
                atom_ss = "X"

        ss = np.zeros(3)
        ss[0] = float(atom_ss in "GHIT")
        ss[1] = float(atom_ss in "BE")
        ss[2] = float(atom_ss not in "GHITBE")
        return ss

    def get_deepsite_features(self, atom, calc_charge=True):
        """Use DeepSites rules for autodock atom types
        """
        element, is_hbond_donor = self.get_autodock_features(atom)
        features = np.zeros(8 if calc_charge else 6, dtype=bool)

        #hydrophobic
        features[0] = (element == 'C') | (element == 'A')

        #aromatic
        features[1] = element == 'A'

        #hbond_acceptor
        features[2] = (element == 'NA') | (element == 'NS') | (element == 'OA') | \
                      (element == 'OS') | (element == 'SA')

        #hbond_donor
        features[3] = is_hbond_donor

        #metal
        features[4] = (element == 'MG') | (element == 'ZN') | (element == 'MN') | \
                      (element == 'CA') | (element == 'FE')

        #occupancies / excluded volume
        features[5] = (element != 'H') & (element != 'HS') & (element != 'HD')

        if calc_charge:
            #positive_ionizable
            charge = self.get_charge(atom)[0]
            features[6] = charge > 0

            #negative_ionizable
            features[7] = charge < 0

        return features.astype(float)

    def get_grid_coords_for_atom_by_kdtree(self, atom, k=4):
        dist = self.get_vdw(atom)[0]
        neighbors = self.voxel_tree.query_ball_point(atom.coord, r=dist)
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_grid_coords_for_residue_by_kdtree(self, residue):
        dist = vdw_aa_radii.get(residue.get_resname(), 3.2)
        center = np.mean([a.get_coord() for a in residue], axis=0)
        neighbors = self.voxel_tree.query_ball_point(center, r=dist)
        return [self.voxel_tree.data[idx] for idx in neighbors]

    def get_grid_coord_for_residue(self, residue, round=False):
        coord = np.mean([a.get_coord() for a in residue], axis=0) #residue["CA"].get_coord() #
        # coord += [self.volume/2]*3
        if round:
            coord = np.around(coord)
        voxel = np.array([
            np.digitize([coord[0]], self.voxels_x)[0]-1,
            np.digitize([coord[1]], self.voxels_y)[0]-1,
            np.digitize([coord[2]], self.voxels_z)[0]-1,
            ])
        return voxel, coord

    def set_volume(self, volume):
        self.volume = volume

    def set_voxel_size(self, voxel_size=None):
        if not self.course_grained:
            self.voxel_size = voxel_size or 1.0
        else:
            self.voxel_size = 10.0

        coords = self.get_coords()
        min_coord = np.floor(np.min(coords, axis=0))-5
        max_coord = np.ceil(np.max(coords, axis=0))+5
        extent_x = np.arange(min_coord[0], max_coord[0], self.voxel_size)
        extent_y = np.arange(min_coord[1], max_coord[1], self.voxel_size)
        extent_z = np.arange(min_coord[2], max_coord[2], self.voxel_size)
        mx, my, mz = np.meshgrid(extent_x, extent_y, extent_z)
        self.voxel_tree = spatial.cKDTree(list(zip(mx.ravel(), my.ravel(), mz.ravel())))

    def get_atoms_from_grids(self, grids, vsize=96, max_radius=40):
        for atom in self.get_atoms():
            coord = self.get_grid_coord(atom, vsize, max_radius)
            if coord.tolist() in grids:
                yield atom

def flip_around_axis(coords, axis = (0.2, 0.2, 0.2)):
    'Flips coordinates randomly w.r.t. each axis with its associated probability'
    for col in xrange(3):
        if np.random.binomial(1, axis[col]):
            coords[:,col] = np.negative(coords[:,col])
    return coords

def rotation_matrix(random = False, theta = 0, phi = 0, z = 0):
    'Creates a rotation matrix'
    # Adapted from: http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
    # Initialization
    if random == True:
        randnums = np.random.uniform(size=(3,))
        theta, phi, z = randnums
    theta = theta * 2.0*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0  # For magnitude of pole deflection.
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    st = np.sin(theta)
    ct = np.cos(theta)
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    M = (np.outer(V, V) - np.eye(3)).dot(R)

    return M, theta, phi, z

if __name__ == "__main__":
  import sys
  assert len(sys.argv) > 1
  structure = Structure(sys.argv[1], panfs=False).extract_chain("A")
  num_atoms = sum(1 for a in structure.structure.get_atoms() if a.get_parent().get_id()[0] == " ")
  for atom in structure.get_atoms():
    features = structure.get_features_for_atom(atom)
    print(atom.get_full_id(), atom.serial_number, "of", num_atoms, features, len(features))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
