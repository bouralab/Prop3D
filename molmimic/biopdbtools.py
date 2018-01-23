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
from sklearn.decomposition import PCA
from scipy.interpolate import RegularGridInterpolator
from Bio import PDB
from Bio import SeqIO
from Bio.PDB.NeighborSearch import NeighborSearch
try:
    from pdb2pqr import mainCommand
except ImportError:
    mainCommand = None
try:
    import freesasa
except ImportError:
    freesasa = None

import h5py

from util import silence_stdout, silence_stderr
from map_residues import map_residues

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


def make_sphere(atom, radius, voxel_size=0.25):
    """Make sphere with given radius centered at the origin"""
    assert radius > 0 and voxel_size > 0

    sphere_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "{}_{}.npy".format(atom, voxel_size))
    try:
        sphere_coords = np.load(sphere_file)
    except IOError:
        weight = 1./voxel_size
        print "using weight", weight
        radius *= max(1, weight)
        shape = 2*radius

        y, x, z = np.ogrid[-radius:radius, -radius:radius, -radius:radius]
        sphere_coords = x*x + y*y + z*z < radius*radius
        sphere_coords = np.array(np.where(sphere_coords==True)).T.astype("float64")
        sphere_coords -= np.mean(sphere_coords, axis=0)
        np.save(sphere_file, sphere_coords)
    return sphere_coords

atom_spheres = defaultdict()
def make_atom_spheres(voxel_size=0.25, default_radius=2):
    global atom_spheres
    #atom_spheres = {atom:make_sphere(vdw) for atom, vdw in vdw_radii.iteritems()}
    atom_spheres = {}
    for atom, vdw in vdw_radii.iteritems():
        atom_spheres[atom] = make_sphere(atom, vdw, voxel_size)
    default_sphere = make_sphere("default", default_radius, voxel_size)
    atom_spheres = defaultdict(lambda:default_sphere, atom_spheres)
    return atom_spheres

maxASA = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0, "K": 201.0, "L": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}

#REMOVE HARD CODED PATH
obsolete_pdbs = pd.read_table(os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "obsolete.dat"), delim_whitespace=True, skiprows=1, names=["status", "data", "old", "new"], header=None,)

pdbparser = PDB.PDBParser()
writer = PDB.PDBIO()

_hydrogen = re.compile("[123 ]*H.*")

class InvalidPDB(RuntimeError):
    pass

class SelectChain(PDB.Select):
    """ Only accept the specified chain and remove hydruogens and hetatms when saving. """
    def __init__(self, chain):
        self.chain = chain

    def accept_model(self, model):
        # model - only keep model 0
        return model.get_id() == 0

    def accept_chain(self, chain):
        return chain.get_id() == self.chain

    def accept_residue(self, residue):
        """Remove HETATMS"""
        hetatm_flag, resseq, icode = residue.get_id()
        return hetatm_flag == " " or residue.get_resname() != "HOH"

    def accept_atom(self, atom):
        """Remove hydrogens"""
        name = atom.get_id()
        return not _hydrogen.match(name) and atom.altloc == " "

def extract_chain(structure, chain, pdb):
    chain = chain.upper()
    chain_pdb = StringIO()
    writer.set_structure(structure)
    writer.save(chain_pdb, select=SelectChain(chain))
    chain_pdb.seek(0)
    new_structure = pdbparser.get_structure(pdb, chain_pdb)
    return new_structure

class Structure(object):
    nFeatures = 59

    def __init__(self, pdb, chain, path=None, id=None, snapshot=True, course_grained=True, input_format="pdb", force_feature_calculation=False):
        if path is None and not os.path.isfile(pdb) and len(pdb) == 4:
            path = "{}/pdb/{}/pdb{}.ent.gz".format(os.environ.get("PDB_SNAPSHOT", "/pdb"), pdb[1:3].lower(), pdb.lower())

            if snapshot:
                if not os.path.isfile(path):
                    path, input_format = download_pdb(pdb)

                if not os.path.isfile(path):
                    raise InvalidPDB()
                    #Obsolete pdb
                    try:
                        _, old_new = obsolete_pdbs.loc[obsolete_pdbs["old"] == pdb.upper()].iterrows().next()
                        if old_new["new"]:
                            pdb = old_new["new"]
                            path = "{}/pdb/{}/pdb{}.ent.gz".format(os.environ.get("PDB_SNAPSHOT", "/pdb"), pdb[1:3].lower(), pdb.lower())
                            if not os.path.isfile(path):
                                raise StopIteration
                        else:
                            raise StopIteration
                    except StopIteration:
                        path, input_format = download_pdb(pdb)
            else:
                path, input_format = download_pdb(pdb)
            pdb = pdb.lower()
        elif path is None:
            path = pdb

        if input_format == "pdb":
            parser = pdbparser
        elif input_format == "mmcif":
            parser = PDB.FastMMCIFParser()
        elif input_format == "mmtf":
            parser = PDB.MMTFParser()
        else:
            raise RuntimeError("Invalid PDB parser (pdb, mmcif, mmtf)")

        try:
            if isinstance(path, str) and os.path.isfile(path):
                if path.endswith(".gz"):
                    with gzip.open(path, 'rb') as f:
                        self.structure = parser.get_structure(pdb, f)
                    path = path[:-2]
                else:
                    self.structure = parser.get_structure(pdb, path)
            elif isinstance(path, Iterable):
                self.structure = parser.get_structure(pdb, path)
            else:
                raise InvalidPDB("Invalid PDB id or file: {} (path={})".format(pdb, path))
        except KeyError:
            #Invalid mmcif file
            raise InvalidPDB("Invalid PDB id or file: {} (path={})".format(pdb, path))

        self.chain = chain.split("_", 1)[0]
        self.structure = extract_chain(self.structure, self.chain, pdb)

        # self.nResidues = sum(1 for _ in self.structure.get_residues())
        # if self.nResidues > 500:
        #     raise InvalidPDB("{}.{} (residue count > 500)".format(pdb, self.chain))

        only_chain = self.structure[0].get_chains().next()
        if self.chain != only_chain.get_id():
            #Reset chain if Bio.PDB changes name after extraction
            self.chain = only_chain.get_id()

        self.pdb = pdb
        self.path = path
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
        self.volume = 256
        self.id = "{}_{}{}".format(self.pdb, self.chain, "_{}".format(id) if id else "")
        self.course_grained = course_grained
        self.nFeatures = 59 if not course_grained else 36

        #Creates self.voxels
        self.set_voxel_size(0.5)

        precalc_features_path = os.environ.get("MOLMIMIC_FEATURES", os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "features2"))
        precalc_features_path = os.path.join(precalc_features_path, "{}.h5".format(self.id))

        if force_feature_calculation or not os.path.isfile(precalc_features_path):
            last_atom = None
            for last_atom in self.structure.get_atoms():
                pass
            if last_atom is None:
                raise InvalidPDB(pdb)
            self.feature_file = h5py.File(precalc_features_path, 'w')
            self.precalc_features = self.feature_file.create_dataset('features', dtype='f', shape=(last_atom.serial_number, Structure.nFeatures), fillvalue=0., compression='gzip', compression_opts=9)
        else:
            try:
                self.precalc_features = h5py.File(precalc_features_path, 'r')["features"]
                self.feature_file = None
            except IOError:
                last_atom = None
                for last_atom in self.structure.get_atoms():
                    pass
                if last_atom is None:
                    raise InvalidPDB(pdb)
                self.feature_file = h5py.File(precalc_features_path, 'w')
                self.precalc_features = self.feature_file.create_dataset('features', dtype='f', shape=(last_atom.serial_number, Structure.nFeatures), fillvalue=0., compression='gzip', compression_opts=9)


        # try:
        #     aa = SeqIO.index(os.path.join(os.environ.get("PDB_SNAPSHOT", "/pdb")), "sequences", "pdb_seqres.txt")
        #     self.aa = str(aa["{}_{}".format(pdb, chain[0].upper())].seq)
        # except KeyError:
        #     self.aa = None
        #     print "Canot read seq file:", "{}/sequences/pdb_seqres.txt".format(os.environ.get("PDB_SNAPSHOT", "/pdb"))

    @staticmethod
    def features_from_string(pdb, chain, resi=None, id=None, input_shape=(96,96,96), batch_index=None, only_aa=False, only_atom=False, course_grained=False, grid=True, return_data=True, return_truth=True, rotate=True, force_feature_calculation=False, expand_atom=True, include_full_protein=False):
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
        s = Structure(
            pdb,
            chain,
            id=id,
            course_grained=course_grained,
            force_feature_calculation=force_feature_calculation)
        s.orient_to_pai()
        if resi:
            binding_site = s.align_seq_to_struc(resi, return_residue=True)
        else:
            binding_site = None

        if grid:
            if rotate:
                s.rotate(1).next()
            features = s.get_features(
                input_shape=input_shape,
                residue_list=binding_site,
                batch_index=batch_index,
                only_aa=only_aa,
                only_atom=only_atom,
                expand_atom=expand_atom,
                include_full_protein=include_full_protein)
        else:
            features = s.get_features_per_atom(binding_site)

        if s.feature_file:
            s.feature_file.close()

        return features


    def get_atoms(self, include_hetatms=False, exclude_atoms=None):
        # for a in self.structure.get_atoms():
        #     hetflag, resseq, icode = a.get_parent().get_id()
        #     if not include_hetatms and hetflag is not ' ':
        #         continue
        #     if exclude_atoms is not None and atom.get_name().strip() in exclude_atoms:
        #         continue
        #     yield a
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

        writer.set_structure(self.structure)
        writer.save(path)

        if file_like:
            path.seek(0)

        if lines:
            path = path.read()

        return path

    def _get_dssp(self):
        if self.dssp is None:
            pdbfd, tmp_pdb_path = tempfile.mkstemp()
            with os.fdopen(pdbfd, 'w') as tmp:
                writer.set_structure(self.structure)
                writer.save(tmp)

            try:
                self.dssp = PDB.DSSP(self.structure[0], tmp_pdb_path, dssp='mkdssp')
            except NameError:
                self.dssp = None

            os.remove(tmp_pdb_path)

        return self.dssp

    def _get_sasa(self):
        if freesasa is None:
            print "SASA not installed! SASA will be 0"
            return None, None
        if self.sasa is None:
            with silence_stdout(), silence_stderr():
                self.sasa_struct = freesasa.structureFromBioPDB(self.structure)
                self.sasa = freesasa.calc(self.sasa_struct)

        return self.sasa, self.sasa_struct

    def _get_pqr(self):
        """Run PDB2PQR to get charge for each atom

        TODO: figure out why unknown residues fro hetatms are absent
        """
        if self.pqr is None:
            pdbfd, tmp_pdb_path = tempfile.mkstemp()
            with os.fdopen(pdbfd, 'w') as tmp:
                # writer.set_structure(self.structure)
                # writer.save(tmp)
                # tmp.seek(0)
                self.save_pdb(tmp, True)

            _, tmp_pqr_path = tempfile.mkstemp()

            with silence_stdout(), silence_stderr():
                mainCommand(["pdb2pqr.py", "--ff=amber", "--whitespace", tmp_pdb_path, tmp_pqr_path])

            os.remove(tmp_pdb_path)

            self.pqr = {}
            with open(tmp_pqr_path) as pqr:
                for line in pqr:
                    for line in pqr:
                        if line.startswith("REMARK"): continue

                        fields = line.rstrip().split()
                        if len(fields) == 11:
                            recordName, serial, atomName, residueName, chainID, residueNumber, X, Y, Z, charge, radius = fields
                        elif len(fields) == 10:
                            recordName, serial, atomName, residueName, residueNumber, X, Y, Z, charge, radius = fields
                        else:
                            print len(fields)
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
            os.remove(tmp_pqr_path)

        return self.pqr

    def _get_cx(self):
        if self.cx is None:
            pdbfd, tmp_pdb_path = tempfile.mkstemp()
            with os.fdopen(pdbfd, 'w') as tmp:
                self.save_pdb(tmp, True)

            with open(tmp_pdb_path) as f:
                cx_f = subprocess.check_output("cx", stdin=f)

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
                        print "    Error, maybe the next line contains it?"
            # self.cx = {int(l[6:11].strip()):float(l[60:66].strip()) \
            #     for l in cx_f.splitlines() if l[:6].strip() in ["ATOM", "HETATM"]}

        return self.cx

    def _mean_coord(self):
        if not self.mean_coord_updated:
            self.mean_coord = np.mean(self.get_coords(), axis=0)
            self.mean_coord_updated = True
        return self.mean_coord

    def get_coords(self, include_hetatms=False, exclude_atoms=None):
        return np.array([a.get_coord() for a in self.get_atoms(
            include_hetatms=include_hetatms, exclude_atoms=exclude_atoms)])

    def orient_to_pai(self, flip_axis=(0.2, 0.2, 0.2)):
        coords = PCA(n_components = 3).fit_transform(self.get_coords())
        coords = flip_around_axis(coords, axis=flip_axis)
        self.update_coords(coords)

    def rotate(self, num=1):
        """Rotate structure in randomly in place"""
        for r in xrange(num):
            M = rotation_matrix(random=True)
            coords = np.dot(self.get_coords(), M)
            self.update_coords(coords)
            yield r

    def update_coords(self, coords):
        for atom, coord in it.izip(self.structure.get_atoms(), coords):
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
                raise InvalidPDB()

    def align_seq_to_struc(self, *seq_num, **kwds):
        return_residue=kwds.get("return_residue", False)

        if return_residue:
            mapped_residues = [self.get_residue_from_resseq(pdbnum) \
                for current_resi, resn, pdbnum, ncbi in map_residues(self.pdb, self.chain, seq_num)]
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

    def get_features(self, input_shape=(96, 96, 96), residue_list=None, batch_index=None, return_data=True, return_truth=True, only_aa=False, only_atom=False, expand_atom=True, include_full_protein=False, voxel_size=0.5):
        # return self.get_features_in_grid(
        #     input_shape=input_shape,
        #     residue_list=residue_list,
        #     batch_index=batch_index,
        #     return_data=return_data,
        #     return_truth=return_truth,
        #     only_aa=only_aa,
        #     expand_atom=expand_atom,
        #     include_full_protein=include_full_protein,
        #     voxel_size=voxel_size)
        if self.course_grained:
            return self.map_residues_to_voxel_space(
                binding_site_residues=residue_list,
                include_full_protein=include_full_protein,
                only_aa=only_aa
            )
        return self.map_atoms_to_voxel_space(
            expand_atom=expand_atom,
            binding_site_residues=residue_list,
            include_full_protein=include_full_protein,
            only_aa=only_aa,
            only_atom=only_atom)

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

            num_voxels = sum([atom_spheres[elem].shape[0]*count for elem, count in elems.iteritems()])
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

    def expand_atom(self, atom):
        """
        """
        try:
            atom_sphere_origin = atom_spheres[atom.element.title()]
        except KeyError:
            raise RuntimeError("invalid atom")

        coord = atom.get_coord()
        coord += [self.volume/2]*3
        center = np.digitize(coord, self.voxels)

        atom_sphere = np.rint(atom_sphere_origin+center).astype("int64")

        return atom_sphere

    def get_neighbors(self, atom, radius=2.):
        if self.kdtree is None:
            self.kdtree = NeighborSearch(self.get_atoms())

        neighbors = self.kdtree.search(atom.get_coord(), radius)
        return neighbors

    def map_atoms_to_voxel_space(self, expand_atom=True, binding_site_residues=None, include_full_protein=False, only_aa=False, only_atom=False):
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
        exclude_atoms = None #["N", "C", "O"] #Ignore backbone

        if binding_site_residues is not None:
            if not include_full_protein:
                atoms = sorted((a for r in binding_site_residues for a in r), key=lambda a: a.get_serial_number())
                atoms = list(self.filter_atoms(atoms, exclude_atoms=exclude_atoms))
                binding_site_atoms = [a.get_serial_number() for a in atoms]
            else:
                atoms = list(self.get_atoms(include_hetatms=True, exclude_atoms=exclude_atoms))
                nAtoms = len(atoms)
                binding_site_atoms = [a.get_serial_number() for r in binding_site_residues for a in r]
        else:
            atoms = list(self.get_atoms(include_hetatms=True, exclude_atoms=exclude_atoms))
            nAtoms = len(atoms)
            binding_site_atoms = []

        if only_atom:
            nFeatures = 5
        elif only_aa:
            nFeatures = 21
        else:
            nFeatures = Structure.nFeatures

        if expand_atom:
            atom_points = self.expand_atom
        else:
            atom_points = lambda a: [self.get_grid_coord(a)]

        data_voxels = defaultdict(lambda: np.zeros(nFeatures))
        truth_voxels = defaultdict(lambda: np.zeros(1))

        for atom in atoms:
            atom_sphere = atom_points(atom)
            features = self.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom)
            truth = int(atom.get_serial_number() in binding_site_atoms)

            for atom_grid in atom_sphere:
                atom_grid = tuple(atom_grid.tolist())
                data_voxels[atom_grid] += features
                truth_voxels[atom_grid] = np.array([truth])

        if binding_site_residues is None:
            return np.array(data_voxels.keys()), np.array(data_voxels.values())
        else:
            truth = np.array([truth_voxels[grid] for grid in data_voxels.keys()])
            return np.array(data_voxels.keys()), np.array(data_voxels.values()), truth

    def map_residues_to_voxel_space(self, binding_site_residues=None, include_full_protein=False, only_aa=False):
        if binding_site_residues is not None:
            if not include_full_protein:
                residues = binding_site_residues
                binding_site_residues = [r.get_id() for r in residues]
            else:
                residues = self.structure.get_residues()
                binding_site_residues = [r.get_id()[1] for r in binding_site_residues]
        else:
            residues = self.structure.get_residues()
            binding_site_residues = []

        if only_aa:
            nFeatures = 21
        else:
            nFeatures = 36

        data_voxels = {}
        truth_voxels = {}
        voxel_coords = {}

        for residue in residues:
            grid, coord = self.get_grid_coord_for_residue(residue)
            grid = tuple(grid.tolist())
            if grid in data_voxels:
                print "Ignoring {}.{} {} {} becuase it falls in grid {} with {}".format(
                    self.pdb,
                    self.chain,
                    residue.get_resname(),
                    residue.get_id()[1],
                    grid,
                    voxel_coords[grid]
                )

                assert 0
                #voxel_coords[grid]
            else:
                data_voxels[grid] = self.get_features_for_residue(residue, only_aa=only_aa)
                truth_voxels[grid] = [int(residue.get_id()[1] in binding_site_residues)]
                voxel_coords[grid] = (residue.get_resname(), residue.get_id()[1], coord)


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

    def get_features_for_atom(self, atom, only_aa=False, only_atom=False, preload=False):
        """Calculate FEATUREs"""
        if not self.feature_file and preload:
            print "preload features"
            features = self.precalc_features[atom.serial_number-1]
            if only_atom:
                return features[13:18]
            elif only_aa:
                return features[35:56]
            else:
                return features
        elif only_atom:
            return self.get_element_type(atom)
        elif only_aa:
            return self.get_residue(atom)
        else:
            features = np.zeros(Structure.nFeatures)

            #atom_type
            features[0:13]  = self.get_atom_type(atom)
            #partial_charge        = None
            features[13:18] = self.get_element_type(atom)
            #hydroxyl              = None
            #amide                 = None
            #amine                 = None
            #carbonyl              = None
            #ring_system           = None
            #peptide               = None
            features[18:19] = self.get_vdw(atom)
            features[19:23] = self.get_charge(atom)
            #charge_with_his       = None
            features[23:27] = self.get_concavity(atom)
            features[27:31] = self.get_hydrophobicity(atom)
            #mobility              = None
            features[31:35] = self.get_accessible_surface_area(atom)

            features[35:56] = self.get_residue(atom)
            # residue_class1        = self.get_residue_class2(atom, one_hot=True)
            # residue_class2        = self.get_residue_class2(atom, one_hot=True)
            features[56:59]  = self.get_ss(atom)

            if self.feature_file:
                self.precalc_features[atom.serial_number-1] = features

            return features

    def get_features_for_residue(self, residue, only_aa=False, preload=False):
        """Calculate FEATUREs"""
        if not self.feature_file and preload:
            print "preload features"
            features = self.precalc_features[atom.serial_number-1]
            if only_aa:
                return features[15:36]
            else:
                return features
        elif only_aa:
            return self.get_residue(residue)
        else:
            features = np.zeros(Structure.nFeatures)

            features[0:4] = self.get_charge(residue)
            features[4:8] = self.get_concavity(residue)
            features[8:12] = self.get_hydrophobicity(residue)
            features[12:15] = self.get_accessible_surface_area(residue)
            features[15:36] = self.get_residue(residue)
            features[36:39]  = self.get_ss(residue)

            if self.feature_file:
                self.precalc_features[atom.serial_number-1] = features

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
        return np.array([vdw.get(atom.get_name().strip()[0], 0.0)])

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
            raise RuntimeErorr("Input must be Atom or Residue")

        pqr = self._get_pqr()
        atom_id = atom.get_full_id()[3:5]

        if atom_id[1][1] != " ":
            #pdb2pqr removes alternate conformations and only uses the first
            atom_id = (atom_id[0], (atom_id[1][0], " "))
        charge_value = pqr.get(atom_id, np.NaN)

        charge = np.zeros(4)
        charge[0] = charge_value
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
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if is_atom:
            try:
                sasa, sasa_struct = self._get_sasa()
                selection = "sele, chain {} and resi {} and name {}".format(self.chain, atom.get_parent().get_id()[1], atom.get_id()[0])
                with silence_stdout(), silence_stderr():
                    selections = freesasa.selectArea([selection], sasa_struct, sasa)
                    atom_area = selections["sele"]
            except (KeyError, AssertionError, AttributeError, TypeError):
                atom_area = np.NaN
        else:
            dssp = self._get_dssp()
            try:
                #("1abc", 0, "A", (" ", 10, "A"), (self.name, self.altloc))
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
            asa = np.zeros(4)
            asa[0] = atom_area
            asa[1] = residue_area
            asa[2] = float(residue_area < 0.2)
            asa[3] = float(residue_area >= 0.2)
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

    def get_grid_coord_for_atom(self, atom, vsize=96, max_radius=40, homothetic_transformation=False, round=False, voxel_size=None):
        """Convert atom cartesian coordiantes to voxel space using a homothetic transformation
        """
        if voxel_size is not None and voxel_size != self.voxel_size:
            self.set_voxel_size(voxel_size)

        #Center at origin
        #print "original:", atom.coord
        new_coord = np.array(atom.coord) - self._mean_coord()
        #print "at origin:", new_coord

        if homothetic_transformation:
            adjusted = (int(vsize/2.-1)/float(max_radius))*new_coord
            new_coord = adjusted + (vsize-1)/2 # Translate center
            voxel = new_coord.astype(int) # Round components
        else:
            #No transformation, just fit inside a (256,256,256) volume
            new_coord += [self.volume/2]*3
            #print "at 126:", new_coord
            if round:
                new_coord = np.around(new_coord)
            voxel = np.digitize(new_coord, self.voxels)


        # if half_angstrom:
        #     is_upper_bound = (np.divmod(new_coord, 1)[1]>=0.5).astype(int)
        #     new_coord = 2*np.floor(new_coord)+is_upper_bound
        #
        # if False:
        #     max_gap = 1./voxel_size
        #     upper_bound = np.divmod(new_coord, 1)[1]
        #     new_coord = max_gap*p.floor(new_coord)+


        #print "voxel:", rounded
        return voxel

    def get_grid_coord_for_residue(self, residue, round=False):
        coord = residue["CA"].get_coord() #np.mean([a.get_coord() for a in residue], axis=0)
        coord += [self.volume/2]*3
        if round:
            coord = np.around(coord)
        #print residue.get_resname(), coord,
        voxel = np.digitize(coord, self.voxels)
        #print voxel
        return voxel, coord

    def set_volume(self, volume):
        self.volume = volume

    def set_voxel_size(self, voxel_size=None):
        if not self.course_grained:
            assert voxel_size is not None
            self.voxel_size = voxel_size
        else:
            self.voxel_size = 3.43
        self.voxels = np.arange(0, self.volume, self.voxel_size)
        make_atom_spheres(self.voxel_size)


        #np.multiply((vsize/2.)-1)/float(max_radius), coords)

    def get_atoms_from_grids(self, grids, vsize=96, max_radius=40):
        for atom in self.get_atoms():
            coord = self.get_grid_coord(atom, vsize, max_radius)
            print coord, grids
            if coord.tolist() in grids:
                yield atom

def download_pdb(id):
    pdbl = PDB.PDBList()
    try:
        fname = pdbl.retrieve_pdb_file(id.upper(), file_format="mmCif")
        if not os.path.isfile(fname):
            raise InvalidPDB(id)
        return fname, "mmcif"
    except IOError:
        raise InvalidPDB(id)

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

    return M

if __name__ == "__main__":
  import sys
  assert len(sys.argv) > 1
  structure = Structure(sys.argv[1], panfs=False).extract_chain("A")
  num_atoms = sum(1 for a in structure.structure.get_atoms() if a.get_parent().get_id()[0] == " ")
  for atom in structure.get_atoms():
    features = structure.get_features_for_atom(atom)
    print atom.get_full_id(), atom.serial_number, "of", num_atoms, features, len(features)
