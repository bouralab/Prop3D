import os
from io import StringIO

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from Bio import PDB
from Bio.PDB.NeighborSearch import NeighborSearch

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from molmimic.util import InvalidPDB, natural_keys
from molmimic.common.ProteinTables import vdw_radii

def number_of_features(only_aa=False, only_atom=False, non_geom_features=False,
  use_deepsite_features=False, course_grained=False):
    if course_grained:
        if non_geom_features:
            return 35
        elif only_aa:
            return 21
        else:
            return 45
    else:
        if non_geom_features:
            return 13
        elif only_atom:
            return 5
        elif only_aa and use_deepsite_features:
            return 21+9
        elif only_aa:
            return 21
        elif use_deepsite_features:
            return 9
        else:
            return 73

def get_feature_names(only_aa=False, only_atom=False, use_deepsite_features=False, course_grained=False):
    if only_atom:
        return ["C", "N", "O", "S", "Unk_element"]
    elif only_aa and use_deepsite_features:
        return PDB.Polypeptide.aa3+[ #DeepSite features
            "hydrophobic_atom",
            "aromatic_atom",
            "hbond_acceptor",
            "hbond_donor",
            "metal",
            "is_hydrogen",
            "is_pos",
            "is_neg",
            "is_conserved"
        ]
    if only_aa:
        return PDB.Polypeptide.aa3
    if use_deepsite_features:
        return [ #DeepSite features
            "hydrophobic_atom",
            "aromatic_atom",
            "hbond_acceptor",
            "hbond_donor",
            "metal",
            "is_hydrogen",
            "is_pos",
            "is_neg",
            "is_conserved"
        ]
    feature_names = [
        "C", "CT", "CA", "N", "N2", "N3", "NA", "O", "O2", "OH", "S", "SH", "Unk_atom",
        "C_elem", "N_elem", "O_elem", "S_elem", "Unk_element",
        "vdw_volume", "charge", "neg_charge", "pos_charge", "neutral_charge",
        "electrostatic_potential", "is_electropositive", "is_electronegative",
        "cx", "is_concave", "is_convex", "is_concave_and_convex",
        "hydrophobicity", "is_hydrophbic", "is_hydrophilic", "hydrophobicity_is_0",
        "atom_asa", "atom_is_buried", "atom_exposed",
        "residue_asa", "residue_buried", "residue_exposed"]
    feature_names += PDB.Polypeptide.aa3+["Unk_residue"]
    feature_names += ["is_helix", "is_sheet", "Unk_SS"]
    feature_names += [ #DeepSite features
        "hydrophobic_atom",
        "aromatic_atom",
        "hbond_acceptor",
        "hbond_donor",
        "metal",
        "is_hydrogen",
    ]
    feature_names += [
        "conservation_normalized",
        "conservation_scale",
        "is_conserved"
    ]

    return feature_names

def get_residue_feature_names():
    feature_names = [
        "charge", "neg_charge", "pos_charge", "neutral_charge",
        "electrostatic_potential", "is_electropositive", "is_electronegative",
        "cx", "is_concave", "is_convex", "is_concave_and_convex",
        "hydrophobicity", "is_hydrophbic", "is_hydrophilic", "hydrophobicity_is_0",
        "residue_asa", "residue_buried", "residue_exposed"]
    feature_names += PDB.Polypeptide.aa3
    feature_names += [
        "Unk_residue", "is_helix", "is_sheet", "Unk_SS",
        "conservation_normalized", "conservation_scale", "is_conserved"
    ]
    return feature_names

class Structure(object):
    def __init__(self, path, pdb, chain, sdi, domNo, input_format="pdb",
      feature_mode="r", features_path=None, residue_feature_mode="r"):
        self.path = path
        if not os.path.isfile(self.path):
            raise InvalidPDB("Cannot find file {}".format(self.path))

        if self.path.endswith(".gz"):
            raise InvalidPDB("Error with {}. Gzipped archives not allowed. Please use constructor or util.get_pdb. File: {}".format(pdb, self.path))

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

        self.pdb = pdb
        self.chain = chain
        self.sdi = sdi
        self.domNo = domNo

        try:
            all_chains = list(self.structure[0].get_chains())
        except (KeyError, StopIteration):
            raise InvalidPDB("Error get chains for {} {}".format(pdb, self.path))

        if len(all_chains) > 1:
            raise InvalidPDB("Only accepts PDBs with 1 chain in {} {}".format(pdb, self.path))

        self.id = "{}_{}_sdi{}_d{}".format(self.pdb, self.chain, self.sdi, self.domNo)

        self.n_residues = max(r.get_id()[1] for r in self.structure.get_residues())
        self.n_atoms = max(a.serial_number for a in self.structure.get_atoms())
        self.n_residue_features = number_of_features(course_grained=True)
        self.n_atom_features = number_of_features(course_grained=False)

        self.features_path = os.environ.get("MOLMIMIC_FEATURES", features_path)
        
        if features_path.endswith("_atom.h5"):
            self.features_path = os.path.dirname(features_path)
            self.atom_features_file = features_file
            self.residue_features_file = os.path.join(self.features_path,
            self.pdb.lower()[1:3], "{}_residue.npy".format(self.id))
        elif features_path.endswith("_residue.npy"):
            self.features_path = os.path.dirname(features_path)
            self.residue_features_file = features_file 
            self.atom_features_file = os.path.join(self.features_path,
            self.pdb.lower()[1:3], "{}_atom.npy".format(self.id))
        else:
            self.atom_features_file = os.path.join(self.features_path,
            self.pdb.lower()[1:3], "{}_atom.npy".format(self.id))
            self.residue_features_file = os.path.join(self.features_path,
            self.pdb.lower()[1:3], "{}_residue.npy".format(self.id))

        if not os.path.isdir(os.path.dirname(self.atom_features_file)):
            os.makedirs(os.path.dirname(self.atom_features_file))

        if False and feature_mode == "r" and not os.path.isfile(self.atom_features_file):
            self.atom_feature_mode = "w+"
        else:
            self.atom_feature_mode = feature_mode

        if False and feature_mode == "r" and not os.path.isfile(self.residue_features_file):
            self.residue_feature_mode = "w+"
        else:
            self.residue_feature_mode = residue_feature_mode

        # self.residue_features = np.memmap(self.residue_features_file,
        #     dtype=np.float, mode=self.residue_feature_mode, shape=(self.n_residues, self.n_residue_features))
        # self.atom_features = np.memmap(self.atom_features_file,
        #     dtype=np.float, mode=self.atom_feature_mode, shape=(self.n_atoms, self.n_atom_features))

        atom_index = [a.serial_number for a in self.structure.get_atoms()]
        if self.atom_feature_mode == "r":
            self.atom_features = pd.read_hdf(self.atom_features_file, "table", mode="r")
            self.atom_features = self.atom_features.reset_index()
            self.atom_features = self.atom_features.reindex(index=atom_index)
        else:
            self.atom_features = pd.DataFrame(columns=get_feature_names(),
                index=atom_index)

        if self.residue_feature_mode == "r":
            self.residue_features = pd.read_hdf(self.residue_features_file, "table", mode="r")
        else:
            residue_index = [r.get_id() for r in self.structure.get_residues()]
            self.residue_features = pd.DataFrame(columns=get_residue_feature_names(),
                index=residue_index)

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

    def get_residue_from_resseq(self, resseq, model=0, chain=None):
        chain = chain or self.chain
        try:
            return self.structure[model][chain][resseq]
        except KeyError as e:
            for r in self.structure[model][chain]:
                if r.get_id() == resseq:
                    return r
                if isinstance(resseq, int) and r.get_id()[1] == resseq:
                    return r
                if isinstance(resseq, str):
                    resseq_parts = natural_keys(resseq)
                    res_seq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
                    try:
                        return self.structure[model][chain][res_seq]
                    except KeyError:
                        return None
            else:
                return None

    def align_seq_to_struc(self, *seq_num, **kwds):
        return_residue = kwds.get("return_residue", False)
        use_mmdb_index = kwds.get("return_residue", False)

        residues = map_residues(self.pdb, self.chain, seq_num, use_mmdb_index=use_mmdb_index)

        if return_residue:
            mapped_residues = [self.get_residue_from_resseq(pdbnum) \
                for current_resi, resn, pdbnum, ncbi in residues]

            if mapped_residues.count(None) > len(mapped_residues)/2.:
                pct_missing = 100*mapped_residues.count(None)/float(len(mapped_residues))
                raise InvalidPDB("Binding Site ({:4}%) missing from structure ({}.{})".format(pct_missing, self.pdb, self.chain))

            mapped_residues = [r for r in mapped_residues if r is not None]
        else:
            mapped_residues = [pdb for current_resi, resn, pdb, ncbi in residues]

        return mapped_residues

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

    def rotate(self, rvs=None, num=1):
        """Rotate structure in randomly in place"""
        for r in range(num):
            if rvs is None:
                M, theta, phi, z = rotation_matrix(random=True)
                #print(M, theta, phi, z)
            else:
                M=rvs
            self.shift_coords_to_origin()
            old_coords = self.get_coords()
            coords = np.dot(self.get_coords(), M)
            self.update_coords(coords)
            assert not np.array_equal(old_coords, self.get_coords()), M
            self.shift_coords_to_volume_center()
            assert not np.array_equal(coords, self.get_coords()), M
            yield r, M

    def update_coords(self, coords):
        for atom, coord in zip(self.structure.get_atoms(), coords):
            atom.set_coord(coord)
        self.mean_coord = None
        self.mean_coord_updated = False

    def calculate_neighbors(self, d_cutoff=100.0, level="R"):
        """
        Calculates intermolecular contacts in a parsed struct object.
        Modified from haddocking/prodigy

        Parameters
        ----------
        struct : Bio.PDB.Structure
            The structure object
        d_cuttoff: float
            Distance to find neighbors

        Returns
        -------
        A list of lists of nearby elements at the specified level: [(a1,b2),]
        """
        atom_list = list(self.structure.get_atoms())
        ns = NeighborSearch(atom_list)
        all_list = ns.search_all(radius=d_cutoff, level=level)

        if not all_list:
            raise ValueError('No contacts found for selection')

        return all_list

    def get_vdw(self, atom):
        return np.array([vdw_radii.get(atom.element.title(), 2.0)])

def flip_around_axis(coords, axis = (0.2, 0.2, 0.2)):
    'Flips coordinates randomly w.r.t. each axis with its associated probability'
    for col in range(3):
        if np.random.binomial(1, axis[col]):
            coords[:,col] = np.negative(coords[:,col])
    return coords

def rotation_matrix(random = False, theta = 0, phi = 0, z = 0, uniform=True):
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

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def get_dihedral(p0, p1, p2, p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

    return M, theta, phi, z
