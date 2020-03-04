import os
import copy
from io import StringIO

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from Bio import PDB
from Bio.PDB.NeighborSearch import NeighborSearch

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from molmimic.util import natural_keys
from molmimic.util.pdb import InvalidPDB
from molmimic.common.ProteinTables import vdw_radii, vdw_aa_radii
from molmimic.common.features import default_atom_feature_df, default_residue_feature_df, \
    atom_features, residue_features

class Structure(object):
    def __init__(self, path, cath_domain, input_format="pdb",
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

        self.cath_domain = self.sdi = cath_domain
        self.pdb = cath_domain[:4]
        self.chain = cath_domain[4]
        self.domNo = cath_domain[5:]

        try:
            self.structure = parser.get_structure(self.cath_domain, self.path)
        except KeyError:
            #Invalid mmcif file
            raise InvalidPDB("Invalid PDB file: {} (path={})".format(self.cath_domain, self.path))

        try:
            all_chains = list(self.structure[0].get_chains())
        except (KeyError, StopIteration):
            raise InvalidPDB("Error get chains for {} {}".format(self.cath_domain, self.path))

        if len(all_chains) > 1:
            raise InvalidPDB("Only accepts PDBs with 1 chain in {} {}".format(self.cath_domain, self.path))

        self.id = "{}{}{:02d}".format(self.pdb, self.chain, int(self.domNo))
        self.n_residue_features = len(residue_features)
        self.n_atom_features = len(atom_features)

        self.features_path = os.environ.get("MOLMIMIC_FEATURES", features_path)

        if features_path.endswith("_atom.h5"):
            self.features_path = os.path.dirname(features_path)
            self.atom_features_file = features_file
            self.residue_features_file = os.path.join(self.features_path,
                "{}_residue.h5".format(self.id))
        elif features_path.endswith("_residue.h5"):
            self.residue_features_file = features_file
            self.atom_features_file = os.path.join(self.features_path,
                "{}_atom.h5".format(self.id))
        else:
            self.atom_features_file = os.path.join(self.features_path,
                "{}_atom.h5".format(self.id))
            self.residue_features_file = os.path.join(self.features_path,
                "{}_residue.h5".format(self.id))

        if not os.path.isdir(os.path.dirname(os.path.abspath(self.atom_features_file))):
            os.makedirs(os.path.abspath(os.path.dirname(self.atom_features_file)))

        if False and feature_mode == "r" and not os.path.isfile(self.atom_features_file):
            self.atom_feature_mode = "w+"
        else:
            self.atom_feature_mode = feature_mode

        if False and feature_mode == "r" and not os.path.isfile(self.residue_features_file):
            self.residue_feature_mode = "w+"
        else:
            self.residue_feature_mode = residue_feature_mode

        atom_index = [self._remove_altloc(a).serial_number for a in self.structure.get_atoms()]
        if self.atom_feature_mode == "r":
            self.atom_features = pd.read_hdf(self.atom_features_file, "table", mode="r")
        else:
            self.atom_features = default_atom_feature_df(len(atom_index)).reindex(atom_index, axis=0)

        if self.residue_feature_mode == "r":
            self.residue_features = pd.read_hdf(self.residue_features_file, "table", mode="r")
        else:
            residue_index = pd.MultiIndex.from_tuples(
                [self._remove_inscodes(r).get_id() for r in \
                self.structure.get_residues()], names=('HET_FLAG', 'resi', 'ins'))
            self.residue_features = default_residue_feature_df(len(residue_index)).reindex(residue_index, axis=0)

    def __abs__(self):
        new = self.copy()
        new.atom_features = new.atom_features.abs()
        return new

    def __sub__(self, other):
        new = self.copy()
        if isinstance(other, Structure):
            new.atom_features -= other.atom_features
        elif isinstance(other, (int, float)):
            new.atom_features -= other
        else:
            raise TypeError
        return new

    def __add__(self, other):
        new = self.copy()
        if isinstance(other, Structure):
            new.atom_features += other.atom_features
        elif isinstance(other, (int, float)):
            new.atom_features += other
        else:
            raise TypeError
        return new

    def __floordiv__(self, other):
        new = self.copy()
        if isinstance(other, Structure):
            new.atom_features /= other.atom_features
        elif isinstance(other, (int, float)):
            new.atom_features /= other
        else:
            raise TypeError
        return new

    def __truediv__(self, other):
        return self.__floordiv__(other)

    def copy(self):
        return copy.copy(self)

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

    def write_features_to_pdb(self, features=None, name=None, coarse_grain=False,
      work_dir=None):
        if work_dir is None:
            work_dir = os.getcwd()

        if features is None:
            features = self.atom_features if not coarse_grain else \
                self.residue_features
        else:
            features = self.atom_features.loc[:, features] if not coarse_grain \
                else self.residue_features.loc[:, features]

        bfactors = [a.bfactor for a in self.structure.get_atoms()]

        path = os.path.join(work_dir, self.cath_domain)
        if name is not None:
            path += "-"+name

        outfiles = {}
        for feature in features.columns:
            self.update_bfactors(features.loc[:, feature]*100)
            outfile = "{}-{}.pdb".format(path, feature)
            self.save_pdb(outfile)
            outfiles[feature] = outfile

        self.update_bfactors(bfactors)

        return outfiles

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
            # if rvs is None or rvs!=np.eye(3):
            #     assert not np.array_equal(old_coords, self.get_coords()), M
            self.shift_coords_to_volume_center()
            # if rvs is None or rvs!=np.eye(3):
            #     assert not np.array_equal(coords, self.get_coords()), M
            yield r, M

    def update_coords(self, coords):
        for atom, coord in zip(self.structure.get_atoms(), coords):
            atom.set_coord(coord)
        self.mean_coord = None
        self.mean_coord_updated = False

    def update_bfactors(self, b_factors):
        for atom, b in zip(self.structure.get_atoms(), b_factors):
            atom.set_bfactor(b)

    def _remove_altloc(self, atom):
        if isinstance(atom, PDB.Atom.Atom):
            return atom
        elif isinstance(atom, PDB.Atom.DisorderedAtom):
            return atom.disordered_get_list()[0]
        else:
            raise RuntimeError("Invalid atom type")

    def _remove_inscodes(self, residue):
        if isinstance(residue, PDB.Residue.Residue):
            return residue
        elif isinstance(residue, PDB.Residue.DisorderedResidue):
            return residue.disordered_get_list()[0]
        else:
            raise RuntimeError("Invalid residue type")

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

    def get_vdw(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            return np.array([vdw_radii.get(atom_or_residue.element.title(), 1.7)])
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            #For coarse graining, not really a vdw radius
            return np.array([vdw_aa_radii.get(atom_or_residue.get_resname(), 3.0)])

    def get_dihedral_angles(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
        else:
            raise RuntimeError("Input must be Atom or Residue")

        if not hasattr(self, "_phi_psi"):
            peptides = PDB.PPBuilder().build_peptides(self.structure[0][self.chain])[0]
            self._phi_psi = {res.get_id():ppl for res, ppl in zip(
                peptides, peptides.get_phi_psi_list())}

        return self._phi_psi.get(residue.get_id(), (None, None))


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
