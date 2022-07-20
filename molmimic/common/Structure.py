import os
import sys
import copy
from io import StringIO

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
from Bio import PDB
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection
from toil.realtimeLogger import RealtimeLogger

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from molmimic.util import natural_keys
from molmimic.util.pdb import InvalidPDB
from molmimic.common.ProteinTables import vdw_radii, vdw_aa_radii
from molmimic.common.features import default_atom_feature_df, default_residue_feature_df, \
    atom_features, residue_features

class Structure(object):
    def __init__(self, path, cath_domain, input_format="pdb",
                 feature_mode="r", features_path=None, residue_feature_mode="r",
                 reset_chain=False, volume=256):
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

        if cath_domain is not None:
            self.cath_domain = self.sdi = cath_domain
            self.pdb = cath_domain[:4]
            self.chain = cath_domain[4]
            self.domNo = cath_domain[5:]
        else:
            self.cath_domain = self.pdb = os.path.splitext(os.path.basename(path))[0]
            self.domNo = "00"
            reset_chain = True

        self.volume = volume

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

        if reset_chain:
            self.chain = all_chains[0].id

        self.id = self.cath_domain #"{}{}{:02d}".format(self.pdb, self.chain, int(self.domNo))
        self.n_residue_features = len(residue_features)
        self.n_atom_features = len(atom_features)

        if features_path is None:
            features_path = os.environ.get("MOLMIMIC_FEATURES", os.getcwd())

        if features_path.endswith("_atom.h5"):
            self.features_path = os.path.dirname(features_path)
            self.atom_features_file = features_path
            self.residue_features_file = os.path.join(self.features_path,
                "{}_residue.h5".format(self.id))
        elif features_path.endswith("_residue.h5"):
            self.residue_features_file = features_path
            self.atom_features_file = os.path.join(self.features_path,
                "{}_atom.h5".format(self.id))
        else:
            self.atom_features_file = os.path.join(features_path,
                "{}_atom.h5".format(self.id))
            self.residue_features_file = os.path.join(features_path,
                "{}_residue.h5".format(self.id))

        self.features_path = os.path.dirname(self.atom_features_file)

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

        if self.atom_feature_mode == "r":
            self.atom_features = pd.read_hdf(self.atom_features_file, "table", mode="r")
        else:
            atom_index = [self._remove_altloc(a).serial_number for a in self.structure.get_atoms()]
            self.atom_features = default_atom_feature_df(len(atom_index)).assign(serial_number=atom_index)
            self.atom_features = self.atom_features.set_index("serial_number")

        if self.residue_feature_mode == "r" and os.path.isfile(self.residue_features_file):
            self.residue_features = pd.read_hdf(self.residue_features_file, "table", mode="r")
        else:
            het, resi, ins = zip(*[self._remove_inscodes(r).get_id() for r in self.structure.get_residues()])
            self.residue_features = default_residue_feature_df(len(het)).assign(HET_FLAG=het, resi=resi, ins=ins)
            self.residue_features = self.residue_features.set_index(["HET_FLAG", "resi", "ins"])

        self.atom_feature_names = copy.copy(atom_features)
        self.residue_feature_names = copy.copy(residue_features)

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

    def normalize_features(self, columns=None):
        new = self.copy()

        if columns is not None:
            if not isinstance(columns, (list, tuple)):
                columns = [columns]
            data = new.atom_features[columns]
        else:
            data = new.atom_features

        min_max_scaler = preprocessing.MinMaxScaler()
        data_scaled = min_max_scaler.fit_transform(data.values)

        if columns is not None:
            new.atom_features.loc[:, columns] = data_scaled
        else:
            new.atom_features.loc[:] = data_scaled

        return new

    def copy(self, empty=False):
        new = copy.deepcopy(self)
        return new

        #Make sure original files do not get overwritten
        new.features_path = os.getcwd()
        unique_id = int.from_bytes(os.urandom(4), sys.byteorder)
        new.atom_features_file = os.path.join(new.features_path , f"{os.path.basename(self.atom_features_file)}_{unique_id}.h5")
        new.residue_feature_file = os.path.join(new.features_path , f"{os.path.basename(self.atom_features_file)}_{unique_id}.h5")
        return new

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        empty = False
        for k, v in self.__dict__.items():
            if k == "atom_features":
                setattr(result, "atom_features", pd.DataFrame(
                    np.nan if empty else np.array(copy.deepcopy(self.atom_features.values.tolist())),
                    index=pd.Index(copy.deepcopy(self.atom_features.index.tolist())),
                    columns=copy.deepcopy(self.atom_features.columns.tolist())))
            elif k == "residue_features":
                setattr(result, "residue_features", pd.DataFrame(
                    np.nan if empty else np.array(copy.deepcopy(self.residue_features.values.tolist())),
                    index=pd.Index(copy.deepcopy(self.residue_features.index.tolist())),
                    columns=copy.deepcopy(self.residue_features.columns.tolist())))


        setattr(result, "features_path", os.getcwd())
        unique_id = int.from_bytes(os.urandom(4), sys.byteorder)
        setattr(result, "atom_features_file", os.path.join(result.features_path , f"{os.path.splitext(os.path.basename(result.atom_features_file))[0]}_{unique_id}.h5"))
        setattr(result, "residue_feature_file", os.path.join(result.features_path , f"{os.path.splitext(os.path.basename(result.atom_features_file))[0]}_{unique_id}.h5"))

        return result

    def get_atoms(self, include_hetatms=False, exclude_atoms=None, include_atoms=None):
        for a in self.filter_atoms(self.structure.get_atoms(), include_hetatms=include_hetatms, exclude_atoms=exclude_atoms, include_atoms=include_atoms):
            yield a

    def get_surface(self, level="R"):
        if self.atom_features["residue_buried"].astype(int).sum() == 0:
            raise RuntimeError("Must calculate features with featurizer before running this")

        surface = self.atom_features[self.atom_features["residue_buried"]==False]

        surface_atoms = [a for a in self.get_atoms() if a in surface.index]

        if level == "A":
            return surface_atoms

        return Selection.unfold_entities(surface_atoms, level)

    def filter_atoms(self, atoms, include_hetatms=False, exclude_atoms=None, include_atoms=None):
        for a in atoms:
            hetflag, resseq, icode = a.get_parent().get_id()
            if not include_hetatms and hetflag != ' ':
                continue
            if exclude_atoms is not None and a.get_name().strip() in exclude_atoms:
                continue
            if include_atoms is not None and a.serial_number not in include_atoms:
                continue
            yield a

    def save_pdb(self, path=None, header=None, file_like=False, rewind=True):
        lines = not path and not file_like
        if path is None:
            path = StringIO()

        if header is not None:
            new_header = ""
            for line in header.splitlines():
                if not line.startswith("REMARK"):
                    line = "REMARK {}\n".format(line)
                new_header += line
            old_header = self.structure.header
            self.structure.header = new_header

        writer = PDB.PDBIO()
        writer.set_structure(self.structure)
        writer.save(path)

        if file_like and rewind:
            path.seek(0)

        if lines:
            path = path.read()

        if header is not None:
            self.structure.header = old_header

        return path

    def write_features(self, features=None, coarse_grained=False, name=None, work_dir=None):
        if work_dir is not None or name is not None:
            if work_dir is None:
                work_dir = self.features_path

            if name is not None:
                residue_features_file = os.path.join(work_dir, name)
                atom_features_file = os.path.join(work_dir, name)
            else:
                residue_features_file = os.path.join(work_dir, os.path.basename(self.residue_features_file))
                atom_features_file = os.path.join(work_dir, os.path.basename(self.atom_features_file))
        else:
            residue_features_file = self.residue_features_file
            atom_features_file = self.atom_features_file

        if features is not None:
            if not isinstance(features, (list, tuple)):
                features = [features]
        else:
            features = self.residue_feature_names if coarse_grained else self.atom_feature_names

        print("Feats to save", features, atom_features_file)

        if coarse_grained:
            self.residue_features = self.residue_features.astype(np.float64)
            self.residue_features = self.residue_features.drop(columns=
                [col for col in self.residue_features if col not in \
                features]) #self.residue_feature_names])
            self.residue_features.to_hdf(residue_features_file, "table")
        else:
            self.atom_features = self.atom_features.astype(np.float64)
            self.atom_features = self.atom_features.drop(columns=
                [col for col in self.atom_features if col not in \
                features]) #self.atom_feature_names])
            self.atom_features.to_hdf(atom_features_file, "table")

    def write_features_to_pdb(self, features_to_use=None, name=None, coarse_grain=False, work_dir=None, other=None):
        if work_dir is None:
            work_dir = os.getcwd()

        if features_to_use is None:
            features = self.atom_features if not coarse_grain else \
                self.residue_features
        else:
            features = self.atom_features.loc[:, features_to_use] if not coarse_grain \
                else self.residue_features.loc[:, features_to_use]

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

    def add_features(self, coarse_grained=False, **features):
        if coarse_grained:
            assert [len(f)==len(self.residue_features) for f in features.values()]
            self.residue_features = self.residue_features.assign(**features)
            self.residue_feature_names += list(features.values())
        else:
            assert [len(f)==len(self.atom_features) for f in features.values()]
            self.atom_features = self.atom_features.assign(**features)
            self.atom_feature_names += list(features.keys())

    def get_pdb_dataframe(self, coarse_grained=False, include_features=False):
        str_type = np.dtype('O', metadata={'vlen': str})
        na = {float:np.nan, str_type:"", int:9999999}

        if coarse_grained:
            df = pd.DataFrame([
                [
                    *residue.id,
                    "".join(map(str, residue.id[1:])).strip(),
                    residue.get_resname(),
                    residue.get_parent().id,
                    np.mean([a.get_bfactor() for a in residue]),  # isotropic B factor
                    *np.mean([a.get_coord() for a in residue], axis=0)
                ] for residue in self.structure.get_residues()],
                columns=["HET_FLAG", "resi", "ins", "residue_id", "residue_name",
                         "chain", "bfactor", "X", "Y", "Z"])
            df = df.set_index(["HET_FLAG", "resi", "ins"])

            na = {float:np.nan, str_type:""}
            for col, dtype in [
              ("residue_id", str_type),
              ("residue_name", str_type),
              ("chain", str_type),
              ("bfactor", str_type),
              ("X", float),
              ("Y", float),
              ("Z", float)]:
                df[col] = df[col].fillna(na[dtype]).astype(dtype)

            if include_features:
                df = pd.merge(df, self.residue_features, left_index=True, right_index=True)
                df = df.reset_index(drop=True)
                #pd.concat((df, self.residue_features), axis=1)
        else:
            df = pd.DataFrame([
                [
                    atom.serial_number,
                    atom.get_fullname(),
                    "".join(map(str, atom.get_parent().id[1:])).strip(),
                    atom.get_parent().get_resname(),
                    atom.get_parent().get_parent().id,
                    atom.get_bfactor(),  # isotropic B factor
                    *atom.coord
                ] for atom in self.structure.get_atoms()],
                columns=["serial_number", "atom_name", "residue_id", "residue_name",
                         "chain", "bfactor", "X", "Y", "Z"])


            for col, dtype in [
              ("serial_number", str_type),
              ("atom_name", str_type),
              ("residue_id", str_type),
              ("residue_name", str_type),
              ("chain", str_type),
              ("bfactor", float),
              ("X", float),
              ("Y", float),
              ("Z", float)]:
                df[col] = df[col].fillna(na[dtype]).astype(dtype)

            if include_features:
                df = pd.merge(df, self.atom_features.reset_index(), on="serial_number")
                #df = pd.concat((df, self.atom_features), axis=1)
        return df

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
            self.mean_coord = np.around(np.mean(self.get_coords(), axis=0), decimals=4)
            self.mean_coord_updated = True
        return self.mean_coord

    def get_max_coord(self):
        return np.max(self.get_coords(), axis=0)

    def get_min_coord(self):
        return np.min(self.get_coords(), axis=0)

    def get_max_length(self, buffer=0, pct_buffer=0):
        length = int(np.ceil(np.linalg.norm(self.get_max_coord()-self.get_min_coord())))
        if pct_buffer!=0:
            length += int(np.ceil(length*pct_buffer))
        else:
            length += buffer
        if length%2 == 1:
            length += 1
        return length

    def shift_coords(self, new_center=None, from_origin=True):
        """if new_center is None, it will shift to the origin"""
        coords = self.get_coords()
        if from_origin or new_center is None:
            mean_coord = self.get_mean_coord()
            coords -= mean_coord
        if new_center is not None:
            coords += np.around(new_center, decimals=4)

        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.mean(coords, axis=0)

    def shift_coords_to_origin(self):
        return self.shift_coords()

    def shift_coords_to_volume_center(self):
        return self.shift_coords(np.array([self.volume/2]*3))

    def resize_volume(self, new_volume, shift=True):
        self.volume = new_volume
        if shift:
            self.shift_coords_to_volume_center()

    def get_coords(self, include_hetatms=False, exclude_atoms=None):
        return np.array([a.get_coord() for a in self.get_atoms(
            include_hetatms=include_hetatms, exclude_atoms=exclude_atoms)]).round(decimals=4)

    def orient_to_pai(self, random_flip=False, flip_axis=(0.2, 0.2, 0.2)):
        self.shift_coords_to_origin()

        coords = PCA(n_components = 3).fit_transform(self.get_coords())
        if random_flip:
            coords = flip_around_axis(coords, axis=flip_axis)

        self.update_coords(coords)

        self.shift_coords_to_volume_center()

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
            coords = np.dot(self.get_coords(), M).round(decimals=4)
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
