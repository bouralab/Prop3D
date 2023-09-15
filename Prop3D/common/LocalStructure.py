import os
import sys
import copy
from io import StringIO
from functools import partial
from collections import defaultdict
from collections.abc import Iterator
from typing import Union, Any, TypeVar, IO, AnyStr

import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.decomposition import PCA

from Bio import PDB
from Bio.PDB import Selection
from Bio.PDB.NeighborSearch import NeighborSearch

from toil.realtimeLogger import RealtimeLogger

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from Prop3D.util import natural_keys
from Prop3D.util.pdb import InvalidPDB
from Prop3D.common.ProteinTables import vdw_radii, vdw_aa_radii
from Prop3D.common.features import all_features

AtomType = TypeVar('AtomType', bound='PDB.Atom')
ResidueType = TypeVar('ResidueType', bound='PDB.Residue')
_Self = TypeVar('_Self', bound='LocalStructure')

class LocalStructure(object):
    """An object to read in a local structure files (PDB, mmCIF) and manipulate proteins in cartiesian space.
    Backend: Bio.PDB. TODO: Convert to use abstract structure

    Parameters
    ----------
    path : str
        Path to structure file on local filesystem
    cath_domain : str
        Name of domain used to save the files later on
    feature_mode : str (r, w, r+, w+):
        Open file for (r)eading or (w)riting. Default reading.
    features_path : str or None
        DEPRECATED. Path to .h5 file for features. Not used anymore.
    residue_feature_mode : str
        DEPRECATED. Open features for reading or writing
    reset_chain : bool
        Get chain name from structure file, not cath_domain. Default is False.
    volume : float
        DEPRACATED. volume for voxels. Use Distributed VoxelizedStructure.
    """
    def __init__(self, path: str, cath_domain: str, input_format: str = "pdb",
                 feature_mode: str = "r", features_path: Union[str, None] = None, residue_feature_mode: str = "r",
                 reset_chain: bool = False, volume: float = 256.) -> None:
        self.path = path
        if not os.path.isfile(self.path):
            raise InvalidPDB("Cannot find file {}".format(self.path))

        if self.path.endswith(".gz"):
            raise InvalidPDB("Error with {}. Gzipped archives not allowed. Please use constructor or util.get_pdb. File: {}".format(pdb, self.path))

        path = Path(path)
        if input_format in ["pdb", "pqr"] or (input_format == "guess" and path.suffix in [".pdb", ".pqr"]):
            parser = PDB.PDBParser()
            self.input_format = "pdb"
        elif input_format == "mmcif"  or (input_format == "guess" and path.suffix == ".cif"):
            parser = PDB.FastMMCIFParser()
            self.input_format = "mmcif"
        elif input_format == "mmtf" or (input_format == "guess" and path.suffix == ".mmtf"):
            parser = PDB.MMTFParser()
            self.input_format = "mmtf"
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
        self.n_residue_features = len(all_features.residue_features)
        self.n_atom_features = len(all_features.atom_features)

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
            self.atom_features = all_features.default_atom_feature_df(len(atom_index)).assign(serial_number=atom_index)
            self.atom_features = self.atom_features.set_index("serial_number")

        if self.residue_feature_mode == "r" and os.path.isfile(self.residue_features_file):
            self.residue_features = pd.read_hdf(self.residue_features_file, "table", mode="r")
        else:
            het, resi, ins = zip(*[self._remove_inscodes(r).get_id() for r in self.structure.get_residues()])
            self.residue_features = all_features.default_residue_feature_df(len(het)).assign(HET_FLAG=het, resi=resi, ins=ins)
            self.residue_features = self.residue_features.set_index(["HET_FLAG", "resi", "ins"])

        self.atom_feature_names = copy.copy(all_features.atom_features)
        self.residue_feature_names = copy.copy(all_features.residue_features)

        self.other_formats = defaultdict(lambda: partial(self.save_pdb, path=f"{self.cath_domain}.pdb") if self.input_format != "pdb" else self.path)

    def __abs__(self) -> _Self:
        """Take the absolue value of all atom features
        """
        new = self.copy()
        new.atom_features = new.atom_features.abs()
        return new

    def __sub__(self, other: _Self) -> _Self:
        """Subtract atom feature values from other structure from this stucture.
        """
        new = self.copy()
        if isinstance(other, LocalStructure):
            new.atom_features -= other.atom_features
        elif isinstance(other, (int, float)):
            new.atom_features -= other
        else:
            raise TypeError
        return new

    def __add__(self, other: _Self) -> _Self:
        """Add atom feature values from other structure from this stucture.
        """
        new = self.copy()
        if isinstance(other, LocalStructure):
            new.atom_features += other.atom_features
        elif isinstance(other, (int, float)):
            new.atom_features += other
        else:
            raise TypeError
        return new

    def __floordiv__(self, other: _Self) -> _Self:
        """Divide atom feature values from other structure from this stucture.
        """
        new = self.copy()
        if isinstance(other, LocalStructure):
            new.atom_features /= other.atom_features
        elif isinstance(other, (int, float)):
            new.atom_features /= other
        else:
            raise TypeError
        return new

    def __truediv__(self, other: _Self) -> _Self:
        """Divide atom feature values from other structure from this stucture.
        """
        return self.__floordiv__(other)

    def normalize_features(self, columns: Union[str, list[str], None] = None) -> _Self:
        """Normalize features using min max scaling

        Parameters
        ----------
        columns: str or list of strs
            Names of feature columns to normalize
        """
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

    def copy(self, empty: bool = False) -> _Self:
        """Create a deep copy of current structure.

        Parameters
        ----------
        (deprecated) empty: bool
            Don't copy features
        """
        new = copy.deepcopy(self)
        return new

        #Make sure original files do not get overwritten
        new.features_path = os.getcwd()
        unique_id = int.from_bytes(os.urandom(4), sys.byteorder)
        new.atom_features_file = os.path.join(new.features_path , f"{os.path.basename(self.atom_features_file)}_{unique_id}.h5")
        new.residue_feature_file = os.path.join(new.features_path , f"{os.path.basename(self.atom_features_file)}_{unique_id}.h5")
        return new

    def __deepcopy__(self, memo: Any) -> _Self:
        """
        """
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

    def get_atoms(self, include_hetatms: bool = False, exclude_atoms: Union[list[AtomType], None] = None, 
                  include_atoms: Union[list[AtomType], None] = None) -> Iterator[AtomType]:
        """Enumerate protein model for all atoms with options to filter

        Parameters
        ----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        for a in self.filter_atoms(self.structure.get_atoms(), include_hetatms=include_hetatms, exclude_atoms=exclude_atoms, include_atoms=include_atoms):
            yield a

    def get_surface(self, level: str = "R") -> Iterator[AtomType]:
        """Returns all surface atoms, using DSSP accessible surface value"
        """
        if self.atom_features["residue_buried"].astype(int).sum() == 0:
            raise RuntimeError("Must calculate features with featurizer before running this")

        surface = self.atom_features[self.atom_features["residue_buried"]==False]

        surface_atoms = [a for a in self.get_atoms() if a in surface.index]

        if level == "A":
            return surface_atoms

        return Selection.unfold_entities(surface_atoms, level)

    def filter_atoms(self, atoms, include_hetatms=False, exclude_atoms: Union[list[AtomType], None] = None, 
                     include_atoms: Union[list[AtomType], None] = None) -> Iterator[AtomType]:
        """Enumerate protein model for all atoms with options to filter

        Parameters
        ----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        for a in atoms:
            hetflag, resseq, icode = a.get_parent().get_id()
            if not include_hetatms and hetflag != ' ':
                continue
            if exclude_atoms is not None and a.get_name().strip() in exclude_atoms:
                continue
            if include_atoms is not None and a.serial_number not in include_atoms:
                continue
            yield a

    def get_residues(self) -> Iterator[ResidueType]:
        for r in self.structure.get_residues():
            yield self._remove_inscodes(r)

    def save_pdb(self, path: Union[str, None] = None, header: Union[str, None] = None, 
                 file_like: Union[str, None] = False, rewind: bool = True) -> Union[str, IO[AnyStr]]:
        """Write PDB to file

        Parameters
        ----------
        path : None or str
            Path to save PDB file. If None, file_like needs to be True.
        header : str or list of strs
            Header string to write to the beginning of each PDB file
        file_like : boolean
            Return a StringIO object of the PDB file, do not write to disk. Default False.
        rewind : boolean
            If returning a file-like object, rewind the beginning of the file

        Returns
        -------
        None or file-like object of PDB file data
        """
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

    def write_features(self, features: Union[list[str], None] = None, coarse_grained: bool = False, 
                       name: Union[str, None] = None, work_dir: Union[str, None] = None) -> None:
        """Write features to a spefic file, depnignd on protein loading class, e.g. HDF

        Parameters
        ----------
        features : str or list of strs
            Features to write
        course_grained : boolean
            Include features only at the residue level. Default False
        name : str
            File name to write file
        work_dir : None or str
            Directory to save file
        """
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

    def write_features_to_pdb(self, features_to_use: Union[list[str], None] = None, name: Union[str, None] = None, 
                              coarse_grain=False, work_dir: bool = None, other: Any = None) -> list[dict[str, str]]:
        """Write features to PDB files in the bfactor column. One feature per PDB file.

        Parameters
        ----------
        features_to_use: str or list of strs
            Features to write
        name : None or str
            File name to write file
        course_grained : boolean
            Include features only at the residue level. Default False
        work_dir : None or str
            Directory to save files.
        other : obj
            Ignored. May be useful for subclasses.

        Returns
        -------
        outfiles : dict 
            File names of each written PDB file for each feature
        """
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

    def add_features(self, coarse_grained: bool = False, **features: Union[pd.Series, np.array, list[float]]) -> None:
        """Add a feature column to dataset

        Parameters
        ----------
        coarse_grained : bool
            Use residue level. Else use atom level. Defualt False.
        **features : 
            use param_name = list of features for each atom or residue 
        """
        if coarse_grained:
            assert [len(f)==len(self.residue_features) for f in features.values()]
            self.residue_features = self.residue_features.assign(**features)
            self.residue_feature_names += list(features.values())
        else:
            assert [len(f)==len(self.atom_features) for f in features.values()]
            self.atom_features = self.atom_features.assign(**features)
            self.atom_feature_names += list(features.keys())

    def get_pdb_dataframe(self, coarse_grained: bool = False, include_featuresd: bool = False) -> pd.DataFrame:
        """Get standard data from PDB as a data frame

        Parameters
        ----------
        coarse_grained : bool
            Use residue level. Else use atom level. Defualt False.
        include_features : list
            Include calculated features into DataFrame. Default is False
        """
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

    def get_residue_from_resseq(self, resseq: Union[str, int], model: int = 0, 
                                chain: Union[str, None] = None) -> Union[ResidueType, None]:
        """Quickly access residue from Bio.PDB by resseq, handling erros

        Parameters
        ----------
        resseq : str or int
            Residue number to search for '22' or '22A'
        model : int
            PDB model number. Default 0.
        chain : None or str
            PDB chain to search. If None, use the chain specific in init. Dfault None.
        """
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

    def align_seq_to_struc(self, *seq_num: Any, **kwds: Any) -> list[ResidueType]:
        """DEPRACATED. Convert PDB, UniProt, and MMDB numberings
        """
        raise NotImplementedError
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

    def get_mean_coord(self) -> np.array:
        """Get the mean XYZ coordinate or center of mass.
        """
        if not self.mean_coord_updated:
            self.mean_coord = np.around(np.mean(self.get_coords(), axis=0), decimals=4)
            self.mean_coord_updated = True
        return self.mean_coord

    def get_max_coord(self) -> np.array:
        """Get the maximum coordinate in each dimanesion
        """
        return np.max(self.get_coords(), axis=0)

    def get_min_coord(self) -> np.array:
        """Get the minimum coordinate in each dimanesion
        """
        return np.min(self.get_coords(), axis=0)

    def get_max_length(self, buffer: float = 0., pct_buffer: float = 0.) -> float:
        """Get the length of the protein to create a volume around

        Parameters
        ----------
        buffer : float
            Amount of space to incldue around volume in Angstroms. Defualt 0
        pct_buffer : float
            Amount of space to incldue around volume in as percentage of the total legnth in Angstroms. Defualt 0

        Returns
        -------
        length : float
            Max length of protein
        """
        length = int(np.ceil(np.linalg.norm(self.get_max_coord()-self.get_min_coord())))
        if pct_buffer!=0:
            length += int(np.ceil(length*pct_buffer))
        else:
            length += buffer
        if length%2 == 1:
            length += 1
        return length

    def shift_coords(self, new_center: Union[np.array, None] = None, 
                     from_origin: Union[np.array, None] = True) -> np.array:
        """Shift coordinates by setting a new center of mass value or shift to the origin
        
        Parameters
        ----------
        new_center : 3-tuple of floats or None
            XYZ cooridnate of new center. if new_center is None, it will shift to the origin. Default is None.
        from_origin : bool
            Start shift from the origin by first subratcting center of mass. Defualt is True.

        Returns
        -------
        The new center coordinate
        """
        coords = self.get_coords()
        if from_origin or new_center is None:
            mean_coord = self.get_mean_coord()
            coords -= mean_coord
        if new_center is not None:
            coords += np.around(new_center, decimals=4)

        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.mean(coords, axis=0)

    def shift_coords_to_origin(self) -> np.array:
        """Center structure at the origin

        Returns
        -------
        The new center coordinate
        """
        return self.shift_coords()

    def shift_coords_to_volume_center(self) -> np.array:
        """Move structure to the center of the volume
        """
        return self.shift_coords(np.array([self.volume/2]*3))

    def resize_volume(self, new_volume, shift=True):
        self.volume = new_volume
        if shift:
            self.shift_coords_to_volume_center()

    def get_coords(self, include_hetatms: bool = False, exclude_atoms: Union[list[ResidueType], None] = None) -> np.array:
        return np.array([a.get_coord() for a in self.get_atoms(
            include_hetatms=include_hetatms, exclude_atoms=exclude_atoms)]).round(decimals=4)

    def orient_to_pai(self, random_flip: bool = False, flip_axis: Union[list[float], np.array] = (0.2, 0.2, 0.2)) -> None:
        """Orient structure to the Principle Axis of Interertia and optionally flip. Modified from EnzyNet

        Parameters
        ----------
        random_flip : bool
            Randomly flip around axis. Defualt is False.
        flip_axis: 3-tuple of floats
            New axis to flip.
        """
        self.shift_coords_to_origin()

        coords = PCA(n_components = 3).fit_transform(self.get_coords())
        if random_flip:
            coords = flip_around_axis(coords, axis=flip_axis)

        self.update_coords(coords)

        self.shift_coords_to_volume_center()

    def rotate(self, rvs: Union[np.array, None] = None, num: int = 1) -> Iterator[tuple[int, np.array]]:
        """Rotate structure by either randomly in place or with a set rotation matrix. 
        Random rotations matrices are drawn from the Haar distribution (the only uniform 
        distribution on SO(3)) from scipy.

        Parameters
        ----------
        rvs : np.array (3x3)
            A rotation matrix. If None, a randome roation matrix is used. Default is None.
        num : int
            Number of rotations to perfom
        return_to : XYZ coordinate
            When finsihed rotating, move structure to this coordinate. Defualt is to the center of mass

        Yields
        ------
        r : int
            Rotation number
        M : np.array (3x3)
            Rotation matrix
        """
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

    def update_coords(self, coords: Union[np.array, list[float]]) -> None:
        """Update XYZ coordinates with a new set of coordinates for the same atoms"""
        for atom, coord in zip(self.structure.get_atoms(), coords):
            atom.set_coord(coord)
        self.mean_coord = None
        self.mean_coord_updated = False

    def update_bfactors(self, b_factors: Union[np.array, list[float]]) -> None:
        """Update bfactors with a new set of bfactors for the same atoms"""
        for atom, b in zip(self.structure.get_atoms(), b_factors):
            atom.set_bfactor(b)

    def _remove_altloc(self, atom: AtomType) -> AtomType:
        """Get the first atom if there are multiple altlocs
        """
        if isinstance(atom, PDB.Atom.Atom):
            return atom
        elif isinstance(atom, PDB.Atom.DisorderedAtom):
            return atom.disordered_get_list()[0]
        else:
            raise RuntimeError("Invalid atom type")

    def _remove_inscodes(self, residue: ResidueType) -> ResidueType:
        """Get the first residue if there are multiple residues for same location
        """
        if isinstance(residue, PDB.Residue.Residue):
            return residue
        elif isinstance(residue, PDB.Residue.DisorderedResidue):
            return residue.disordered_get_list()[0]
        else:
            raise RuntimeError("Invalid residue type")

    def calculate_neighbors(self, d_cutoff: float = 100.0, level: str="R") -> Union[list[AtomType], list[ResidueType]]:
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

    def get_vdw(self, atom_or_residue: Union[AtomType, ResidueType]) -> np.array:
        """Get van der walls radii for an atom or residue
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            return np.array([vdw_radii.get(atom_or_residue.element.title(), 1.7)])
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            #For coarse graining, not really a vdw radius
            return np.array([vdw_aa_radii.get(atom_or_residue.get_resname(), 3.0)])

    def get_dihedral_angles(self, atom_or_residue: Union[AtomType, ResidueType]) -> Union[tuple[float, float], tuple[None, None]]:
        """Get deidral angle for atom (mapped up to residue) or the residue
        """
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


def flip_around_axis(coords : np.array, axis: Union[tuple[float, float, float], np.array] = (0.2, 0.2, 0.2)) -> np.array:
    'Flips coordinates randomly w.r.t. each axis with its associated probability'
    for col in range(3):
        if np.random.binomial(1, axis[col]):
            coords[:,col] = np.negative(coords[:,col])
    return coords

def rotation_matrix(random: bool = False, theta: Union[int, float] = 0, phi: Union[int, float] = 0, 
                    z: Union[int, float] = 0, uniform: bool =True) -> tuple[np.array, float, float, float]:
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

def unit_vector(vector: np.array) -> np.array:
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1: np.array, v2: np.array) -> np.array:
    """ Get the angle in radians between vectors 'v1' and 'v2'

    Examples:
        angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
        angle_between((1, 0, 0), (1, 0, 0))
            0.0
        angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def get_dihedral(p0: np.array, p1: np.array, p: np.array, p3: np.array) -> float:
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
