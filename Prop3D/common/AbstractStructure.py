import os
import copy
import numpy as np
import pandas as pd
from scipy.stats import special_ortho_group

from Prop3D.common.ProteinTables import three_to_one

class AbstractStructure(object):
    """A base structure class that holds default methods for dealing with structures.
    Subclasses can be created to use different protein libraires, e.g. BioPython, BioTite,
    our own HDF/HSDS distributed stucture.

    Parameters:
    -----------
    name: str
        Name of proten, e.g. PDB or CATH id
    file_mode: str (r, w, r+, w+, etc)
        Open file for reading or writing. Defualt is just reading, no methods will affect underlying file
    coarse_grained: boolean
        Use a residue only model instead of an all atom model. Defualt False. Warning, not fully implemented.
    """
    def __init__(self, name, file_mode="r", coarse_grained=False):
        assert hasattr(self, "features"), "Must sublass and instialize features as a recarray or pd.DataFrame"
        self.name = name
        self.file_mode = file_mode
        self.coarse_grained = coarse_grained

        self.coords = self.get_coords()

        if isinstance(self.features, pd.DataFrame):
            self.is_dataframe = True
        else:
            self.is_dataframe = False

        self.mean_coord_updated = False

    def copy(self, empty=False):
        """Create a deep copy of current structure.

        Parameters:
        -----------
        (deprecated) empty: bool
            Don't copy features
        """
        new = copy.deepcopy(self)
        return new

    def deep_copy_feature(self, feature_name):
        """Subclass this method to handle custom copying of specific features

        Parameters:
        -----------
        feature_name: str
            Feature name to copy

        Raises
        ------
        NotImeplementedError if no method to handle feature
        """
        raise NotImplementedError

    def __deepcopy__(self, memo):
        """
        """
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        empty = False
        for k, v in self.__dict__.items():
            try:
                new_value = self.deep_copy_feature(k, memo)
            except NotImplementedError:
                new_value = copy.deepcopy(v, memo)

            setattr(result, k, new_value)
            
        return result

    def __abs__(self):
        """Take the absolue value of all features
        """
        new = self.copy()
        new.features = new.features.abs()
        return new

    def __sub__(self, other):
        """Subtract feature values from other structure from this stucture.
        """
        new = self.copy()
        if isinstance(other, AbstractStructure):
            assert self.coarse_grained == other.coarse_grained, "Must both be at atom or residue level"
            new.features -= other.features
        elif isinstance(other, (int, float)):
            new.features -= other
        else:
            raise TypeError
        return new

    def __add__(self, other):
        """Add feature values from other structure from this stucture.
        """
        new = self.copy()
        if isinstance(other, AbstractStructure):
            assert self.coarse_grained == other.coarse_grained, "Must both be at atom or residue level"
            new.features += other.features
        elif isinstance(other, (int, float)):
            new.features += other
        else:
            raise TypeError
        return new

    def __floordiv__(self, other):
        """Divide feature values from other structure from this stucture.
        """
        new = self.copy()
        if isinstance(other, AbstractStructure):
            assert self.coarse_grained == other.coarse_grained, "Must both be at atom or residue level"
            new.features /= other.features
        elif isinstance(other, (int, float)):
            new.features /= other
        else:
            raise TypeError
        return new

    def __truediv__(self, other):
        """Divide feature values from other structure from this stucture.
        """
        return self.__floordiv__(other)

    def normalize_features(self, columns=None):
        """Normalize features using min max scaling

        Parameters:
        -----------
        columns: str or list of strs
            Names of feature columns to normalize
        """
        new = self.copy()

        if columns is not None:
            if not isinstance(columns, (list, tuple)):
                columns = [columns]
            data = new.features[columns]
        else:
            data = new.features

        min_max_scaler = preprocessing.MinMaxScaler()
        data_scaled = min_max_scaler.fit_transform(data.values \
            if self.is_dataframe else data)

        if self.is_dataframe:
            if columns is not None:
                new.features.loc[:, columns] = data_scaled
            else:
                new.features.loc[:] = data_scaled
        else:
            if columns is not None:
                new.features[columns] = data_scaled.astype(data.dtype)
            else:
                new.features.loc[:] = data_scaled.astype(data.dtype)

        return new

    def get_atoms(self, include_hetatms=False, exclude_atoms=None, include_atoms=None):
        """Subclass to enumerate protein model for all atoms with options to filter

        Parameters:
        -----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        raise NotImplementedError

    def filter_atoms(atoms=None, include_hetatms=False, exclude_atoms=None, include_atoms=None):
        """Subclass to enumerate protein model for all atoms with options to filter

        Parameters:
        -----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        yield from self.get_atoms(atoms, include_hetatms=include_hetatms,
            exclude_atoms=exclude_atoms, include_atoms=include_atoms)

    def get_surface(self):
        """Returns all surface atoms, using DSSP accessible surface value"
        """
        surface = self.features[self.features["residue_buried"]==0]
        return surface

    def get_bfactors(self):
        """Get bfactors for all atoms
        """
        raise NotImplementedError

    def save_pdb(self, path=None, header=None, file_like=False, rewind=True):
        """Write PDB to file

        Parameters:
        -----------
        path : None or str
            Path to save PDB file. If None, file_like needs to be True.
        header : str or list of strs
            Header string to write to the beginning of each PDB file
        file_like : boolean
            Return a StringIO object of the PDB file, do not write to disk. Default False.
        rewind : boolean
            If returning a file-like object, rewind the beginning of the file

        Returns:
        --------
        None or file-like object of PDB file data
        """
        raise NotImplementedError

    def write_features(self, features=None, coarse_grained=False, name=None, work_dir=None):
        """Subclass to write features to a spefic file, depnignd on protein loading class, e.g. HDF

        Parameters:
        -----------
        features : str or list of strs
            Features to write
        course_grained : boolean
            Include features only at the residue level. Default False
        name : str
            File name to write file
        work_dir : None or str
            Directory to save file
        """
        raise NotImplementedError

    def write_features_to_pdb(self, features_to_use=None, name=None, coarse_grain=False, work_dir=None, other=None):
        """Write features to PDB files in the bfactor column. One feature per PDB file.

        Parameters:
        -----------
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
        File names of each written PDB file for each feature
        """
        if work_dir is None:
            work_dir = os.getcwd()

        if features_to_use is None:
            features = self.features
            if self.is_dataframe:
                features_to_use = self.features.columns
            else:
                features_to_use = self.features.dtype.names
        else:
            features = self.features[features_to_use]

        bfactors = self.get_bfactors()

        path = os.path.join(work_dir, self.name.replace("/", "_"))
        if name is not None:
            path += "-"+name

        outfiles = {}
        for feature in features_to_use:
            self.update_bfactors(features[feature]*100)
            outfile = "{}-{}.pdb".format(path, feature)
            self.save_pdb(outfile)
            outfiles[feature] = outfile

        self.update_bfactors(bfactors)

        return outfiles

    def add_features(self, coarse_grained=False, **features):
        """Subclass to add a feature column to dataset
        """
        raise NotImplementedError

    def get_coords(self, include_hetatms=False, exclude_atoms=None):
        """Subclass to return all XYZ coordinates fro each atom

        Parameters:
        -----------
        include_hetatms : bool
            Include heteroatoms or not. Default is False.
        exclude_atoms : bool
            Atoms to exlude while getting coordinates

        Returns:
        --------
        XYZ coordinates from each atom in a specified format of subclass
        """
        raise NotImplementedError
    
    def get_sequence(self):
        """Get amino acid sequence from structure
        """
        sequence = ""
        prev_resi = None
        for a in self.get_atoms():
            resi = a["residue_id"]
            if resi != prev_resi:
                resn = three_to_one(a["residue_name"].decode("utf-8"))
                sequence += resn
        
        return sequence

    def update_coords(self, coords):
        """Sublcass to create method to update XYZ coordinates with a new set of coordinates for the same atoms"""
        self.coords = coords
        self.mean_coord = None
        self.mean_coord_updated = False

    def get_mean_coord(self):
        """Get the mean XYZ coordinate or center of mass.
        """
        if not self.mean_coord_updated:
            self.mean_coord = np.around(np.nanmean(self.coords, axis=0), decimals=4)
            self.mean_coord_updated = True
        return self.mean_coord

    def get_max_coord(self):
        """Get the maximum coordinate in each dimanesion
        """
        return np.nanmax(self.coords, axis=0)

    def get_min_coord(self):
        """Get the minimum coordinate in each dimanesion
        """
        return np.nanmin(self.coords, axis=0)

    def get_max_length(self, buffer=0, pct_buffer=0):
        """Get the length of the protein to create a volume around

        Parameters:
        -----------
        buffer : float
            Amount of space to incldue around volume in Angstroms. Defualt 0
        pct_buffer : float
            Amount of space to incldue around volume in as percentage of the total legnth in Angstroms. Defualt 0

        Returns:
        --------
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

    def shift_coords(self, new_center=None, from_origin=True):
        """Shift coordinates by setting a new center of mass value or shift to the origin
        
        Parameters:
        -----------
        new_center : 3-tuple of floats or None
            XYZ cooridnate of new center. if new_center is None, it will shift to the origin. Default is None.
        from_origin : bool
            Start shift from the origin by first subratcting center of mass. Defualt is True.

        Returns:
        --------
        The new center coordinate
        """
        coords = self.coords
        if from_origin or new_center is None:
            mean_coord = self.get_mean_coord()
            coords -= mean_coord
        if new_center is not None:
            coords += np.around(new_center, decimals=4)

        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.nanmean(coords, axis=0)

    def shift_coords_to_origin(self):
        """Center structure at the origin

        Returns:
        --------
        The new center coordinate
        """
        return self.shift_coords()

    def orient_to_pai(self, random_flip=False, flip_axis=(0.2, 0.2, 0.2)):
        """Orient structure to the Principle Axis of Interertia and optionally flip. Modified from EnzyNet

        Parameters:
        -----------
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

    def rotate(self, rvs=None, num=1, return_to=None):
        """Rotate structure by either randomly in place or with a set rotation matrix. 
        Random rotations matrices are drawn from the Haar distribution (the only uniform 
        distribution on SO(3)) from scipy.

        Parameters:
        -----------
        rvs : np.array (3x3)
            A rotation matrix. If None, a randome roation matrix is used. Default is None.
        num : int
            Number of rotations to perfom
        return_to : XYZ coordinate
            When finsihed rotating, move structure to this coordinate. Defualt is to the center of mass

        Yields:
        -------
        r : int
            Rotation number
        M : np.array (3x3)
            Rotation matrix
        """
        for r in range(num):
            if rvs is None:
                M, theta, phi, z = special_ortho_group.rvs(3)
            else:
                M=rvs
            self.shift_coords_to_origin()
            old_coords = self.get_coords()
            coords = np.dot(self.coords, M).round(decimals=4)
            self.update_coords(coords)

            self.shift_coords(self.get_mean_coord() if return_to is None else return_to)

            yield r, M

    def update_bfactors(self, b_factors):
        """Sublcass to create method to update bfactors with a new set of bfactors for the same atoms"""
        raise NotImplementedError

    def calculate_neighbors(self, d_cutoff=100.0):
        """Subclass to find all nearest neighbors within a given radius.

        Parameters:
        -----------
        d_cutoff : float
            Distance cutoff to find neighbors. Deualt is 100 Angtroms
        """
        raise NotImplementedError

    def get_vdw(self, atom_or_residue):
        """Get van der walls radii for an atom or residue
        """
        raise NotImplementedError

    def get_dihedral_angles(self, atom_or_residue):
        """Get deidral angle for atom (mapped up to residue) or the residue
        """
        raise NotImplementedError

    def get_secondary_structures_groups(self, verbose=False):
        """Get groups of adjecent atom rows the belong to the same secondary structure.
        We use DSSP to assing secondary structures to each reisdue mapped down to atoms. 
        If any segment was <4 residues, they were merged the previous and next groups

        Returns:
        --------
        ss_groups : list of pd.DataFrames of atoms in each group
        loop_for_ss : dict of pd.DataFrames of atoms in the loop following each ss group
        original_order : dict
        ss_type : dict
        leading_trailing_residues : dict
        number_ss : int
        """
        assert not self.coarse_grained

        if not self.is_dataframe:
            ss_type = pd.DataFrame.from_records(
                self.data[["serial_number", "is_helix", "is_sheet", "Unk_SS"]],
                index="serial_number")
            get_id = lambda x: x["residue_id"]
        else:
            ss_type = self.features[["is_helix", "is_sheet", "Unk_SS"]]
            get_id = lambda x: x.get_id()

        ss_type = ss_type.rename(columns={"is_helix":"H", "is_sheet":"E", "Unk_SS":"X"})
        ss_type = ss_type.idxmax(axis=1)

        ss_groups = ss_type.groupby([(ss_type != ss_type.shift()).cumsum()-1])

        #Merge group shorter than 4 residues
        for i, ss_group in ss_groups:
            if 0<i<ss_groups.ngroups-1:
                this_group = ss_group.iloc[0]
                prev_group = ss_groups.get_group(i-1).iloc[0]
                next_group = ss_groups.get_group(i+1).iloc[0]
                this_group_atoms = tuple(self.get_atoms(include_atoms=ss_group.index))
                this_group_residues = tuple(self.unfold_entities(this_group_atoms, "R"))

                if verbose:
                    this_group_id = tuple(np.unique(r["residue_id"])[0] for r in this_group_residues)
                    print(this_group_id[0], this_group_id[-1], this_group)

                if len(this_group_residues)<4 and this_group != "X":
                    ss_type.loc[ss_group.index] = "X"
                if len(this_group_residues)<3 and prev_group == next_group:
                    ss_type.loc[ss_group.index] = prev_group
                elif len(this_group_residues)<3 and next_group == "X" and this_group != prev_group:
                    ss_type.loc[ss_group.index] = "X"

                if this_group=="H" and prev_group=="E" and next_group=="E" and len(this_group_residues)<5:
                    ss_type.loc[ss_group.index] = "X"

                if len(this_group_residues)<3: #OB=5
                    if prev_group == next_group:
                        ss_type.loc[ss_group.index] = prev_group
                    else:
                        pass

                    if this_group=="H" and prev_group=="E" and next_group=="E":
                        ss_type.loc[ss_group.index] = "X"
                    elif this_group=="E" and prev_group=="H" and next_group=="H":
                        ss_type.loc[ss_group.index] = "X"

                if len(this_group_residues)>10 and this_group=="X":
                    ss_type.loc[ss_group.index] = "H"

        #Regroup with correct SS
        ss_atom_groups = ss_type.groupby([(ss_type != ss_type.shift()).cumsum()-1])

        ss_groups = []
        loop_for_ss = {}
        original_order = {}
        ss_type = {}
        leading_trailing_residues = {}

        for i, ss_group in ss_atom_groups:
            #Get all atoms from SS and loops
            ss_atoms = tuple(self.get_atoms(include_atoms=ss_group.index))
            ss_residues = tuple(self.unfold_entities(ss_atoms, "R"))
            ss_residues_id = tuple(np.unique(r["residue_id"])[0] for r in ss_residues)

            if ss_group.iloc[0] != "X":
                ss_groups.append(ss_residues)
                original_order[ss_residues_id] = len(ss_groups)
                ss_type[ss_residues_id] = ss_group.iloc[0]
            elif len(ss_groups)>0 and ss_group.iloc[0] == "X":
                loop_for_ss[ss_residues_id] = ss_residues

        first_group = ss_atom_groups.get_group(0)
        if first_group.iloc[0] == "X":
            loop_atoms = tuple(self.get_atoms(include_atoms=first_group.index))
            leading_trailing_residues[1] = tuple(self.unfold_entities(loop_atoms, "R"))

        last_group = ss_atom_groups.get_group(ss_atom_groups.ngroups-1)
        if last_group.iloc[0] == "X":
            loop_atoms = tuple(self.get_atoms(include_atoms=last_group.index))
            leading_trailing_residues[len(ss_groups)] = tuple(self.unfold_entities(loop_atoms, "R"))

        number_ss = len(ss_groups)
        if verbose:
            print()
            for ss in ss_groups:
                ss_id = tuple(np.unique(r["residue_id"])[0] for r in ss)
                _ssid = list(sorted(map(int, ss_id)))
                print(_ssid[0], _ssid[-1], ss_type[ss_id])

        return ss_groups, loop_for_ss, original_order, ss_type, leading_trailing_residues, number_ss

    def remove_loops(self, verbose=False):
        raise NotImplementedError


class Structure(AbstractStructure):
    def __init__(self, path, name, input_format="pdb",feature_mode="r", features_path=None, residue_feature_mode="r", reset_chain=False, coarse_grained=False):
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

        if name is not None:
            self.name = name
            if len(name) == 7:
                self.cath_domain = cath_domain
                self.pdb = cath_domain[:4]
                self.chain = cath_domain[4]
                self.domNo = cath_domain[5:]
            else:
                self.pdb = name
                self.domNo = "00"
                reset_chain = True
        else:
            self.name = self.pdb = os.path.splitext(os.path.basename(path))[0]
            self.domNo = "00"
            reset_chain = True

        self.volume = volume

        try:
            self.structure = parser.get_structure(self.name, self.path)
        except KeyError:
            #Invalid mmcif file
            raise InvalidPDB("Invalid PDB file: {} (path={})".format(self.name, self.path))

        try:
            all_chains = list(self.structure[0].get_chains())
        except (KeyError, StopIteration):
            raise InvalidPDB("Error get chains for {} {}".format(self.name, self.path))

        if len(all_chains) > 1:
            raise InvalidPDB("Only accepts PDBs with 1 chain in {} {}".format(self.name, self.path))

        if reset_chain:
            self.chain = all_chains[0].id

        self.id = self.name #"{}{}{:02d}".format(self.pdb, self.chain, int(self.domNo))
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


        self.features = self.atom_features
        super().__init__(name, file_mode=feature_mode, features_path=features_path,
            coarse_grained=coarse_grained)

    def copy(self, empty=False):
        new = super().copy(empty=empty)

        #Make sure original files do not get overwritten
        new.features_path = os.getcwd()
        unique_id = int.from_bytes(os.urandom(4), sys.byteorder)
        new.atom_features_file = os.path.join(new.features_path , f"{os.path.basename(self.atom_features_file)}_{unique_id}.h5")
        new.residue_feature_file = os.path.join(new.features_path , f"{os.path.basename(self.atom_features_file)}_{unique_id}.h5")
        return new

    def deep_copy_feature(self, feature_name):
        if feature_name == "atom_features":
            return pd.DataFrame(
                np.nan if empty else np.array(copy.deepcopy(self.atom_features.values.tolist())),
                index=pd.Index(copy.deepcopy(self.atom_features.index.tolist())),
                columns=copy.deepcopy(self.atom_features.columns.tolist()))
        elif feature_name == "residue_features":
            return pd.DataFrame(
                np.nan if empty else np.array(copy.deepcopy(self.residue_features.values.tolist())),
                index=pd.Index(copy.deepcopy(self.residue_features.index.tolist())),
                columns=copy.deepcopy(self.residue_features.columns.tolist()))
        else:
            raise KeyError(feature_name)

    def __deepcopy__(self, memo):
        result = super().__deepcopy__(memo)

    def get_bfactors(self):
        return [a.bfactor for a in self.structure.get_atoms()]

    def get_atoms(self, atoms, include_hetatms=False, exclude_atoms=None, include_atoms=None):
        for a in atoms:
            hetflag, resseq, icode = a.get_parent().get_id()
            if not include_hetatms and hetflag != ' ':
                continue
            if exclude_atoms is not None and a.get_name().strip() in exclude_atoms:
                continue
            if include_atoms is not None and a.serial_number not in include_atoms:
                continue
            yield a

    def unfold_entities(self, entity_list, target_level="A"):
        return _unfold_entities(entity_list, target_level=level)

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

        path = os.path.join(work_dir, self.name)
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
        assert [len(f)==len(self.features) for f in features.values()]
        self.features = self.features.assign(**features)
        self.atom_feature_names += list(features.keys())

    def update_coords(self, coords):
        super().update_coords(coords)
        for atom, coord in zip(self.structure.get_atoms(), coords):
            atom.set_coord(coord)
