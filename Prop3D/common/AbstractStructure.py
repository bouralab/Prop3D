import os
import copy
from sys import stdin
from typing import Any, Union, IO, AnyStr, TypeVar, Optional

from collections.abc import Iterator

import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA
from scipy.stats import special_ortho_group

from Prop3D.common.ProteinTables import three_to_one
from Prop3D.util import natural_keys

_Self = TypeVar('_Self', bound='AbstractStructure')

class AbstractStructure(object):
    """A base structure class that holds default methods for dealing with structures.
    Subclasses can be created to use different protein libraires, e.g. BioPython, BioTite,
    our own HDF/HSDS distributed stucture.

    Parameters
    ----------
    name: str
        Name of proten, e.g. PDB or CATH id
    file_mode: str (r, w, r+, w+, etc)
        Open file for reading or writing. Defualt is just reading, no methods will affect underlying file
    coarse_grained: boolean
        Use a residue only model instead of an all atom model. Defualt False. Warning, not fully implemented.
    """
    def __init__(self, name: str, file_mode: str = "r", coarse_grained: bool = False) -> None:
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

    def copy(self, empty: bool = False) -> _Self:
        """Create a deep copy of current structure.

        Parameters
        ----------
        (deprecated) empty: bool
            Don't copy features
        """
        new = copy.deepcopy(self)
        return new

    def deep_copy_feature(self, feature_name: str) -> Any:
        """Subclass this method to handle custom copying of specific features

        Parameters
        ----------
        feature_name: str
            Feature name to copy

        Raises
        ------
        NotImeplementedError if no method to handle feature
        """
        raise NotImplementedError

    def __deepcopy__(self, memo: dict) -> Any:
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

    def __abs__(self) -> _Self:
        """Take the absolue value of all features
        """
        new = self.copy()
        new.features = new.features.abs()
        return new

    def __sub__(self, other: _Self) -> _Self:
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

    def __add__(self, other: _Self) -> _Self:
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

    def __floordiv__(self, other: _Self) -> _Self:
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

    def __truediv__(self, other: _Self) -> _Self:
        """Divide feature values from other structure from this stucture.
        """
        return self.__floordiv__(other)

    def normalize_features(self, columns: Union[str, list[str]] = None) -> _Self:
        """Normalize features using min max scaling

        Parameters
        ----------
        columns: str or list of strs
            Names of feature columns to normalize

        Returns
        -------
        A copy of this AbstractStructure with normalized features in the dataframe
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

    def get_atoms(self, include_hetatms: bool = False, exclude_atoms: bool = None, include_atoms: bool = None) -> Iterator[Any]:
        """Subclass to enumerate protein model for all atoms with options to filter

        Parameters
        ----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        raise NotImplementedError

    def filter_atoms(atoms: list[Any] = None, include_hetatms: bool = False, exclude_atoms: list[Any] = None, include_atoms: list[Any] = None)  -> Iterator[Any]:
        """Subclass to enumerate protein model for all atoms with options to filter

        Parameters
        ----------
        include_hetatms : boolean
            Inclue hetero atoms or not. Default is False.
        exlude_atoms : list
            List of atoms to skip during enumeration (depends on model if id or pyton object)
        inlude_atoms : list
            List of atoms to inllude during enumeration (depends on model if id or pyton object)
        """
        yield from self.get_atoms(atoms, include_hetatms=include_hetatms,
            exclude_atoms=exclude_atoms, include_atoms=include_atoms)

    def get_surface(self) -> Any:
        """Returns all surface atoms, using DSSP accessible surface value"
        """
        surface = self.features[self.features["residue_buried"]==0]
        return surface

    def get_bfactors(self) -> Any:
        """Get bfactors for all atoms
        """
        raise NotImplementedError

    def save_pdb(self, path: Union[str, None] = None, header: Union[str, list[str], None] = None, file_like: bool = False, rewind: bool = True) -> Union[str, IO[AnyStr]]:
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
        raise NotImplementedError

    def write_features(self, features: Union[str, list[str], None], coarse_grained: bool = False, name: Union[str, None] = None, work_dir: Union[str, None] = None) -> None:
        """Subclass to write features to a spefic file, depnignd on protein loading class, e.g. HDF

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
        raise NotImplementedError

    def write_features_to_pdb(self, features_to_use: Union[str, list[str], None], name: Union[str, None] = None, coarse_grain: bool = False, work_dir: Union[str, None] = None, other: Any = None):
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

    def add_features(self, coarse_grained: bool = False, **features):
        """Subclass to add a feature column to dataset
        """
        raise NotImplementedError

    def get_coords(self, include_hetatms: bool = False, exclude_atoms: list[Any] = None) -> Any:
        """Subclass to return all XYZ coordinates fro each atom

        Parameters
        ----------
        include_hetatms : bool
            Include heteroatoms or not. Default is False.
        exclude_atoms : list of atoms
            Atoms to exlude while getting coordinates

        Returns
        -------
        XYZ coordinates from each atom in a specified format of subclass
        """
        raise NotImplementedError
    
    def get_sequence(self) -> str:
        """Get amino acid sequence from structure
        """
        sequence = ""
        prev_resi = None
        for a in self.get_atoms():
            resi = a["residue_id"]
            if resi != prev_resi:
                resn = three_to_one(a["residue_name"].decode("utf-8"))
                sequence += resn
                prev_resi = resi
        
        return sequence

    def update_coords(self, coords: np.array) -> None:
        """Sublcass to create method to update XYZ coordinates with a new set of coordinates for the same atoms"""
        self.coords = coords
        self.mean_coord = None
        self.mean_coord_updated = False

    def get_mean_coord(self) -> np.array:
        """Get the mean XYZ coordinate or center of mass.
        """
        if not self.mean_coord_updated:
            self.mean_coord = np.around(np.nanmean(self.coords, axis=0), decimals=4)
            self.mean_coord_updated = True
        return self.mean_coord

    def get_max_coord(self) -> np.array:
        """Get the maximum coordinate in each dimanesion
        """
        return np.nanmax(self.coords, axis=0)

    def get_min_coord(self) -> np.array:
        """Get the minimum coordinate in each dimanesion
        """
        return np.nanmin(self.coords, axis=0)

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
        coords = self.coords
        if from_origin or new_center is None:
            mean_coord = self.get_mean_coord()
            coords -= mean_coord
        if new_center is not None:
            coords += np.around(new_center, decimals=4)

        self.update_coords(coords)
        self.mean_coord_updated = False
        return np.nanmean(coords, axis=0)

    def shift_coords_to_origin(self) -> float:
        """Center structure at the origin

        Returns
        -------
        The new center coordinate
        """
        return self.shift_coords()

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

    def rotate(self, rvs: Union[np.array, None] = None, num: int = 1, return_to: Union[tuple[float], np.array, None] = None) -> Iterator[tuple[int, np.array]]:
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
                M = special_ortho_group.rvs(3)
                #print("rotating randomly")
            else:
                M=rvs
                #assert 0
            self.shift_coords_to_origin()
            old_coords = self.get_coords()
            coords = np.dot(self.coords, M).round(decimals=4)
            self.update_coords(coords)

            self.shift_coords(self.get_mean_coord() if return_to is None else return_to)
            
            yield r, M

    def update_bfactors(self, b_factors: list[Any]) -> None:
        """Sublcass to create method to update bfactors with a new set of bfactors for the same atoms"""
        raise NotImplementedError

    def calculate_neighbors(self, d_cutoff: float = 100.0) -> Any:
        """Subclass to find all nearest neighbors within a given radius.

        Parameters
        ----------
        d_cutoff : float
            Distance cutoff to find neighbors. Deualt is 100 Angtroms
        """
        raise NotImplementedError

    def get_vdw(self, element_name: str, residue: bool = False) -> float:
        """Get van der walls radii for an atom or residue
        """
        if not residue:
            return np.array([vdw_radii.get(element_name.title(), 1.7)])
        else:
            #For coarse graining, not really a vdw radius
            return np.array([vdw_aa_radii.get(element_name, 3.0)])

        raise NotImplementedError

    def get_dihedral_angles(self, atom_or_residue: Any) -> float:
        """Get deidral angle for atom (mapped up to residue) or the residue
        """
        raise NotImplementedError

    def get_secondary_structures_groups(self, verbose: bool = False, ss_min_len: int = 3, is_ob: bool=False, assume_correct: Optional[pd.Series]=None) -> tuple[list[pd.DataFrame], dict[tuple[str], pd.DataFrame], dict[tuple[str], int], dict[tuple[str], str], dict[int, list[Any]], int]:
        """Get groups of adjecent atom rows the belong to the same secondary structure.
        We use DSSP to assing secondary structures to each reisdue mapped down to atoms. 
        If any segment was <4 residues, they were merged the previous and next groups

        Returns
        -------
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

        if assume_correct is None:
            ss_groups = ss_type.groupby((ss_type != ss_type.shift()).cumsum()-1)
            ob_has_h = False
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
                    if len(this_group_residues)<3 and prev_group == next_group and not is_ob:
                        ss_type.loc[ss_group.index] = prev_group
                    elif len(this_group_residues)<3 and next_group == "X" and this_group != prev_group:
                        ss_type.loc[ss_group.index] = "X"

                    if this_group=="H" and prev_group=="E" and next_group=="E" and len(this_group_residues)<5:
                        ss_type.loc[ss_group.index] = "X"

                    if (not is_ob and len(this_group_residues)<ss_min_len) or (is_ob and len(this_group_residues)<5): #OB=5
                        if prev_group == next_group:
                            ss_type.loc[ss_group.index] = prev_group
                        else:
                            pass

                        if this_group=="H" and prev_group=="E" and next_group=="E":
                            ss_type.loc[ss_group.index] = "X"
                        elif this_group=="E" and prev_group=="H" and next_group=="H":
                            ss_type.loc[ss_group.index] = "X"

                    if this_group=="X" and (len(this_group_residues)>10 or (is_ob and not ob_has_h and len(this_group_residues)>5)):
                        ss_type.loc[ss_group.index] = "H"
                        ob_has_h = True
        else:
            ss_type = assume_correct

        #Regroup with correct SS
        ss_atom_groups = ss_type.groupby([(ss_type != ss_type.shift()).cumsum()-1])
        import pdb; pdb.set_trace()
        ss_groups = []
        loop_for_ss = {}
        original_order = {}
        ss_type_final = {}
        leading_trailing_residues = {}

        prev_ss_id = None
        for i, ss_group in ss_atom_groups:
            #Get all atoms from SS and loops
            ss_atoms = tuple(self.get_atoms(include_atoms=ss_group.index))
            ss_residues = tuple(self.unfold_entities(ss_atoms, "R"))
            ss_residues_id = tuple(np.unique(r["residue_id"])[0] for r in ss_residues)

            if ss_group.iloc[0] != "X":
                ss_groups.append(ss_residues)
                original_order[ss_residues_id] = len(ss_groups)
                ss_type_final[ss_residues_id] = ss_group.iloc[0]
                prev_ss_id = ss_residues_id
            elif prev_ss_id is not None and len(ss_groups)>0 and ss_group.iloc[0] == "X":
                #loop_for_ss[ss_residues_id] = ss_residues
                loop_for_ss[prev_ss_id] = ss_residues

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
                print(_ssid[0], _ssid[-1], ss_type_final[ss_id])

        should_continue = None
        while not isinstance(should_continue, bool):
            should_continue = input("Is this correct? [Y/n]")
            if should_continue.lower() in ["", "y", "n"]:
                should_continue = not (should_continue == "n")
        

        if not should_continue:
            print("Enter correct segments per line (start stop type '1 9 E'):")
            #segment_resi = []
            ss_type.loc[:] = "X"
            for line in stdin:
                try:
                    start, stop, ss_type_str = line.rstrip().split()
                except Exception:
                    if line.rstrip() == "":
                        break
                    else:
                        raise RuntimeError("Must be (start stop type '1 9 E')")
                
                start_key = natural_keys(start, use_int=True)
                stop_key = natural_keys(stop, use_int=True)

                assert stop_key>start_key, (start_key, stop_key)

                all_keys = [natural_keys(x.decode("ascii"), use_int=True) for x in self.data["residue_id"]]
                use_idx = [i for i, k in enumerate(all_keys) if start_key<=k<=stop_key]
                
                if ss_type_str == "E":
                    ss_type.loc[self.data[use_idx]["serial_number"]] = "E"
                    # self.data[use_idx]["is_sheet"] = 1.
                    # self.data[use_idx]["is_helix"] = 0.
                    # self.data[use_idx]["Unk_SS"] = 0.
                    #segment_resi += use_idx
                elif ss_type_str == "H":
                    ss_type.loc[self.data[use_idx]["serial_number"]] = "H"
                    # self.data[use_idx]["is_sheet"] = 0.
                    # self.data[use_idx]["is_helix"] = 1.
                    # self.data[use_idx]["Unk_SS"] = 0.
                    #segment_resi += use_idx
                else:
                    ss_type.loc[self.data[use_idx]["serial_number"]] = "X"
                    # self.data[use_idx]["is_sheet"] = 0.
                    # self.data[use_idx]["is_helix"] = 0.
                    # self.data[use_idx]["Unk_SS"] = 1.
            
            # loop_resi = list(sorted(set(range(len(self.data)))-set(segment_resi)))
            # self.data[loop_resi]["is_sheet"] = 0.
            # self.data[loop_resi]["is_helix"] = 0.
            # self.data[loop_resi]["Unk_SS"] = 1.

            # import pdb; pdb.set_trace()
        
            return self.get_secondary_structures_groups(verbose=verbose, ss_min_len=ss_min_len, is_ob=is_ob, assume_correct=ss_type)

        return ss_groups, loop_for_ss, original_order, ss_type_final, leading_trailing_residues, number_ss

    def remove_loops(self, verbose: bool = False) -> None:
        raise NotImplementedError
