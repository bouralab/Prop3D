from collections import defaultdict
from typing import Union
from collections.abc import Iterator

import numpy as np
from scipy import spatial
import numpy.lib.recfunctions

from Prop3D.common.DistributedStructure import DistributedStructure
from Prop3D.common.ProteinTables import vdw_aa_radii
from Prop3D.common.features import all_features

class DistributedVoxelizedStructure(DistributedStructure):
    """A structure class to deal with structures originating from a distributed
    HSDS instance. This class also handles proteins in voxelized volumes

    Parameters
    ----------
    path : str
        Path to h5 file in HSDS endpoint
    key : str
        Key to access speficic protein inside the HDF file
    cath_domain_dataset : str
        The CATH superfamily if endpoint is setup to use CATH (use '/' instead of '.')
    coarse_grained : boolean
        Use a residue only model instead of an all atom model. Defualt False. Warning, not fully implemented.
    volume : float
        Size in Angstroms^3 of the entire volume. Defualt is 264. 
    voxel_size : float
        Resolution to map atoms to inside volume. Default is 1.0.
    rotate : None, bool, np.array
        Rotate moleule randomly (set to True) or using a given roation matrix. If None, no rotations will be perfomed. Defualt None.
    use_features : list of str
        Features to include for trianing a model.
    predict_features : list of str
    replace_na : bool
        Replace "Not a number" values with defualt values. Defualt is False.
    ligand : bool
        Not used often. Only used in simple_fft_scoring_features for specyng if protein is interacting partner. Defualt is False.
    """
    def __init__(self, path: str, key: str, cath_domain_dataset: str, coarse_grained: bool = False,
      volume: float = 264., voxel_size: float = 1.0, rotate: Union[bool, np.array, None] = None, use_features: Union[list[str], None] = None, predict_features: Union[list[str], None] = None,
      replace_na: bool = False, ligand: bool = False) -> None:
        super().__init__(path, key, cath_domain_dataset, coarse_grained=coarse_grained)

        self.mean_coord = np.zeros(3)
        self.mean_coord_updated = False

        self.volume = volume
        self.voxel_size = voxel_size
        self.voxel_tree = None
        self.atom_tree = None

        self.use_features = use_features if use_features is not None else self.feature_names
        self.predict_features = predict_features

        if self.predict_features is not None and use_features is not None:
            assert len(set(self.predict_features).intersection(set(self.use_features)))==0, \
                "Cannot train on and predict the same features"

        self.replace_na = replace_na

        self.ligand = ligand

        self.features = self.features[self.use_features]

        if rotate is None or (isinstance(rotate, bool) and not rotate):
            self.shift_coords_to_volume_center()
            self.set_voxel_size(self.voxel_size)
        elif isinstance(rotate, str) and rotate == "pai":
            self.orient_to_pai()
            self.shift_coords_to_volume_center()
            self.set_voxel_size(self.voxel_size)
        elif (isinstance(rotate, str) and rotate == "random") or (isinstance(rotate, bool) and rotate):
            next(self.rotate())
        elif isinstance(rotate, np.ndarray):
            next(self.rotate(rotate))
        else:
            raise RuntimeError("Invalid rotation option. Must be None or False for no rotation, 'pai' to orient to princple axis, 'random' for random rotation matrix, or an actual roation matrix")

            # if r is None:
            #     setattr(result, k, copy.deepcopy(v, memo))

            # if "features" in k:
            #     setattr(result, k, self.deep_copy_feature(k))
            # else:
            #     print(k, v)
            #     setattr(result, k, copy.deepcopy(v, memo))
    
    def create_full_volume(self, input_shape: Union[np.array, list[int], None] = None) -> np.array:
        """Create a dense representation of the protein

        Parameters
        ----------
        input_shape : 3-tuple
            New volume size. If None, use given volume size. Defualt is None.
        
        Returns
        -------
        A dense grid with the protein in the center
        """
        if input_shape is None:
            input_shape = [self.volme]*3

        truth_grid = np.zeros(list(input_shape)+[1])
        for atom in self.get_atoms():
            for grid in self.get_vdw_grid_coords_for_atom(atom[["X", "Y", "Z"]]):
                truth_grid[grid[0], grid[1], grid[2], 0] = 1
        return truth_grid

    def shift_coords_to_volume_center(self) -> np.array:
        """Shift coordinatesto the center of the volume

        Returns
        -------
        The new center coordinate
        """
        return self.shift_coords(np.array([self.volume/2]*3))

    def resize_volume(self, new_volume: float, shift: bool = True) -> None:
        """Increase or decrease the volume and move the protein to its new center

        Parameters
        ----------
        new_volume : 3-tuple
            New volume size. If None, use given volume size. Defualt is None
        shift : bool
            Move protein to volume's new center. Defualt is True
        """
        self.volume = new_volume
        if shift:
            self.shift_coords_to_volume_center()
        self.set_voxel_size(self.voxel_size)

    def rotate(self, rvs: Union[np.array, None] = None, num: int = 1, return_to: Union[tuple[float], np.array, None] = None) -> Iterator[tuple[int, np.array]]:
        """Rotate structure by either randomly in place or with a set rotation matrix. 
        Random rotations matrices are drawn from the Haar distribution (the only uniform 
        distribution on SO(3)) from scipy.

        This method automatically recenters the protein in the current volume size.

        Parameters
        ----------
        rvs : np.array (3x3)
            A rotation matrix. If None, a randome roation matrix is used. Default is None.
        num : int
            Number of rotations to perfom
        return_to : None or XYZ coordinate
            When finsihed rotating, move structure to this coordinate. Defualt is to the center volume center

        Yields
        ------
        r : int
            Rotation number
        M : np.array (3x3)
            Rotation matrix
        """
        if return_to is None:
            return_to=[self.volume/2]*3
        for r in super().rotate(rvs=rvs, num=num, return_to=return_to):
            self.set_voxel_size(self.voxel_size)
            yield r

    def orient_to_pai(self, random_flip: bool = False, flip_axis: Union[tuple[float], np.array] = (0.2, 0.2, 0.2)) -> None:
        """Orient structure to the Principle Axis of Interertia and optionally flip. Modified from EnzyNet.

        This method automatically recenters the protein in the current volume size.

        Parameters
        ----------
        random_flip : bool
            Randomly flip around axis. Defualt is False.
        flip_axis: 3-tuple of floats
            New axis to flip.
        """
        super().orient_to_pai(random_flip=random_flip, flip_axis=flip_axis)
        self.shift_coords_to_volume_center()

    def get_features_per_atom(self, residue_list: list[str]) -> np.array:
        """Get features for each atom, but not organized in grid
        
        Parameters
        ----------
        residue_list : list
            A list of residue ids
        
        Returns
        -------
        A numpy array of atoms in selected residues
        """
        return self.data[self.data[:,"residue_id"].isin(residue_list)]

    def get_features(self, residue_list: Union[list[str], None] = None, only_aa: bool = False, only_atom: bool = False,
      non_geom_features: bool = False, use_deepsite_features: bool = False, expand_atom: bool = False,
      undersample: bool = False, autoencoder: bool = False):
        """Deprecated. not used. See map_atoms_to_voxel_space and map_residues_to_voxel_space
        """
        if self.coarse_grained:
            return self.map_residues_to_voxel_space(
                truth_residues=residue_list,
                include_full_protein=include_full_protein,
                only_aa=only_aa,
                non_geom_features=non_geom_features,
                undersample=undersample
            )
        return self.map_atoms_to_voxel_space(
            expand_atom=expand_atom,
            truth_residues=residue_list,
            include_full_protein=include_full_protein,
            only_aa=only_aa,
            only_atom=only_atom,
            use_deepsite_features=use_deepsite_features,
            non_geom_features=non_geom_features,
            undersample=undersample)

    def map_atoms_to_voxel_space(self, truth_residues: Union[list[str], None] = None,
      only_surface: bool = False, autoencoder: bool = False, return_voxel_map: bool = False,
      return_serial: bool = False, return_b: bool = False, nClasses: int = 2, simple_fft: Union[str, None] = None,
      verbose: bool = False, use_raw_atom_coords: bool = False) -> tuple[
          np.array, 
          np.array, 
          np.array, 
          dict[int, list[tuple[int,int,int]]], 
          dict[tuple[int,int,int], list[int]],
          dict[tuple[int,int,int], float]
          ]:
        """Map atoms to sparse voxel space.

        Parameters
        ----------
        truth_residues : list of residue_ids or None
            If a binding site (or other site of interest) is known, add the list of residue_ids
        only_surface : bool
            Only return voxels for atoms present on the surface of the protein. Default False.
        autoencoder : bool
            Use same features for input and output. Default is False.
        return_voxel_map : bool
            Return a mapping to back from voxels to atoms. Defualt is False.
        reutrn_serial : bool
            Return atoms serials present in each voxel. Defualt is False.
        return_b : bool
            Return all bfactors used in each voxel. Default is False.
        nClasses : int [1,2]
            Rerely used. If only predicited one feature, create nClasses number of features. E.g. for 2: a feature would be [is_false, is_true]
        simple_fft : None or str ["simple", "zdock"]
            Type of fft features to use. Defualt is None.
        verbose : bool
            Print out logs while running. Defualt is False.
        use_raw_atom_coords : bool
            Instead of mapping atoms to voxel, attach features to raw atom coords. Default is False.

        Returns
        -------
        coords : np.array((nVoxels,3))
        feats : np.array((nVoxels,nFeatures))
        truth : np.array((nVoxels,nFeatures))
        voxel_map : Dictionary of voxels to atoms
        serial : Dictionary of serials to voxels
        b : Dictionary of voxels to b factors
        """
        assert not self.coarse_grained, "Cannot be used with the coarse graned model"
        assert [isinstance(truth_residues, (list, tuple)), autoencoder, isinstance(self.predict_features, (list, tuple))].count(True) == 1, \
            "Only truth_residues or autoencoder can be set"

        if truth_residues is not None:
            predicting_features = False
        else:
            predicting_features = isinstance(self.predict_features, (list, tuple))

        data_voxels = defaultdict(lambda: np.zeros(len(self.use_features)))
        truth_voxels = {}

        voxel_map = {}

        b_factors_voxels = {}
        serial_number_voxels = defaultdict(list)

        skipped = 0
        skipped_inside = []

        if nClasses == 2:
            true_value_ = np.array([0.,1.])
            neg_value_ = np.array([1.,0.])
        elif nClasses == 1:
            true_value_ = np.array([1.])
            neg_value_ = np.array([0.])
        elif nClasses == "sfams":
            raise RuntimeError("Sfams not implemented")
        else:
            true_value_ = np.array([1.])
            neg_value_ = np.array([0.])

        data = self.data #[self.use_features]

        if self.replace_na:
            for feature in self.use_features:
                ind = data[feature] == np.nan
                data[feature][ind] = all_features.default_atom_features[feature]

        for atom_index in range(len(self.data)):
            atom = data[atom_index]

            if only_surface and atom["residue_buried"]==1:
                continue

            if autoencoder or predicting_features:
                truth = True
            elif truth_residues is None:
                truth = False
            else:
                truth = atom["residue_id"] in truth_residues

            features = atom[self.use_features]

            features = numpy.lib.recfunctions.structured_to_unstructured(features)

            if simple_fft is not None:
                features = self.simple_fft_scoring_features(atom, mode=simple_fft)

            #Handle truth values if its not an autoencoder
            if predicting_features:
                truth_value = atom[self.predict_features].tolist()
            elif truth_residues is not None:
                truth_value = true_value_.copy() if truth else neg_value_.copy()

            if use_raw_atom_coords:
                grid_coords = [tuple(self.coords[atom_index])]
            else:
                grid_coords = [tuple(g) for g in self.get_vdw_grid_coords_for_atom(atom, atom_index)]

            voxel_map[atom["serial_number"]] = grid_coords

            for atom_grid in grid_coords:
                try:
                    data_value = np.maximum(features, data_voxels[atom_grid])
                    data_voxels[atom_grid] = data_value
                except ValueError:
                    print(data_voxels[atom_grid].shape, features.shape)
                    raise
                if not autoencoder:
                    truth_voxels[atom_grid] = np.maximum(
                        truth_value, truth_voxels.get(atom_grid, truth_value))

                b_factors_voxels[atom_grid] = np.maximum(
                    atom["bfactor"], b_factors_voxels.get(atom_grid, 0))
                serial_number_voxels[atom_grid].append(atom["serial_number"])

        outputs = None

        try:
            coords, feats = zip(*data_voxels.items())
            outputs = [np.array(coords), np.array(feats)]

            if (truth_residues is not None or predicting_features) and not autoencoder:
                truth = np.array([truth_voxels[grid] for grid in coords])
            else:
                truth = None

            outputs.append(truth)

        except Exception as e:
            print(e)
            raise

        if return_voxel_map:
            outputs.append(voxel_map)
        else:
            outputs.append(None)

        if return_serial:
            outputs.append([serial_number_voxels[grid] for grid in coords])
        else:
            outputs.append(None)

        if return_b:
            outputs.append(np.array([b_factors_voxels[grid] for grid in coords]))

        return outputs

    def map_residues_to_voxel_space(self, truth_residues: Union[list[str], None] = None,
      only_surface: bool = False, autoencoder: bool = False, return_voxel_map: bool = False,
      return_serial: bool = False, return_b: bool = False, nClasses: int = 2, simple_fft: Union[str, None] = None,
      verbose: bool = False) -> None:
        raise NotImplementedError
        return self.map_atoms_to_voxel_space(self, truth_residues=truth_residues,
            only_surface=only_surface, autoencoder=autoencoder,
            return_voxel_map=return_voxel_map, return_serial=return_serial,
            return_b=return_b, nClasses=nClasses, simple_fft=simple_fft,
            verbose=verbose)

    def simple_fft_scoring_features(self, atom_or_residue: Union[int, str], mode: str = "simple", b: int = 3) -> np.array:
        """If voxelized proteins will be used in FFT docking type algorithms, use these features.
        
        Rp=−1  on a surface layer and Rp=1 on the core of the receptor,
        Lp=1 on the entire ligand, and Rp=Lp=0 everywhere else. It is clear that
        this scoring function, which is essentially the one used by
        Katchalski-Katzir et al. (5), reaches its minimum on a conformation in
        which the ligand maximally overlaps with the surface layer of the receptor,
        thus providing optimal shape complementarity. https://doi.org/10.1073/pnas.1603929113
        
        Parameters
        ----------
        atom_or_residue : int
            Index of atom or residue in data table
        mode : str or bool ["simple", "zdock", True]
            Type of features to use. Defualt is simple.
        b : int
            Not used.
        """

        if not self.coarse_grained:
            residue_buried = self.atom_features[atom_or_residue, "residue_rasa"]<0.5
            charge = self.features[atom_or_residue, "charge"]
            electrostatic_potential = self.features[atom_or_residue, "electrostatic_potential"]
        else:
            residue_buried = self.features[atom_or_residue, "residue_buried"]
            charge = self.features[atom_or_residue, "charge"]
            electrostatic_potential = self.eatures[atom_or_residue, "electrostatic_potential"]

        if mode in [True, "simple"]:
            if not self.ligand:
                if residue_buried:
                    features = np.array([-15, 0])
                else:
                    features = np.array([1, 0])
            else:
                features = np.array([1, 0])
                # if residue_buried:
                #     features = np.array([-15, 0])
                # else:
                #     features = np.array([1, 0])

            return features

        elif mode == "zdock":
            psc_elec = np.array([
                3.5**2 if residue_buried else 3.5, #Same for ligand
                charge if not self.ligand else 0
            ])

            return psc_elec

        else:
            print("Mode is", mode, mode in [True])

    def get_vdw_grid_coords_for_atom(self, atom: np.array, atom_index: int) -> Iterator[list[tuple[int, int, int]]]:
        """Get all grid coordinates for atom that itersect with its van der walls volume

        Parameters
        ----------
        atom : np.array
            Full atom features from 1 row of the entire feature dateset. Must include columns X,Y,Z
        atom_index : int
            deprecated.

        Yeilds
        -------
        Grid coordinates for atom
        """
        dist = self.get_vdw(atom)
        coord = np.around(atom[["X", "Y", "Z"]].tolist(), decimals=4)
        neighbors = self.voxel_tree.query_ball_point(coord, r=dist)
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_closest_grid_coord_for_atom(self, atom: np.array) -> Iterator[list[tuple[int, int, int]]]:
        coord = np.around(atom[["X", "Y", "Z"]].tolist(), decimals=4)
        _, neighbors = self.voxel_tree.query([coord])
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_vdw_grid_coords_for_residue(self, residue: np.array) -> Iterator[list[tuple[int, int, int]]]:
        dist = vdw_aa_radii.get(residue["residue_name"], 3.2)
        coords = [np.around(a[["X", "Y", "Z"]].tolist(), decimals=4) for a in residue]
        center = np.nanmean(coords, axis=0)
        neighbors = self.voxel_tree.query_ball_point(center, r=dist)
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_closest_grid_coord_for_residue(self, residue: np.array) -> Iterator[list[tuple[int, int, int]]]:
        """Yields all grid points that are near any atom inside the given residue
        """
        center = np.nanmean([a[["X", "Y", "Z"]] for a in residue], axis=0)
        _, neighbors = self.voxel_tree.query([center])
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def set_voxel_size(self, voxel_size: float = 1.0, full_grid: bool = True) -> None:
        """Set the voxel size or resultion of the mapping the structure to the volume.

        This method calculates the entire grid used to search over while mapping atoms.

        Parameters
        ----------
        voxel_size : float
            New voxel size in angstroms
        full_grid : bool
            Create coords for the entire volume. If False, only save grid point closest to the protein. Default True.
        """
        self.voxel_size = voxel_size

        coords = self.get_coords()
        min_coord = np.floor(np.nanmin(coords, axis=0))-5
        max_coord = np.ceil(np.nanmax(coords, axis=0))+5
        max_dist = np.linalg.norm(max_coord-min_coord)

        if full_grid:
            min_coord_ = np.zeros(3)
            max_coord_ = np.array([self.volume]*3)
            max_dist_ = np.linalg.norm(max_coord_-min_coord_)

            fail = False
            if np.any(min_coord<min_coord_):
                print(f"Min coordinate outside grid: {min_coord} < {min_coord_}")
                fail=True
            if np.any(max_coord>=max_coord_):
                print(f"Max coordinate outside grid: {max_coord} > {max_coord_}")
                fail=True
            assert not fail
            max_coord = max_coord
            min_coord = min_coord

        extent_x = np.arange(min_coord[0], max_coord[0], self.voxel_size)
        extent_y = np.arange(min_coord[1], max_coord[1], self.voxel_size)
        extent_z = np.arange(min_coord[2], max_coord[2], self.voxel_size)
        mx, my, mz = np.meshgrid(extent_x, extent_y, extent_z)

        self.voxel_tree = spatial.cKDTree(list(zip(mx.ravel(), my.ravel(), mz.ravel())))
        #spatial.cKDTree(self.get_coords())

    def convert_voxels(self, grid: np.array, radius: float = 2.75, level: str = "A") -> Union[np.array, list[np.array]]:
        """Convert grid points to atoms

        Parameters
        ----------
        grid : 3-tuple
            grid point
        radius : float
            Find atoms within a certain radius in angstroms from the grid point. Default 2.7.
        level : str
            If 'R', map atoms back to residues
        """
        if not hasattr(self, "atom_tree") or self.atom_tree is None:
            self.atom_tree = spatial.cKDTree(self.coords)

        idx = self.atom_tree.query_ball_point(grid, radius)

        if level == "R":
            return list(self.unfold_entities(self.data[idx]))

        return self.data[idx]

    def get_overlapping_voxels(self) -> Iterator[tuple[list[tuple[int, int, int]],list[tuple[int, int, int], list[tuple[int, int, int]]]]]:
        """For all pairs on interacting atoms (<5 Angstroms), get the overlapping voxels

        Yields
        ------
        a1 : int
            atom serial 1
        a2 : int
            atom serial 2
        overlap : list
            list of overlapping voxels coorindates    
        """
        neighbor_atoms = self.calculate_neighbors(d_cutoff=5.0, level="A")
        for a1, a2 in neighbor_atoms:
            v1 = set(self.get_vdw_grid_coords_for_atom(a1))
            v2 = set(self.get_vdw_grid_coords_for_atom(a2))
            overlap = v1.intersection(v2)
            yield a1, a2, overlap
