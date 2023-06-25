from collections import defaultdict

import numpy as np
from scipy import spatial
import numpy.lib.recfunctions

from Prop3D.common.DistributedStructure import DistributedStructure
from Prop3D.common.ProteinTables import vdw_aa_radii
from Prop3D.common.features import default_atom_features, default_residue_features

class DistributedVoxelizedStructure(DistributedStructure):
    def __init__(self, path, key, cath_domain_dataset, coarse_grained=False,
      volume=264, voxel_size=1.0, rotate=None, use_features=None, predict_features=None,
      replace_na=False, ligand=False):
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
    
    def create_full_volume(self, input_shape=(96, 96, 96)):
        truth_grid = np.zeros(list(input_shape)+[1])
        for atom in self.get_atoms():
            for grid in self.get_vdw_grid_coords_for_atom(atom["X", "Y", "Z"]):
                truth_grid[grid[0], grid[1], grid[2], 0] = 1
        return truth_grid

    def shift_coords_to_volume_center(self):
        return self.shift_coords(np.array([self.volume/2]*3))

    def resize_volume(self, new_volume, shift=True):
        self.volume = new_volume
        if shift:
            self.shift_coords_to_volume_center()

    def rotate(self, rvs=None, num=1, return_to=None):
        if return_to is None:
            return_to=[self.volume/2]*3
        for r in super().rotate(rvs=rvs, num=num, return_to=return_to):
            self.set_voxel_size(self.voxel_size)
            yield r

    def orient_to_pai(self, random_flip=False, flip_axis=(0.2, 0.2, 0.2)):
        super().orient_to_pai(random_flip=random_flip, flip_axis=flip_axis)
        self.shift_coords_to_volume_center()

    def get_features_per_atom(self, residue_list):
        """Get features for eah atom, but not organized in grid"""
        return self.data[self.data[:,"residue_id"].isin(residue_list)]

    def get_features(self, residue_list=None, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, expand_atom=False,
      undersample=False, autoencoder=False):
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

    def map_atoms_to_voxel_space(self, truth_residues=None,
      only_surface=False, autoencoder=False, return_voxel_map=False,
      return_serial=False, return_b=False, nClasses=2, simple_fft=None,
      verbose=False, use_raw_atom_coords=False):
        """Map atoms to sparse voxel space.

        Parameters
        ----------
        truth_residues : list of Bio.PDB.Residue objects or None
            If a binding is known, add the list of Bio.PDB.Residue objects, usually
            obtained by Structure.align_seq_to_struc()
        include_full_protein : boolean
            If true, all atoms from the protein are used. Else, just the atoms from the
            defined binding site. Only makes sense if truth_residues is not None
        Returns
        -------
        indices : np.array((nVoxels,3))
        data : np.array((nVoxels,nFeatures))
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
                data[feature][ind] = default_atom_features[feature]

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
                truth_value = atom[self.self.predict_features]
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

            if truth_residues is not None and not autoencoder:
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

    def map_residues_to_voxel_space(self, truth_residues=None, only_surface=True,
      autoencoder=False, return_voxel_map=False, return_serial=False, return_b=False,
      nClasses=2, simple_fft=None, verbose=False):
        return map_atoms_to_voxel_space(self, truth_residues=truth_residues,
            only_surface=only_surface, autoencoder=autoencoder,
            return_voxel_map=return_voxel_map, return_serial=return_serial,
            return_b=return_b, nClasses=nClasses, simple_fft=simple_fft,
            verbose=verbose)

    def simple_fft_scoring_features(self, atom_or_residue, mode="simple", b=3):
        """Rp=−1  on a surface layer and Rp=1 on the core of the receptor,
        Lp=1 on the entire ligand, and Rp=Lp=0 everywhere else. It is clear that
        this scoring function, which is essentially the one used by
        Katchalski-Katzir et al. (5), reaches its minimum on a conformation in
        which the ligand maximally overlaps with the surface layer of the receptor,
        thus providing optimal shape complementarity. https://doi.org/10.1073/pnas.1603929113"""

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

    def get_vdw_grid_coords_for_atom(self, atom, atom_index):
        dist = self.get_vdw(atom)
        coord = np.around(self.coords[atom_index], decimals=4)
        neighbors = self.voxel_tree.query_ball_point(coord, r=dist)
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_closest_grid_coord_for_atom(self, atom):
        _, neighbors = self.voxel_tree.query([atom.coord])
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_vdw_grid_coords_for_residue(self, residue):
        dist = vdw_aa_radii.get(residue.get_resname(), 3.2)
        center = np.nanmean([a.get_coord() for a in residue], axis=0)
        neighbors = self.voxel_tree.query_ball_point(center, r=dist)
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    def get_closest_grid_coord_for_residue(self, residue):
        center = np.nanmean([a.get_coord() for a in residue], axis=0)
        _, neighbors = self.voxel_tree.query([center])
        for idx in neighbors:
            yield self.voxel_tree.data[idx]

    # def rotate(self, rvs=None, num=1):
    #     for r, M in super().rotate(rvs=rvs, num=num):
    #         self.set_voxel_size(self.voxel_size)
    #         yield r, M

    def resize_volume(self, new_volume, shift=True):
        super().resize_volume(new_volume, shift=shift)
        self.set_voxel_size(self.voxel_size)

    def set_voxel_size(self, voxel_size=None, full_grid=True):
        self.voxel_size = voxel_size or 1.0

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

    def convert_voxels(self, grid, radius=2.75, level="A"):
        """Convert grid points to atoms
        """
        if self.atom_tree is None:
            self.atom_tree = spatial.cKDTree(list(self.get_atoms()))

        idx = self.atom_tree.query_ball_point(grid, radius)
        return self.data[idx]

    def get_overlapping_voxels(self):
        neighbor_atoms = self.calculate_neighbors(d_cutoff=5.0, level="A")
        for a1, a2 in neighbor_atoms:
            v1 = set(self.get_vdw_grid_coords_for_atom(a1))
            v2 = set(self.get_vdw_grid_coords_for_atom(a2))
            overlap = v1.intersection(v2)
            yield a1, a2, overlap
