import os
from io import StringIO
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import spatial
from Bio import PDB
from Bio.PDB.NeighborSearch import NeighborSearch

from molmimic.common.Structure import Structure, number_of_features
from molmimic.common.ProteinTables import vdw_radii, vdw_aa_radii, surface_areas
from molmimic.common.featurizer import atom_feature_names, residue_feature_names

class ProteinVoxelizer(Structure):
    def __init__(self, path, cath_domain, input_format="pdb",
      volume=264, voxel_size=1.0, rotate=True, features_path=None,
      residue_feature_mode=None, use_features=None):
        super().__init__(path, cath_domain,
            input_format=input_format, feature_mode="r", features_path=features_path, residue_feature_mode=residue_feature_mode)

        self.mean_coord = np.zeros(3)
        self.mean_coord_updated = False

        self.volume = volume
        self.voxel_size = voxel_size
        self.voxel_tree = None
        self.atom_tree = None

        self.use_features = use_features

        if rotate is None:
            self.shift_coords_to_volume_center()
        else:
            next(self.rotate(rotate))

    def get_flat_features(self, resi=None):
        features = [s.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features) \
            for atom in s.structure.get_atoms()]

        return features

    def create_full_volume(self, input_shape=(96, 96, 96)):
        truth_grid = np.zeros(list(input_shape)+[1])
        for atom in self.get_atoms():
            grid = self.get_grid_coord(atom, vsize=input_shape[0])
            truth_grid[grid[0], grid[1], grid[2], 0] = 1
        return truth_grid

    def get_features_per_atom(residue_list):
        """Get features for eah atom, but not organized in grid"""
        features = [self.get_features_for_atom(a) for r in residue_list for a in r]
        return features

    def get_features(self, residue_list=None, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, expand_atom=False,
      undersample=False, autoencoder=False):
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

    def map_atoms_to_voxel_space(self, binding_site_residues=None,
      include_full_protein=True, only_aa=False, only_atom=False,
      use_deepsite_features=False, non_geom_features=False,
      only_surface=True, autoencoder=False, return_voxel_map=False, nClasses=2):
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
        assert [isinstance(binding_site_residues, (list, tuple)), autoencoder].count(True) == 1, \
            "Only binding_site_residues or autoencoder can be set"
        if binding_site_residues is not None:
            if not include_full_protein:
                atoms = sorted((a for r in binding_site_residues for a in r),
                    key=lambda a: a.get_serial_number())
                atoms = list(self.filter_atoms(atoms))
                binding_site_atoms = [a.get_serial_number() for a in atoms]
                non_binding_site_atoms = []
                nAtoms = len(atoms)
            else:
                atoms = list(self.get_atoms(include_hetatms=True))
                nAtoms = len(atoms)
                binding_site_atoms = [a.get_serial_number() for r in binding_site_residues for a in r]
                non_binding_site_atoms = []
                atoms = list(atoms)

        else:
            atoms = list(self.get_atoms(include_hetatms=True))
            nAtoms = len(atoms)
            binding_site_atoms = []
            non_binding_site_atoms = []


        nFeatures = number_of_features(
            only_aa=only_aa,
            only_atom=only_atom,
            non_geom_features=non_geom_features,
            use_deepsite_features=use_deepsite_features,
            course_grained=False)

        data_voxels = defaultdict(lambda: np.zeros(nFeatures if self.use_features is None else len(self.use_features)))
        truth_voxels = {}

        voxel_map = {}

        skipped = 0
        skipped_inside = []

        if nClasses == 2:
            true_value_ = np.array([0.,1.])
            neg_value_ = np.array([1.,0.])
        elif nClasses == 1:
            true_value_ = np.array([1.])
            neg_value_ = np.array([0.])
        elif nClasses == "sfams":
            pass

        for atom in atoms:
            atom = self._remove_altloc(atom)

            if autoencoder:
                truth = True
            elif binding_site_residues is None:
                truth = False
            else:
                truth = atom.get_serial_number() in binding_site_atoms

            if not truth and undersample and atom.get_serial_number() not in non_binding_site_atoms:
                skipped += 1
                continue

            if only_surface:
                features, is_buried = self.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features, warn_if_buried=True)
                if not truth and is_buried:
                    skipped += 1
                    skipped_inside.append((atom, features, is_buried))
                    continue
            else:
                features = self.get_features_for_atom(atom, only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features)

            features = features.astype(float)

            truth_value = true_value_.copy() if truth else neg_value_.copy()

            for atom_grid in self.get_grid_coords_for_atom_by_kdtree(atom):
                atom_grid = tuple(atom_grid.tolist())
                try:
                    data_value = np.maximum(features, data_voxels[atom_grid])
                    data_voxels[atom_grid] = data_value
                except ValueError:
                    print(nFeatures, data_voxels[atom_grid].shape, features.shape)
                    raise
                truth_voxels[atom_grid] = truth_value
                voxel_map[atom_grid] = atom.get_parent().get_id()


        #assert len(list(data_voxels)) > 0, (self.pdb, self.chain, self.sdi, data_voxels, list(data_voxels), np.array(list(data_voxels)), skipped, nAtoms, skipped_inside)

        outputs = None

        try:
            if binding_site_residues is None and not autoencoder:
                outputs = [np.array(list(data_voxels)), np.array(list(data_voxels.values()))]
            else:
                data = np.array(list(data_voxels.values())) #Always in same order as put in
                if autoencoder:
                    truth = data
                else:
                    truth = np.array([truth_voxels[grid] for grid in data_voxels])
                outputs = [np.array(list(data_voxels)), data, truth]
        except Exception as e:
            print(e)
            raise

        outputs.append(np.array(list(voxel_map.values())) if return_voxel_map else None)
        return outputs

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
                features = self.get_features_for_residue(residue, only_aa=only_aa, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features)
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
            truth = np.array([truth_voxels[grid] for grid in list(data_voxels.keys())])
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

    def get_features_for_atom(self, atom, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, warn_if_buried=False):
        """Calculate FEATUREs"""
#         from molmimic.parsers.Consurf import run_consurf
#         if not hasattr(self, "_consurf"):
#             self._consurf = run_consurf(self, self.pdb, self.chain)
#             print(self._consurf.keys())
#             assert 0, self._consurf.keys()

        # from molmimic.parsers.eppic import run_eppic
        # if not hasattr(self, "_eppic"):
        #     self._eppic = run_eppic(self.pdb, self.chain)

#         if not hasattr(self, "_infocontent"):
#             from Bio.Align import AlignInfo
#             from Bio.Align import MultipleSeqAlignment
#             self._align = MultipleSeqAlignment("~/SageMaker/cdd/fasta/cd00096.FASTA")
#             self._infocontent = AlignInfo(self._align).information_content()

        if isinstance(atom, PDB.Atom.DisorderedAtom):
            #All altlocs have been removed so onlt one remains
            atom = atom.disordered_get_list()[0]

        try:
            features = self.atom_features.loc[atom.serial_number]
            is_buried = bool(features["atom_is_buried"])


#             from Bio.SubsMat import MatrixInfo as matlist
#             matrix = matlist.blosum62
#             next(pairwise2.align.globaldx("KEVLA", "EVL", matrix))

            # try:
            #     con = np.array([self._eppic[atom.get_parent().get_id()]>0.6]).astype(float)
            # except Exception as e:
            #     assert 0, (atom.get_parent().get_id(), self._eppic.keys(), e)
            #     raise

            if only_aa and use_deepsite_features:
                feats = features[
                    atom_feature_names["get_residue"] + \
                    atom_feature_names["get_deepsite_features"] + \
                    atom_feature_names["get_charge_and_electrostatics"][1:3] +\
                    atom_feature_names["get_evolutionary_conservation_score"][-1:]]

                if warn_if_buried:
                    return feats, is_buried
                else:
                    return feats

            if use_deepsite_features:
                feats = features[
                    atom_feature_names["get_deepsite_features"] + \
                    atom_feature_names["get_charge_and_electrostatics"][1:3] +\
                    atom_feature_names["get_evolutionary_conservation_score"][-1:]]
                if warn_if_buried:
                    return feats, is_buried
                else:
                    return feats

            if only_atom:
                feats = features[atom_feature_names["get_atom_type"]]
                if warn_if_buried:
                    return feats, is_buried
                else:
                    return feats

            elif only_aa:
                feats = features[atom_feature_names["get_residue"]]
                if warn_if_buried:
                    return feats, is_buried
                else:
                    return feats
            elif non_geom_features:
                feats = features[
                    atom_feature_names["get_element_type"] + \
                    atom_feature_names["get_charge_and_electrostatics"] +\
                    atom_feature_names["get_hydrophobicity"] +\
                    [float(is_buried)]]
                if warn_if_buried:
                    return feats, is_buried
                else:
                    return feats
            elif self.use_features is not None:
                feats = features[self.use_feats]
                if warn_if_buried:
                    return feats, is_buried
                else:
                    return feats
            else:
                print("ALL FEATS")
                if warn_if_buried:
                    return features, is_buried
                else:
                    return features
        except ValueError as e:
            # print e
            # pass
            raise

    def get_features_for_residue(self, residue, only_aa=False, non_geom_features=False,
      use_deepsite_features=False):
        """Calculate FEATUREs"""
        try:
            features = self.residue_features[residue.get_id()[1]-1]
            if non_geom_features:
                return np.concatenate((
                    features[15:36],
                    features[0:4],
                    features[8:12],
                    ))
            elif only_aa:
                return features[15:36]
            else:
                return features
        except ValueError:
            pass

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

    def rotate(self, rvs=None, num=1):
        for r, M in super().rotate(rvs=rvs, num=num):
            self.set_voxel_size(self.voxel_size)
            yield r, M


    def set_voxel_size(self, voxel_size=None):
        self.voxel_size = voxel_size or 1.0

        coords = self.get_coords()
        min_coord = np.floor(np.min(coords, axis=0))-5
        max_coord = np.ceil(np.max(coords, axis=0))+5
        extent_x = np.arange(min_coord[0], max_coord[0], self.voxel_size)
        extent_y = np.arange(min_coord[1], max_coord[1], self.voxel_size)
        extent_z = np.arange(min_coord[2], max_coord[2], self.voxel_size)
        mx, my, mz = np.meshgrid(extent_x, extent_y, extent_z)
        self.voxel_tree = spatial.cKDTree(list(zip(mx.ravel(), my.ravel(), mz.ravel())))

    def convert_voxels(self, grid, radius=2.75, level="A"):
        """Convert grid points to atoms
        """
        if self.atom_tree is None:
            self.atom_tree = NeighborSearch(list(self.structure.get_atoms()))

        return self.atom_tree.search(grid, radius, level=level)

def voxels_to_unknown_structure(coords, data, bucket_size=10):
    from Bio.PDB.kdtrees import KDTree
    assert bucket_size > 1
    assert self.coords.shape[1] == 3
    kdt = KDTree(coords, bucket_size)
    neighbors_max = kdt.neighbor_search(1.55) #2.75
    neighbors_min = kdt.neighbor_search(1.53) #1.3
    
    neighbors_min = [sorted((n1.index1, n2.index2)) for n in neighbors_min]
    
    neighbor_atoms = defaultdict(list)
    
    for neighbor in neighbors_max:
        i1 = neighbor.index1
        i2 = neighbor.index2
        
        if sorted((i1, i2)) in neighbors_min:
            continue
        
        voxel1 = data[i1]
        voxel2 = data[i2]
        
        atom1_isC = voxel1["hydrophic"] and not voxel1["aromatic"]
        atom2_isC = voxel2["hydrophic"] and not voxel2["aromatic"]
        
        if atom1_isC and atom2_isC:
            #Could a CA, CB, or CG?
            pass
        
        if atomtype_1 == "C" and atomtype_2 == "C":
            pass

    #Refold using rosetta CA-CA??

def voxels_to_structure_from_starting_structure(coords, data, starting_structure, starting_data, features=None, diff=True, prefix=None, work_dir=None):
    reconstructed_features = []
    for atom in starting_structure.get_atoms():
        atom = self._remove_altloc(atom)
        
        atom_orig_feats = []
        atom_reconstructed_feats = []
        for atom_grid in self.get_grid_coords_for_atom_by_kdtree(atom):
            idx = starting_data.index(atom_grid)
            atom_orig_feats.append(starting_data[idx])
            atom_reconstructed_feats.append(coords[atom_grid])
        
        atom_orig_feats = np.avarage(np.concatenate(atom_orig_feats))
        atom_reconstructed_feats = np.avarage(np.concatenate(atom_reconstructed_feats))
        
        reconstructed_features_ = atom_reconstructed_feats-atom_orig_feats if diff else atom_reconstructed_feats
        reconstructed_features.append(reconstructed_features_)
    
    reconstructed_features = np.concatenate(reconstructed_features, axis=0)
    
    if features is None:
        features = range(starting_data.shape[1])
        
    if prefix is None:
        prefix = "reconstructed_{}".format(starting_structure.cathdomain)
        
    if work_dir is None:
        work_dir = os.getcwd()
        
    outfiles = []
    
    for i, feature in enumerate(features):
        starting_structure.update_bfactors(features[:, i])
        outfile = os.path.join(work_dir, "{}_{}.pdb".format(prefix, feature))
        starting_structure.write_pdb(outfile)
        outfiles.append(outfile)
        
    return outfiles    

if __name__ == "__main__":
  import sys
  assert len(sys.argv) > 1
  structure = Structure(sys.argv[1], panfs=False).extract_chain("A")
  num_atoms = sum(1 for a in structure.structure.get_atoms() if a.get_parent().get_id()[0] == " ")
  for atom in structure.get_atoms():
      print(atom.get_full_id(), atom.serial_number, "of", num_atoms, features, len(features))
      features = structure.get_features_for_atom(atom)
