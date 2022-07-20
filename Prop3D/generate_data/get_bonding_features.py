def get_bonding_features(self):
    voxel_structure = Voxelizer(self.pdb, self.cath_domain, input_format=self.input_format,
      volume=self.volume, voxel_size=self.voxel_size, rotate=False, features_path=self.features_path,
      use_features=self.use_features, ligand=self.ligand)
    coord_atoms = {}
    coordsFile = os.path.join(self.work_dir, f"{self.cath_domain}.bondedVoxels.csv")
    with open(coordsFile, "w") as f:
        for a1, a2, voxels in voxel_structure.get_overlapping_voxels():
            for v in voxels:
                print(f"{v[0]},{v[1]},{v[2]}", file=f)
                coord_atoms.append((a1.get_full_id(), a2.get_full_id()))

    electrostatics = Electrostatics(work_dir=self.work_dir, job=self.job)
    vertex_charges = electrostatics.get_electrostatics_at_coordinates_from_pdb(
        coordsFile, self.path, force_field="parse", noopt=True, apbs_input=True,
        whitespace=True)

    for (atom1, atom2), charge in zip(coord_atoms, vertex_charges):
        pass
