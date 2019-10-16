import os

import pandas as pd
import numpy as np
from sklearn.gaussian_process.kernels import RBF
from Bio import PDB
import freesasa

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from toil.realtimeLogger import RealtimeLogger

from molmimic.generate_data.util import InvalidPDB, natural_keys, silence_stdout, silence_stderr
from molmimic.parsers.FreeSASA import run_freesasa_biopython
from molmimic.parsers.Electrostatics import run_pdb2pqr_APBS
from molmimic.parsers.CX import run_cx
from molmimic.parsers.Consurf import run_consurf
from molmimic.parsers.eppic import run_eppic
from molmimic.parsers.DSSP import run_dssp
from molmimic.parsers.PDBQT import get_autodock_features

from molmimic.common.Structure import Structure, number_of_features, \
    get_feature_names, angle_between, get_dihedral
from molmimic.common.ProteinTables import surface_areas, hydrophobicity_scales

atom_feature_names = {
    "get_atom_type": [
        "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F",
        "MG", "P", "SA", "S", "CL", "CA", "MN", "FE", "ZN", "BR", "I", "Unk_atom"],
    "get_element_type": [
        "C_elem", "N_elem", "O_elem", "S_elem", "H_elem", "F_elem",
        "MG_elem", "P_elem", "S_elem", "CL_elem", "CA_elem", "MN_elem",
        "FE_elem", "ZN_elem", "BR_elem", "I_elem", "Unk_elem"],
    "get_vdw": ["vdw_volume"],
    "get_charge_and_electrostatics": [
        "charge", "neg_charge", "pos_charge", "neutral_charge",
        "electrostatic_potential", "is_electropositive", "is_electronegative"],
    "get_concavity": ["cx", "is_concave", "is_convex", "is_concave_and_convex"],
    "get_hydrophobicity": [
        "hydrophobicity", "is_hydrophbic", "is_hydrophilic", "hydrophobicity_is_0"],
    "get_accessible_surface_area": [
        "atom_asa", "atom_is_buried", "atom_exposed",
        "residue_asa", "residue_buried", "residue_exposed"],
    "get_residue": PDB.Polypeptide.aa3+["Unk_residue"],
    "get_ss": ["is_helix", "is_sheet", "Unk_SS"],
    "get_deepsite_features": [ #DeepSite features
        "hydrophobic_atom",
        "aromatic_atom",
        "hbond_acceptor",
        "hbond_donor",
        "metal"],
    "get_evolutionary_conservation_score": [
        "consurf_conservation_normalized",
        "consurf_conservation_scale",
        "consurf_is_conserved",
        "eppic_conservation",
        "eppic_is_conserved",
        "is_conserved"]
}
atom_feature_names_flat = [col_name for _, col_names for atom_feature_names.items() for col_name in col_names]

residue_feature_names = {
    "get_charge_and_electrostatics": [
        "charge", "neg_charge", "pos_charge", "neutral_charge",
        "electrostatic_potential", "is_electropositive", "is_electronegative"],
    "get_concavity": [
        "cx", "is_concave", "is_convex", "is_concave_and_convex"],
    "get_hydrophobicity": [
        "hydrophobicity", "is_hydrophbic", "is_hydrophilic", "hydrophobicity_is_0"],
    "get_accessible_surface_area": [
        "residue_asa", "residue_buried", "residue_exposed"],
    "get_residue": PDB.Polypeptide.aa3+["Unk_residue"],
    "get_ss": ["is_helix", "is_sheet", "Unk_SS"],
    "get_evolutionary_conservation_score": [
        "consurf_conservation_normalized",
        "consurf_conservation_scale",
        "consurf_is_conserved",
        "eppic_conservation",
        "eppic_is_conserved",
        "is_conserved"]
}
residue_feature_names_flat = [col_name for _, col_names for residue_feature_names.items() for col_name in col_names]


class ProteinFeaturizer(Structure):
    def __init__(self, path, cath_domain, job, work_dir,
      input_format="pdb", force_feature_calculation=False, update_features=None):
        feature_mode = "w+" if force_feature_calculation else "r"
        super(ProteinFeaturizer, self).__init__(
            path, cath_domain,
            input_format=input_format,
            feature_mode=feature_mode,
            residue_feature_mode=feature_mode,
            features_path=work_dir)
        self.job = job
        self.work_dir = work_dir
        self.update_features = update_features

    def calculate_flat_features(self, course_grained=False, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False):
        if course_grained:
            features = [self.calculate_features_for_residue(
                self._remove_inscodes(r), only_aa=only_aa,
                non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for r in self.structure.get_residues()]
            if self.residue_feature_mode == "w+" or self.update_features is not None:
                self.write_features(course_grained=True)
            return features, self.residue_features_file
        else:
            features = [self.calculate_features_for_atom(
                self._remove_altloc(atom), only_aa=only_aa,
                only_atom=only_atom, non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for atom in self.structure.get_atoms()]
            if self.atom_feature_mode == "w+" or self.update_features is not None:
                self.write_features()
            return features, self.atom_features_file

    def write_features(self, course_grained=False):
        if course_grained:
            self.residue_features = self.residue_features.astype(np.float64)
            self.residue_features = self.residue_features.drop(columns=
                [col for col in self.residue_features if c not in \
                residue_feature_names_flat])
            self.residue_features.to_hdf(self.residue_features_file, "table")
        else:
            self.atom_features = self.atom_features.astype(np.float64)
            self.atom_features = self.atom_features.drop(columns=
                [col for col in self.atom_features if c not in \
                atom_feature_names_flat])
            self.atom_features.to_hdf(self.atom_features_file, "table")

    def get_features_per_atom(residue_list):
        """Get features for eah atom, but not organized in grid"""
        features = [self.get_features_for_atom(self._remove_altloc(a)) for r in residue_list for a in r]
        return features

    def get_features_per_residue(self, residue_list):
        features = [self.get_features_for_residue(self._remove_inscodes(r)) for r in residue_list]
        return features

    def calculate_features_for_atom(self, atom, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, warn_if_buried=False):
        if self.update_features is not None:
            for feat_type, feat_names in atom_feature_names.items():
                if feat_type in self.update_features:
                    getattr(self, feat_type)(atom)
                else:
                    use_feats = [feat_name in self.update_features for feat_name \
                        in feat_names]
                    if any(use_feats):
                        getattr(self, feat_type)(atom)
            if warn_if_buried:
                if "is_buried" not in self.atom_features.columns:
                    is_burried = self.get_accessible_surface_area(atom, save=False)
                else:
                    is_buried = self.atom_features.loc[atom.serial_number, "is_buried"]
                return self.atom_features, bool(is_buried["is_buried"])
            else:
                return self.atom_features

        if use_deepsite_features:
            self.get_deepsite_features(atom)
            if warn_if_buried:
                is_burried = self.get_accessible_surface_area(atom, save=False)
        elif only_atom:
            self.get_element_type(atom)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(atom, save=False)
        elif only_aa:
            self.get_residue(atom)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(atom, save=False)
        elif non_geom_features:
            self.get_element_type(atom)
            self.get_charge_and_electrostatics(atom)
            self.get_hydrophobicity(atom)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(atom, save=False)
        else:
            self.get_atom_type(atom)
            self.get_element_type(atom)
            self.get_vdw(atom)
            self.get_charge_and_electrostatics(atom)
            self.get_concavity(atom)
            self.get_hydrophobicity(atom)
            self.get_accessible_surface_area(atom)
            self.get_residue(atom)
            self.get_ss(atom)
            self.get_deepsite_features(atom, calc_charge=False,
                calc_conservation=False)
            self.get_evolutionary_conservation_score(atom)

            is_buried = self.atom_features.loc[atom.serial_number, "atom_is_buried"]

        # RealtimeLogger.info("Finished atom {} {}".format(atom, atom.serial_number))
        # RealtimeLogger.info("Feats {}".format(features))
        # RealtimeLogger.info("Feats shape {}".format(features.shape))
        # RealtimeLogger.info("atom_features shape {}".format(self.atom_features.shape))

        #self.atom_features.loc[atom.serial_number, :] = features

        #RealtimeLogger.info("atom_features is {}".format(self.atom_features))

        if warn_if_buried:
            return self.atom_features, bool(is_buried["atom_is_buried"])
        else:
            return self.atom_features

    def calculate_features_for_residue(self, residue, only_aa=False, non_geom_features=False,
      use_deepsite_features=False, warn_if_buried=False):
        """Calculate FEATUREs"""
        if self.update_features is not None:
            for feat_type, feat_names in residue_feature_names.items():
                if feat_type in self.update_features:
                    getattr(self, feat_type)(atom)
                else:
                    use_feats = [feat_name in self.update_features for feat_name \
                        in feat_names]
                    if any(use_feats):
                        getattr(self, feat_type)(atom)
            if warn_if_buried:
                if "is_buried" not in self.residue_features.columns:
                    is_burried = self.get_accessible_surface_area(residue, save=False)
                else:
                    is_buried = self.residue_features.loc[residue.get_id(), "is_buried"]
                return self.residue_features, bool(is_buried["is_buried"])
            else:
                return self.atom_features

        if non_geom_features:
            self.get_residue(residue)
            self.get_charge_and_electrostatics(residue)
            self.get_hydrophobicity(residue)
            self.get_evolutionary_conservation_score(residue)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(residue, save=False)
        elif only_aa:
            self.get_residue(residue)
            if warn_if_buried:
                is_buried = self.get_accessible_surface_area(residue, save=False)
        else:
            self.get_charge_and_electrostatics(residue)
            self.get_concavity(residue)
            self.get_hydrophobicity(residue)
            self.get_accessible_surface_area(residue)
            self.get_residue(residue)
            self.get_ss(residue)
            self.get_evolutionary_conservation_score(residue)
            is_buried = self.residue_features.loc[[residue.get_id()], "residue_buried"]

        if warn_if_buried:
            return self.residue_features, bool(is_buried["residue_buried"])
        else:
            return self.residue_features

    def get_atom_type(self, atom):
        """Get Autodock atom type"""

        if not hasattr(self, "_autodock"):
            self._autodock = get_autodock_features(self.path, work_dir=self.work_dir, job=self.job)

        try:
            atom_type, h_bond_donor = self._autodock[atom.serial_number]
        except KeyError:
            atom_type = "Unk_atom"

        if atom_type == "  ":
            atom_type = "Unk_atom"

        try:
            self.atom_features.loc[atom.serial_number, atom.get_name().strip()] = 1.
        except KeyError:
            self.atom_features.loc[atom.serial_number, "Unk_atom"] = 1.

        return self.atom_features.loc[atom.serial_number,
            atom_feature_names["get_atom_type"]]

    def get_element_type(self, atom):
        elems = "CNOS"
        elem_col = "{}_elem".format(atom.element)

        if elem_col in atom_feature_names["get_element_type"]:
            self.atom_features.loc[atom.serial_number, elem_col] = 1.
        else:
            self.atom_features.loc[atom.serial_number, "Unk_elem"] = 1.

        return self.atom_features.loc[atom.serial_number,
            atom_feature_names["get_element_type"]]

    def get_charge_and_electrostatics(self, atom_or_residue, only_charge=False,
      only_bool=False):
        if self.update_features is not None:
            if "electrostatic_potential" not in self.update_features and \
              "is_electropositive" not in self.update_features and \
              "is_electronegative" not in self.update_features:
                 only_charge = True

        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            atoms = pd.concat([self.get_charge_and_electrostatics(
                self._remove_altloc(a)) for a in residue], axis=0)
            charge_value = atoms["charge"].sum()
            electrostatic_pot_value = atoms["electrostatic_potential"].sum()

            charge = [
                charge_value,
                float(charge_value < 0),
                float(charge_value > 0),
                float(charge_value == 0)]
            cols = residue_feature_names["get_charge_and_electrostatics"][:4]
            if not only_charge:
                charge += [
                    electrostatic_pot_value,
                    float(electrostatic_pot_value < 0),
                    float(electrostatic_pot_value > 0)]
                cols += residue_feature_names["get_charge_and_electrostatics"][4:]

            if only_bool:
                charge = charge[1:3]
                cols = cols[1:3]
                if not only_charge:
                    charge += charge[5:7]
                    cols += charge[5:7]

            idx = [residue.get_id()]
            self.residue_features.loc[idx, cols] = charge
            return self.residue_features.loc[idx, cols]
        else:
            raise RuntimeError("Input must be Atom or Residue: {}".format(type(atom_or_residue)))

        if not hasattr(self, "_pqr"):
            try:
                self._pqr = run_pdb2pqr_APBS(self.path, pdb2pqr_whitespace=True,
                    work_dir=self.work_dir, job=self.job, only_charge=only_charge)
            except (SystemExit, KeyboardInterrupt):
                raise
            except Exception as e:
                self._pqr = {}
                RealtimeLogger.info("ELECTROSTATICS failed ({}): {}".format(type(e), e))

        atom_id = atom.get_full_id()[3:5]

        if atom_id[1][1] != " ":
            #pdb2pqr removes alternate conformations and only uses the first
            atom_id = (atom_id[0], (atom_id[1][0], " "))

        try:
            charge_value, electrostatic_pot_value = self._pqr[atom_id]
        except KeyError:
            charge_value, electrostatic_pot_value = np.NaN, np.NaN

        charge = [
            charge_value,
            float(charge_value < 0),
            float(charge_value > 0),
            float(charge_value == 0)]
        cols = atom_feature_names["get_charge_and_electrostatics"][:4]
        if not only_charge:
            charge += [
                electrostatic_pot_value,
                float(electrostatic_pot_value < 0),
                float(electrostatic_pot_value > 0)]
            cols += atom_feature_names["get_charge_and_electrostatics"][4:]

        if only_bool:
            charge = charge[1:3]
            cols = cols[1:3]
            if not only_charge:
                charge += charge[5:7]
                cols += charge[5:7]

        charge = np.array(charge)
        idx = atom.serial_number
        self.atom_features.loc[idx, cols] = charge

        return self.atom_features.loc[idx, cols]

    def get_concavity(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            concavity_value = pd.concat([self.get_concavity(
                self._remove_altloc(a)) for a in residue], axis=0)["cx"].mean()
            #concave, convex, or both
            concavity = np.array([
                concavity_value,
                float(concavity_value <= 2),
                float(concavity_value > 5),
                float(2 < concavity_value <= 5)])

            idx = [residue.get_id()]
            cols = residue_feature_names["get_concavity"]
            self.residue_features.loc[idx, cols] = concavity

            return self.residue_features.loc[idx, cols]
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if not hasattr(self, "_cx"):
            self._cx = run_cx(self.path, work_dir=self.work_dir, job=self.job)

        concavity_value = self._cx.get(atom.serial_number, np.NaN)
        #concave, convex, or both
        concavity = np.array([
            concavity_value,
            float(concavity_value <= 2),
            float(concavity_value > 5),
            float(2 < concavity_value <= 5)])

        idx = atom.serial_number
        cols = atom_feature_names["get_concavity"]
        self.atom_features.loc[idx, cols] = concavity

        return self.atom_features.loc[idx, cols]

    def get_hydrophobicity(self, atom_or_residue, scale="kd"):
        assert scale in list(hydrophobicity_scales.keys())
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
            use_atom = True
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            use_atom = False
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        try:
            resname = PDB.Polypeptide.three_to_one(residue.get_resname())
            hydrophobicity_value = hydrophobicity_scales[scale][resname]
        except KeyError:
            hydrophobicity_value = np.nan

        hydrophobicity = np.array([
            hydrophobicity_value,
            float(hydrophobicity_value < 0),
            float(hydrophobicity_value > 0),
            float(hydrophobicity_value == 0)])

        if use_atom:
            idx = atom.serial_number
            cols = atom_feature_names["get_hydrophobicity"]
            self.atom_features.loc[idx, cols] = hydrophobicity
            return self.atom_features.loc[idx, cols]
        else:
            idx = [residue.get_id()]
            cols = residue_feature_names["get_hydrophobicity"]
            self.residue_features.loc[idx, cols] = hydrophobicity
            return self.residue_features.loc[idx, cols]

    def get_accessible_surface_area(self, atom_or_residue, save=True):
        """Returns the ASA value from freesasa (if inout is Atom) and the DSSP
        value (if input is Atom or Residue)

        Returns
        -------
        If input is residue a 3-vector is returned, otherwise a 4-vector is returned
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            return pd.concat((
                self.get_accessible_surface_area_atom(atom_or_residue, save=save),
                self.get_accessible_surface_area_residue(atom_or_residue, save=save)),
                axis=0)
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            return self.get_accessible_surface_area_residue(atom_or_residue, save=save)
        else:
            raise RuntimeError("Input must be Atom or Residue")

    def get_accessible_surface_area_atom(self, atom_or_residue, save=True):
        """Returns the ASA value from freesasa (if inout is Atom) and the DSSP
        value (if input is Atom or Residue)

        Returns
        -------
        If input is residue a 3-vector is returned, otherwise a 4-vector is returned
        """
        if not isinstance(atom_or_residue, PDB.Atom.Atom):
            raise RuntimeErorr("Input must be Atom")

        atom = atom_or_residue
        residue = atom_or_residue.get_parent()

        if not hasattr(self, "_sasa"):
            self._sasa = run_freesasa_biopython(self.path)

        sasa, sasa_struct = self._sasa

        total_area = surface_areas.get(atom.element.title(), 1.0)
        try:
            selection = "sele, chain {} and resi {} and name {}".format(self.chain, atom.get_parent().get_id()[1], atom.get_id()[0])
            with silence_stdout(), silence_stderr():
                selections = freesasa.selectArea([selection], sasa_struct, sasa)
                atom_area = selections["sele"]
            fraction = atom_area/total_area
        except (KeyError, AssertionError, AttributeError, TypeError):
            raise
            atom_area = np.NaN
            fraction = np.NaN

        asa = np.array([
            atom_area,
            float(fraction <= 0.2), #buried
            float(fraction > 0.2) #exposed
        ])
        idx = atom.serial_number
        cols = atom_feature_names["get_accessible_surface_area"][:3]

        if save:
            self.atom_features.loc[idx, cols] = asa
            return self.atom_features.loc[idx, cols]
        else:
            return pd.Series(asa, index=cols)

    def get_accessible_surface_area_residue(self, atom_or_residue, save=True):
        """Returns the ASA value from freesasa (if inout is Atom) and the DSSP
        value (if input is Atom or Residue)

        Returns
        -------
        If input is residue a 3-vector is returned, otherwise a 4-vector is returned
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeError("Input must be Atom or Residue")

        if not hasattr(self, "_dssp"):
            self._dssp = run_dssp(self.structure, self.path,
                work_dir=self.work_dir, job=self.job)

        try:
            residue_area = float(self._dssp[residue.get_full_id()[2:]][3])
        except (KeyError, AssertionError, AttributeError, TypeError, ValueError):
            try:
                #Remove HETATMs
                residue_area = float(self._dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][3])
            except (KeyError, AssertionError, AttributeError, TypeError, ValueError):
                residue_area = np.NaN

        asa = np.array([
            residue_area,
            float(residue_area < 0.2),
            float(residue_area >= 0.2)])

        if is_atom:
            idx = atom.serial_number
            cols = atom_feature_names["get_accessible_surface_area"][3:]
            self.atom_features.loc[idx, cols] = asa
            if save:
                return self.atom_features.loc[idx, cols]
            else:
                return pd.Series(asa, index=cols)
        else:
            idx = [residue.get_id()]
            cols = residue_feature_names["get_accessible_surface_area"]
            self.residue_features.loc[idx, cols] = asa
            if save:
                return self.residue_features.loc[idx, cols]
            else:
                return pd.Series(asa, index=cols)

    def get_residue(self, atom_or_residue):
        """
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        # residues = [0.]*(len(PDB.Polypeptide.aa1)+1) #one hot residue
        #
        # try:
        #     residues[PDB.Polypeptide.three_to_index(residue.get_resname())] = 1.
        # except (ValueError, KeyError) as e:
        #     residues[-1] = 1.
        # return np.array(residues)

        col = residue.get_resname() if residue.get_resname() in PDB.Polypeptide.aa3 \
            else "Unk_element"

        if is_atom:
            self.atom_features.loc[atom.serial_number, col] = 1.
            return self.atom_features.loc[atom.serial_number,
                atom_feature_names["get_residue"]]
        else:
            self.residue_features.loc[[residue.get_id()], col] = 1.
            return self.residue_features.loc[[residue.get_id()],
                residue_feature_names["get_residue"]]

    def get_ss(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if not hasattr(self, "_dssp"):
            self._dssp = run_dssp(self.structure, self.path,
                work_dir=self.work_dir, job=self.job)

        try:
            atom_ss = self._dssp[residue.get_full_id()[2:]][2]
        except (KeyError, AssertionError, AttributeError, TypeError):
            try:
                #Remove HETATMs
                atom_ss = self._dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][2]
            except (KeyError, AssertionError, AttributeError, TypeError):
                atom_ss = "X"

        ss = np.array([
            float(atom_ss in "GHIT"),
            float(atom_ss in "BE"),
            float(atom_ss not in "GHITBE")])

        if is_atom:
            idx = atom.serial_number
            cols = atom_feature_names["get_ss"]
            self.atom_features.loc[idx, cols] = ss
            return self.atom_features.loc[idx, cols]
        else:
            idx = [residue.get_id()]
            cols = residue_feature_names["get_ss"]
            self.residue_features.loc[idx, cols] = ss
            return self.residue_features.loc[idx, cols]

    def get_deepsite_features(self, atom, calc_charge=True, calc_conservation=True):
        """Use DeepSites rules for autodock atom types
        """
        if not hasattr(self, "_autodock"):
            self._autodock = get_autodock_features(self.path, work_dir=self.work_dir, job=self.job)

        try:
            atom_type, h_bond_donor = self._autodock[atom.serial_number]
        except KeyError:
            atom_type = "  "

        idx = atom.serial_number
        cols = atom_feature_names["get_deepsite_features"]

        features = np.array([ #zeros(nFeatures, dtype=bool)
            #hydrophobic
            (atom_type == 'C') | (atom_type == 'A'),
            #aromatic
            atom_type == 'A',
            #hbond_acceptor
            (atom_type == 'NA') | (atom_type == 'NS') | (atom_type == 'OA') | \
            (atom_type == 'OS') | (atom_type == 'SA'),
            #hbond_donor
            (atom_type == 'HS') | (atom_type == 'HD'),
            #metal
            (atom_type == 'MG') | (atom_type == 'ZN') | (atom_type == 'MN') | \
            (atom_type == 'CA') | (atom_type == 'FE')], dtype=float)

        self.atom_features.loc[idx, cols] = features

        if calc_charge:
            charge = self.get_charge_and_electrostatics(atom, only_charge=True, only_bool=True)
            cols += charge.columns.tolist()

        if calc_conservation:
            cons = self.get_evolutionary_conservation_score(atom, only_bool=True)
            cols += cons.columns.tolist()

        return self.atom_features.loc[idx, cols]

    def get_evolutionary_conservation_score(self, atom_or_residue, consurf=True,
      eppic=True, only_bool=False):
        if self.update_features is not None:
            if "is_conserved" in self.update_features:
                consurf = eppic = True
            else:
                if "consurf_conservation_normalized" in self.update_features or \
                  "consurf_conservation_scale" in self.update_features or \
                  "consurf_is_conserved" in self.update_features:
                    consurf = True
                else:
                    consurf = False

                if "eppic_conservation" in self.update_features or \
                  "eppic_is_conserved" in self.update_features:
                    eppic = True
                else:
                    eppic = False

        if consurf and not hasattr(self, "_consurf"):
            self._consurf = run_consurf(self.pdb, self.chain, work_dir=self.work_dir)

        if eppic and not hasattr(self, "_eppic"):
            self._eppic = run_eppic(self.pdb, self.chain, work_dir=self.work_dir)
            #RealtimeLogger.info("EPPIC is {}".format(self._eppic))

        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeError("Input must be Atom or Residue")

        idx = residue.get_id()

        consurf_cols = []
        consurf_is_conserved = 0
        if consurf:
            try:
                consurf_normalized_score, consurf_conservation_score = \
                    self._consurf[residue.get_id()]
            except KeyError as e:
                consurf_normalized_score = consurf_conservation_score = np.NaN

            consurf_cols = ["consurf_conservation_normalized",
                "consurf_conservation_scale", "consurf_is_conserved"]
            consurf_is_conserved = float(consurf_conservation_score>=6)
            consurf_value = np.array([
                consurf_normalized_score, consurf_conservation_score,
                consurf_is_conserved
            ])

            if not only_bool:
                if is_atom:
                    idx = atom.serial_number
                    self.atom_features.loc[idx, consurf_cols] = consurf_value
                else:
                    idx = [residue.get_id()]
                    self.residue_features.loc[idx, consurf_cols] = consurf_value

        eppic_cols = []
        eppic_is_conserved = 0
        if eppic:
            try:
                eppic_conservation = self._eppic[residue.get_id()]
            except Exception as e:
                eppic_conservation = np.NaN
                RealtimeLogger.info("EPPIC key not found {} {}".format(e, self._eppic.keys()))

            eppic_cols = ["eppic_conservation", "eppic_is_conserved"]
            eppic_is_conserved = float(eppic_conservation<0.5)
            eppic_value = np.array([eppic_conservation, eppic_is_conserved])

            if not only_bool:
                if is_atom:
                    idx = atom.serial_number
                    self.atom_features.loc[idx, eppic_cols] = eppic_value
                else:
                    idx = residue.get_id()
                    self.residue_features.loc[idx, eppic_cols] = eppic_value

        if is_atom:
            idx = atom.serial_number
            is_conserved = float(self.atom_features.loc[idx, ["consurf_is_conserved",
                "eppic_is_conserved"]].sum() > 0) if not only_bool else \
                float(consurf_is_conserved+eppic_is_conserved > 0)
            self.atom_features.loc[idx, "is_conserved"] = is_conserved
            return self.atom_features.loc[idx, consurf_cols+eppic_cols+["is_conserved"]]
        else:
            idx = [residue.get_id()]
            is_conserved = self.residue_features.loc[idx, ["consurf_is_conserved",
                "eppic_is_conserved"]].sum() > 0 if not only_bool else \
                float(consurf_is_conserved+eppic_is_conserved > 0)
            self.residue_features.loc[idx, "is_conserved"] = is_conserved
            return self.residue_features.loc[idx, consurf_cols+eppic_cols+["is_conserved"]]

    def calculate_graph(self, d_cutoff=100.):
        import networkx as nx
        structure_graph = nx.Graph()
        for r1, r2 in self.calculate_neighbors(d_cutoff=d_cutoff):
            structure_graph.add_edge(r1.get_id(), r2.get_id(),
                attr_dict=self.get_edge_features(r1, r2))

        edge_file = os.path.join(self.work_dir, "{}.edges.gz".format(self.id))
        nx.write_edgelist(structure_graph, edge_file)
        return edge_file

    def get_edge_features(self, r1, r2):
        r1_pos = np.array([self._remove_altloc(a).get_coord() for a in r1]).mean(axis=0)
        r2_pos = np.array([self._remove_altloc(a).get_coord() for a in r2]).mean(axis=0)

        distance = np.linalg.norm(r1_pos-r2_pos)
        angle = angle_between(r1_pos, r2_pos)

        #Add in features from gregarious paper
        #Left handed or right handed: direction of cross product

        if "O" in r1 and "N" in r1 and "O" in r2 and "N" in r2:
            a = np.array(r1["O"].get_coord())-np.array(r1["N"].get_coord())
            b = np.array(r2["O"].get_coord())-np.array(r2["N"].get_coord())

            omega = angle_between(a, b)

            c_a, c_b = np.array(r1["O"].get_coord()), np.array(r2["O"].get_coord())
            mid_a = np.array([r1["O"].get_coord(), r1["N"].get_coord()]).mean(axis=0)
            mid_b = np.array([r2["O"].get_coord(), r2["N"].get_coord()]).mean(axis=0)

            theta = get_dihedral(c_a, mid_a, mid_b, c_b)

            c = mid_b-mid_a

            scalar_triple_product = c*np.cross(a,b)
        else:
            omega = np.nan
            theta = np.nan

        #chirality = int(theta > 0)

        return {
            "distance":distance,
            "angle":angle,
            "omega":omega,
            "theta":theta,
            #"chirality":chirality
        }
