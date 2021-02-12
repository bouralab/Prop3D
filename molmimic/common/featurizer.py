import os

import pandas as pd
import numpy as np
from sklearn.gaussian_process.kernels import RBF
from Bio import PDB
import freesasa

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from toil.realtimeLogger import RealtimeLogger

from molmimic.util.pdb import InvalidPDB
from molmimic.util import natural_keys, silence_stdout, silence_stderr
from molmimic.util.iostore import IOStore
from molmimic.parsers import mgltools
from molmimic.parsers.FreeSASA import run_freesasa_biopython
from molmimic.parsers.Electrostatics import APBS, Pdb2pqr
from molmimic.parsers.cx import CX
from molmimic.parsers.dssp import DSSP
from molmimic.parsers.eppic import EPPICApi, EPPICLocal

from molmimic.common.Structure import Structure, angle_between, get_dihedral
from molmimic.common.ProteinTables import hydrophobicity_scales
from molmimic.common.features import atom_features, residue_features, \
    atom_features_by_category, residue_features_by_category

class ProteinFeaturizer(Structure):
    def __init__(self, path, cath_domain, job, work_dir,
      input_format="pdb", force_feature_calculation=False, update_features=None, features_path=None, **kwds):
        feature_mode = "w+" if force_feature_calculation else "r"
        if features_path is None and update_features is not None:
            features_path = work_dir
        super(ProteinFeaturizer, self).__init__(
            path, cath_domain,
            input_format=input_format,
            feature_mode=feature_mode,
            residue_feature_mode=feature_mode,
            features_path=features_path,
            **kwds)
        self.job = job
        self.work_dir = work_dir
        self.update_features = update_features

    def calculate_flat_features(self, coarse_grained=False, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False):
        if coarse_grained:
            features = [self.calculate_features_for_residue(
                self._remove_inscodes(r), only_aa=only_aa,
                non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for r in self.structure.get_residues()]
            if self.residue_feature_mode == "w+" or self.update_features is not None:
                self.write_features(coarse_grained=True)
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

    def calculate_flat_residue_features(self, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False):
        return self.calculate_flat_features(coarse_grained=True,
            only_aa=only_aa, only_atom=only_atom,
            non_geom_features=non_geom_features,
            use_deepsite_features=use_deepsite_features)

    def write_features(self, coarse_grained=False):
        if coarse_grained:
            self.residue_features = self.residue_features.astype(np.float64)
            self.residue_features = self.residue_features.drop(columns=
                [col for col in self.residue_features if col not in \
                residue_features])
            self.residue_features.to_hdf(self.residue_features_file, "table")
        else:
            self.atom_features = self.atom_features.astype(np.float64)
            self.atom_features = self.atom_features.drop(columns=
                [col for col in self.atom_features if col not in \
                atom_features])
            self.atom_features.to_hdf(self.atom_features_file, "table")

    def get_features_per_atom(self, residue_list):
        """Get features for eah atom, but not organized in grid"""
        features = [self.get_features_for_atom(self._remove_altloc(a)) for r in residue_list for a in r]
        return features

    def get_features_per_residue(self, residue_list):
        features = [self.get_features_for_residue(self._remove_inscodes(r)) for r in residue_list]
        return features

    def calculate_features_for_atom(self, atom, only_aa=False, only_atom=False,
      non_geom_features=False, use_deepsite_features=False, warn_if_buried=False):
        if self.update_features is not None:
            for feat_type, feat_names in atom_features_by_category.items():
                if feat_type in self.update_features:
                    getattr(self, feat_type)(atom)
                else:
                    use_feats = [feat_name in self.update_features for feat_name \
                        in feat_names]
                    if any(use_feats):
                        getattr(self, feat_type)(atom)
            if warn_if_buried:
                if "residue_buried" not in self.atom_features.columns:
                    is_burried = self.get_accessible_surface_area(atom, save=False)
                else:
                    is_buried = self.atom_features.loc[atom.serial_number, "residue_buried"]
                return self.atom_features, bool(is_buried["residue_buried"])
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
            self.get_deepsite_features(atom, calc_charge=False, calc_conservation=False)
            self.get_evolutionary_conservation_score(atom)

            is_buried = self.atom_features.loc[atom.serial_number, "residue_buried"]

        # RealtimeLogger.info("Finished atom {} {}".format(atom, atom.serial_number))
        # RealtimeLogger.info("Feats {}".format(features))
        # RealtimeLogger.info("Feats shape {}".format(features.shape))
        # RealtimeLogger.info("atom_features shape {}".format(self.atom_features.shape))

        #self.atom_features.loc[atom.serial_number, :] = features

        #RealtimeLogger.info("atom_features is {}".format(self.atom_features))

        if warn_if_buried:
            return self.atom_features, bool(is_buried["residue_buried"])
        else:
            return self.atom_features

    def calculate_features_for_residue(self, residue, only_aa=False, non_geom_features=False,
      use_deepsite_features=False, warn_if_buried=False):
        """Calculate FEATUREs"""
        if self.update_features is not None:
            for feat_type, feat_names in residue_features_categories.items():
                if feat_type in self.update_features:
                    getattr(self, feat_type)(atom)
                else:
                    use_feats = [feat_name in self.update_features for feat_name \
                        in feat_names]
                    if any(use_feats):
                        getattr(self, feat_type)(atom)
            if warn_if_buried:
                if "residue_buried" not in self.residue_features.columns:
                    is_buried = self.get_accessible_surface_area(residue, save=False)
                else:
                    is_buried = self.residue_features.loc[residue.get_id(), "residue_buried"]
                return self.residue_features, bool(is_buried["residue_buried"])
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
            self.get_vdw(residue)
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
            prep = mgltools.PrepareReceptor(job=self.job, work_dir=self.work_dir)
            self._autodock = prep.get_autodock_atom_types(self.path)

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
            atom_features_by_category["get_atom_type"]]

    def get_element_type(self, atom):
        elems = "CNOS"
        elem_col = "{}_elem".format(atom.element)

        if elem_col in atom_features_by_category["get_element_type"]:
            self.atom_features.loc[atom.serial_number, elem_col] = 1.
        else:
            self.atom_features.loc[atom.serial_number, "Unk_elem"] = 1.

        return self.atom_features.loc[atom.serial_number,
            atom_features_by_category["get_element_type"]]

    def get_vdw(self, atom_or_residue):
        vdw = super().get_vdw(atom_or_residue)

        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
            self.atom_features.loc[atom.serial_number, "vdw_radii"] = vdw
            return self.atom_features.loc[atom.serial_number, "vdw_radii"]
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            #For coarse graining, not really a vdw radius
            residue = atom_or_residue
            idx = residue.get_id()
            self.residue_features.loc[idx, "vdw_radii"] = vdw
            return self.residue_features.loc[idx, ["vdw_radii"]]

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
                float(charge_value > 0)]
            cols = residue_features_by_category["get_charge_and_electrostatics"][:3]
            if not only_charge:
                charge += [
                    electrostatic_pot_value,
                    float(electrostatic_pot_value < 0)]
                cols += residue_features_by_category["get_charge_and_electrostatics"][3:]

            if only_bool:
                charge = charge[1:3]
                cols = cols[1:3]
                if not only_charge:
                    charge += charge[4:]
                    cols += charge[4:]

            idx = residue.get_id()
            self.residue_features.loc[idx, cols] = charge
            return self.residue_features.loc[idx, cols]
        else:
            raise RuntimeError("Input must be Atom or Residue: {}".format(type(atom_or_residue)))

        if not hasattr(self, "_pqr") or (not only_charge and len(list(self._pqr.values())[0])==1):
            try:
                if only_charge:
                    pdb2pqr = Pdb2Pqr(work_dir=self.work_dir, job=self.job)
                    self._pqr = pdb2pqr.get_charge_from_pdb_file(self.path, with_charge=False)
                else:
                    apbs = APBS(work_dir=self.work_dir, job=self.job)
                    self._pqr = apbs.get_atom_potentials_from_pdb(self.path)
            except (SystemExit, KeyboardInterrupt):
                raise
            except Exception as e:
                raise
                self._pqr = {}
                RealtimeLogger.info("ELECTROSTATICS failed ({}): {}".format(type(e), e))

        atom_id = atom.get_full_id()[3:5]

        if atom_id[1][1] != " ":
            #pdb2pqr removes alternate conformations and only uses the first
            atom_id = (atom_id[0], (atom_id[1][0], " "))

        if only_charge:
            charge_value = self._pqr.get(atom_id, np.nan)
            electrostatic_pot_value = np.nan
        else:
            try:
                charge_value, electrostatic_pot_value = self._pqr[atom_id]
            except KeyError:
                charge_value, electrostatic_pot_value = np.NaN, np.NaN

        charge = [
            charge_value,
            float(charge_value < 0),
            float(charge_value > 0)]
        cols = atom_features_by_category["get_charge_and_electrostatics"][:3]
        if not only_charge:
            charge += [
                electrostatic_pot_value,
                float(electrostatic_pot_value < 0)]
            cols += atom_features_by_category["get_charge_and_electrostatics"][3:]

        if only_bool:
            charge = charge[1:3]
            cols = cols[1:3]
            if not only_charge:
                charge += charge[4:]
                cols += charge[4:]

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

            concavity = np.array([
                concavity_value,
                float(concavity_value <= 2)])

            idx = residue.get_id()
            cols = residue_features_by_category["get_concavity"]
            self.residue_features.loc[idx, cols] = concavity

            return self.residue_features.loc[idx, cols]
        else:
            raise RuntimeErorr("Input must be Atom or Residue")

        if not hasattr(self, "_cx"):
            cx = CX(work_dir=self.work_dir, job=self.job)
            self._cx = cx.get_concavity(self.path)

        concavity_value = self._cx.get(atom.serial_number, np.NaN)


        concavity = np.array([
            concavity_value,
            float(concavity_value <= 2) #is concave
            ])

        idx = atom.serial_number
        cols = atom_features_by_category["get_concavity"]
        self.atom_features.loc[idx, cols] = concavity

        return self.atom_features.loc[idx, cols]

    def get_hydrophobicity(self, atom_or_residue):
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
            hydrophobicity = hydrophobicity_scales["kd"].get(resname, np.nan)
            biological = hydrophobicity_scales["biological"].get(resname, np.nan)
            octanal = hydrophobicity_scales["octanal"].get(resname, np.nan)
        except KeyError:
            hydrophobicity, biological, octanal = np.nan, np.nan, np.nan


        result = np.array([
            hydrophobicity,
            float(hydrophobicity > 0),
            biological,
            octanal
            ])

        if use_atom:
            idx = atom.serial_number
            cols = atom_features_by_category["get_hydrophobicity"]
            self.atom_features.loc[idx, cols] = hydrophobicity
            return self.atom_features.loc[idx, cols]
        else:
            idx = residue.get_id()
            cols = residue_features_by_category["get_hydrophobicity"]
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

    def get_accessible_surface_area_atom(self, atom, save=True):
        """Returns the ASA value from freesasa (if inout is Atom) and the DSSP
        value (if input is Atom or Residue)

        Returns
        -------
        If input is residue a 3-vector is returned, otherwise a 4-vector is returned
        """
        if not isinstance(atom, PDB.Atom.Atom):
            raise RuntimeErorr("Input must be Atom")

        if not hasattr(self, "_sasa"):
            self._sasa = run_freesasa_biopython(self.path)

        sasa, sasa_struct = self._sasa

        try:
            selection = "sele, chain {} and resi {} and name {}".format(
                self.chain, atom.get_parent().get_id()[1], atom.get_id()[0])
            with silence_stdout(), silence_stderr():
                selections = freesasa.selectArea([selection], sasa_struct, sasa)
                atom_area = selections["sele"]
        except (KeyError, AssertionError, AttributeError, TypeError):
            raise
            atom_area = np.NaN

        if save:
            idx = atom.serial_number
            self.atom_features.loc[idx, "atom_asa"] = atom_area
            return self.atom_features.loc[idx, ["atom_asa"]]
        else:
            return pd.Series([atom_area], index=["atom_asa"])

    def get_accessible_surface_area_residue(self, atom_or_residue, acc_threshold=0.2, save=True):
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
            raise RuntimeError("Input must be Atom or Residue, not {}".format(atom_or_residue))

        if not hasattr(self, "_dssp"):
            dssp = DSSP(work_dir=self.work_dir, job=self.job)
            self._dssp = dssp.get_dssp(self.structure, self.path)

        residue_key = [residue.parent.get_id(), residue.get_id()] #[self.chain, residue.get_id()]
        try:
            residue_rasa = float(self._dssp[tuple(residue_key)][3])
        except KeyError as e1:
            residue_key[1] = tuple([" "]+list(residue_key[1])[1:])
            try:
                residue_rasa = float(self._dssp[tuple(residue_key)][3])
            except KeyError as e2:
                residue_rasa = np.nan

        asa = np.array([
            residue_rasa,
            float(residue_rasa < acc_threshold)])

        if is_atom:
            idx = atom.serial_number
            cols = atom_features_by_category["get_accessible_surface_area"][1:]
            self.atom_features.loc[idx, cols] = asa
            if save:
                return self.atom_features.loc[idx, cols]
            else:
                return pd.Series(asa, index=cols)
        else:
            idx = residue.get_id()
            cols = residue_features_by_category["get_accessible_surface_area"]
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
            self.atom_features.loc[atom.serial_number, PDB.Polypeptide.aa3] = 0.
            self.atom_features.loc[atom.serial_number, col] = 1.
            return self.atom_features.loc[atom.serial_number,
                atom_features_by_category["get_residue"]]
        else:
            self.residue_features.loc[residue.get_id(), PDB.Polypeptide.aa3] = 0.
            self.residue_features.loc[residue.get_id(), col] = 1.
            return self.residue_features.loc[residue.get_id(),
                residue_features_by_category["get_residue"]]

    def get_ss(self, atom_or_residue):
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
            dssp = DSSP(work_dir=self.work_dir, job=self.job)
            self._dssp = dssp.get_dssp(self.structure, self.path)

        try:
            atom_ss = self._dssp[residue.get_full_id()[2:]][2]
        except (KeyError, AssertionError, AttributeError, TypeError):
            try:
                #Remove HETATMs
                atom_ss = self._dssp[(residue.get_full_id()[2], (' ', residue.get_full_id()[3][1], ' '))][2]
            except (KeyError, AssertionError, AttributeError, TypeError):
                atom_ss = "X"

        phi, psi = self.get_dihedral_angles(residue)
        if phi is None:
            phi = np.nan
        if psi is None:
            psi = np.nan

        ss = np.array([
            phi,
            np.sin(phi),
            np.cos(phi),
            psi,
            np.sin(psi),
            np.cos(psi),
            float(atom_ss in "GHIT"),
            float(atom_ss in "BE"),
            float(atom_ss not in "GHITBE")])

        if is_atom:
            idx = atom.serial_number
            cols = atom_features_by_category["get_ss"]
            self.atom_features.loc[idx, cols] = ss
            return self.atom_features.loc[idx, cols]
        else:
            idx = residue.get_id()
            cols = residue_features_by_category["get_ss"]
            self.residue_features.loc[idx, cols] = ss
            return self.residue_features.loc[idx, cols]

    def get_deepsite_features(self, atom, calc_charge=True, calc_conservation=True):
        """Use DeepSites rules for autodock atom types
        """
        if not hasattr(self, "_autodock"):
            prep = mgltools.PrepareReceptor(job=self.job, work_dir=self.work_dir)
            self._autodock = prep.get_autodock_atom_types(self.path, verify=True)

        try:
            atom_type, h_bond_donor = self._autodock[atom.serial_number]
        except KeyError:
            atom_type = "  "

        idx = atom.serial_number
        cols = atom_features_by_category["get_deepsite_features"]

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

    def get_evolutionary_conservation_score(self, atom_or_residue, eppic=True,
      only_bool=False, run_eppic_for_domain_on_failure=False):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            use_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
            cols = atom_features_by_category["get_evolutionary_conservation_score"]
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            use_atom = False
            residue = atom_or_residue
            cols = residue_features_by_category["get_evolutionary_conservation_score"]
        else:
            raise RuntimeError("Input must be Atom or Residue")

        idx = residue.get_id()

        if self.update_features is not None and \
          ("get_evolutionary_conservation_score" not in self.update_features or \
          "is_conserved" not in self.update_features or \
          "eppic_entropy" not in self.update_features):
            #No need to update
            return self.atom_features.loc[atom.serial_number, cols] if is_atom else \
                self.residue_features.loc[idx, cols]

        if not hasattr(self, "_eppic"):
            pdbe_store = IOStore.get("aws:us-east-1:molmimic-pdbe-service")
            eppic_store = IOStore.get("aws:us-east-1:molmimic-eppic-service")

            try:
                eppic_api = EPPICApi(self.pdb, eppic_store, pdbe_store,
                    use_representative_chains=False, work_dir=self.work_dir)
                self._eppic = eppic_api.get_entropy_scores(self.chain)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                if run_eppic_for_domain_on_failure:
                    eppic_local = EPPICLocal(work_dir=self.work_dir, job=self.job)
                    self._eppic = eppic_local.get_entropy_scores(self.path)
                else:
                    self._eppic = {}

        result = pd.Series(np.empty(len(cols)), index=cols, dtype=np.float64)
        result["eppic_entropy"] = self._eppic.get(residue.get_id(), np.nan)
        result["is_conserved"] = float(result["eppic_entropy"]<.5)

        if only_bool:
            return result["is_conserved"]

        if use_atom:
            idx = atom.serial_number
            self.atom_features.loc[idx, result.index] = result
        else:
            idx = residue.get_id()
            self.residue_features.loc[idx, result.index] = result

        return result

    def calculate_graph(self, d_cutoff=100.):
        import networkx as nx
        structure_graph = nx.Graph()
        for r1, r2 in self.calculate_neighbors(d_cutoff=d_cutoff):
            structure_graph.add_edge(r1.get_id(), r2.get_id(),
                attr_dict=self.get_edge_features(r1, r2))

        edge_file = os.path.join(self.work_dir, "{}.edges.gz".format(self.id))
        nx.write_edgelist(structure_graph, edge_file)
        return structure_graph, edge_file

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
