import os
from itertools import groupby
from typing import Union, TypeVar

import pandas as pd
import numpy as np
from sklearn.gaussian_process.kernels import RBF
from Bio import PDB
import freesasa
import networkx as nx
        

import warnings
warnings.simplefilter('ignore', PDB.PDBExceptions.PDBConstructionWarning)

from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from Prop3D.util.pdb import InvalidPDB
from Prop3D.util import natural_keys, silence_stdout, silence_stderr
from Prop3D.generate_data.data_stores import data_stores
from Prop3D.parsers import mgltools
from Prop3D.parsers.FreeSASA import run_freesasa_biopython
from Prop3D.parsers.Electrostatics import APBS, Pdb2pqr
from Prop3D.parsers.cx import CX
from Prop3D.parsers.dssp import DSSP
from Prop3D.parsers.eppic import EPPICApi, EPPICLocal
from Prop3D.parsers.frustratometeR import FrustratometeR

from Prop3D.common.LocalStructure import LocalStructure, AtomType, ResidueType, angle_between, get_dihedral
from Prop3D.common.ProteinTables import hydrophobicity_scales
from Prop3D.common.features import default_features, custom_features, all_features

class ProteinFeaturizer(LocalStructure):
    """An object to calculate biophysical properties from a single protein in a local structure file

    Parameters
    ----------
    path : str
        Path to local structure file
    cath_domain : str
        Name of protein domain
    job : toil.job or None
        Toil Job, needed for some application to run from Toil
    work_dir : str
        Where to save all temporary files from all run software. If None, use cwd.
    input_format : str
        Input format, "pdb", "mmCIF", what ever Bio.PDB understands 
    force_feature_calculation : bool
        Recalculate all features even already present. Defualt is False, do not ovewrite. 
    update_features : list
        List of features names or feature groups to update (while keeping the rest the same). Defualt is None, update all features
    features_path : path or None
        Path to save feature file. If None, use cwd.
    """
    def __init__(self, path: str, cath_domain: str, job: Union[Job, None], work_dir: Union[str, None],
      input_format:str = "pdb", force_feature_calculation: bool = False, update_features: list[str] = None, features_path: str = None, **kwds) -> None:
        feature_mode = "w+" if force_feature_calculation else "r"
        if features_path is None: # and update_features is not None:
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

    def calculate_flat_features(self, coarse_grained: bool = False, only_aa: bool = False, only_atom: bool = False,
      non_geom_features: bool = False, use_deepsite_features: bool = False, write: bool = True) -> tuple[list[pd.DataFrame], str]:
        """Calculate features for each atom (or residue)

        Parameters
        ----------
        coarse_grained : bool
            Use residue features. Default false (Atoms features)
        only_aa : bool
            Only calculate features at the residue level. Defualt false.
        only_atom : bool
            Only calculate features at the atom level. Defualt false.
        non_geom_features : bool
            Only calculate features at that are non-geometric Defualt false.
        use_deepsite_features : bool
            Copy features first used by DeepSite by classifying autodock names
        write : bool
            Save features to file after calculating.

        Returns
        -------
        features : pd.DataFrame
            All caculated features
        feature_file : Str
            Path to where features were written to (if chosen to write)
        """
        if coarse_grained:
            features = [self.calculate_features_for_residue(
                self._remove_inscodes(r), only_aa=only_aa,
                non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for r in self.structure.get_residues()]
            if write and (self.residue_feature_mode == "w+" or self.update_features is not None):
                self.write_features(coarse_grained=True)
            return features, self.residue_features_file
        else:
            features = [self.calculate_features_for_atom(
                self._remove_altloc(atom), only_aa=only_aa,
                only_atom=only_atom, non_geom_features=non_geom_features,
                use_deepsite_features=use_deepsite_features) \
                for atom in self.structure.get_atoms()]
            if write and (self.atom_feature_mode == "w+" or self.update_features is not None):
                self.write_features()
            return features, self.atom_features_file

    def calculate_flat_residue_features(self, only_aa: bool = False, only_atom: bool = False,
      non_geom_features: bool = False, use_deepsite_features: bool = False, write: bool = True) -> tuple[list[pd.DataFrame], str]:
        """See calculate_flat_features"""
        return self.calculate_flat_features(coarse_grained=True,
            only_aa=only_aa, only_atom=only_atom,
            non_geom_features=non_geom_features,
            use_deepsite_features=use_deepsite_features, write=write)

    def get_features_per_atom(self, residue_list: list[ResidueType]) -> list[pd.DataFrame]:
        """Get features for each atom in a list of residues"""
        features = [self.get_features_for_atom(self._remove_altloc(a)) for r in residue_list for a in r]
        return features

    def get_features_per_residue(self, residue_list: list[ResidueType]) -> list[pd.DataFrame]:
        """Get features for each atom in a list of residues"""
        features = [self.get_features_for_residue(self._remove_inscodes(r)) for r in residue_list]
        return features

    def calculate_features_for_atom(self, atom: PDB.Atom, only_aa: bool = False, only_atom: bool = False,
      non_geom_features: bool = False, use_deepsite_features: bool = False, warn_if_buried: bool = False) -> Union[pd.DataFrame, tuple[pd.DataFrame, bool]]:
        """Calculate features for a single atom
        
        Parameters
        ----------
        atom : Bio.PDB.Atom
            Atom that needs features
        only_aa : bool
            Only calculate features at the residue level. Defualt false.
        only_atom : bool
            Only calculate features at the atom level. Defualt false.
        non_geom_features : bool
            Only calculate features at that are non-geometric Defualt false.
        use_deepsite_features : bool
            Copy features first used by DeepSite by classifying autodock names
        warn_if_buried : bool
            returns a bool if residue is buried or not from DSSP

        Returns
        -------
        atom_features : pd.DataFrame
            All calcuated atom features
        is_burried : bool
            Is residue buried or not
        """
        if self.update_features is not None:
            for feat_type, feat_names in all_features.atom_features_by_category.items():
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
            self.get_frustration(atom)
            
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
            self.get_frustration(atom)
            self.get_custom_features(atom)

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

    def calculate_features_for_residue(self, residue: ResidueType, only_aa: bool = False, non_geom_features: bool = False,
                                       use_deepsite_features: bool = False, warn_if_buried: bool = False) -> Union[pd.DataFrame, tuple[pd.DataFrame, bool]]:
        """Calculate features for a single atom
        
        Parameters
        ----------
        residue : Bio.PDB.Residue
            Residue that needs features
        only_aa : bool
            Only calculate features at the residue level. Defualt false.
        non_geom_features : bool
            Only calculate features at that are non-geometric Defualt false.
        use_deepsite_features : bool
            Copy features first used by DeepSite by classifying autodock names
        warn_if_buried : bool
            returns a bool if residue is buried or not from DSSP

        Returns
        -------
        atom_features : pd.DataFrame
            All calcuated atom features
        is_burried : bool
            Is residue buried or not
        """
        if self.update_features is not None:
            for feat_type, feat_names in residue_features_by_category.items():
                if feat_type in self.update_features:
                    getattr(self, feat_type)(residue)
                else:
                    use_feats = [feat_name in self.update_features for feat_name \
                        in feat_names]
                    if any(use_feats):
                        getattr(self, feat_type)(residue)
            if warn_if_buried:
                if "residue_buried" not in self.residue_features.columns:
                    is_buried = self.get_accessible_surface_area(residue, save=False)
                else:
                    is_buried = self.residue_features.loc[residue.get_id(), "residue_buried"]
                return self.residue_features, bool(is_buried["residue_buried"])
            else:
                return self.residue_features

        if non_geom_features:
            self.get_residue(residue)
            self.get_charge_and_electrostatics(residue)
            self.get_hydrophobicity(residue)
            self.get_evolutionary_conservation_score(residue)
            self.get_frustration(residue)
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
            self.get_frustration(residue)
            self.get_custom_features(residue)
            is_buried = self.residue_features.loc[[residue.get_id()], "residue_buried"]

        if warn_if_buried:
            return self.residue_features, bool(is_buried["residue_buried"])
        else:
            return self.residue_features

    def get_atom_type(self, atom: AtomType) -> pd.DataFrame:
        """Get Autodock atom type for Bio.PDB.Atom"""

        if not hasattr(self, "_autodock"):
            prep = mgltools.PrepareReceptor(job=self.job, work_dir=self.work_dir)
            self._autodock = prep.get_autodock_atom_types(self.path)

        try:
            atom_type, h_bond_donor = self._autodock[int(atom.serial_number)]
        except KeyError:
            atom_type = "Unk_atom"

        if atom_type not in default_features.atom_features_by_category["get_atom_type"]: # == "  ":
            atom_type = "Unk_atom"

        values = default_features.default_atom_features[default_features.atom_features_by_category["get_atom_type"]]
        values[atom_type] = 1.0

        self.atom_features.loc[atom.serial_number, default_features.atom_features_by_category["get_atom_type"]] = values #1.0

        return self.atom_features.loc[atom.serial_number,
            default_features.atom_features_by_category["get_atom_type"]]

    def get_element_type(self, atom: AtomType) -> pd.DataFrame:
        """Get element name for Bio.PDB.Atom"""
        elems = "CNOS"
        elem_col = "{}_elem".format(atom.element)

        values = default_features.default_atom_features[default_features.atom_features_by_category["get_element_type"]]

        if elem_col in default_features.atom_features_by_category["get_element_type"]:
            values[elem_col] = 1.0
        else:
            values["Unk_elem"] = 1.0

        self.atom_features.loc[atom.serial_number, default_features.atom_features_by_category["get_element_type"]] = values

        return self.atom_features.loc[atom.serial_number,
            default_features.atom_features_by_category["get_element_type"]]

    def get_vdw(self, atom_or_residue: Union[AtomType, ResidueType]) -> pd.DataFrame:
        """Get Van der Waals radius for Bio.PDB.Atom"""
        vdw = super().get_vdw(atom_or_residue)

        if isinstance(atom_or_residue, AtomType):
            atom = atom_or_residue
            self.atom_features.loc[atom.serial_number, "vdw_radii"] = vdw
            return self.atom_features.loc[atom.serial_number, "vdw_radii"]
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            #For coarse graining, not really a vdw radius
            residue = atom_or_residue
            idx = residue.get_id()
            self.residue_features.loc[idx, "vdw_radii"] = vdw
            return self.residue_features.loc[idx, ["vdw_radii"]]

    def get_charge_and_electrostatics(self, atom_or_residue: Union[AtomType, ResidueType], only_charge: bool = False,
                                      only_bool: bool = False, calculate: bool = True) -> pd.DataFrame:
        """Run pdb2par and APBS on an atom or residue to get charge and electrostatic information

        Parameters
        ----------
        atom_or_residue : Bio.PDB.Atom or Bio.PDB.Residue
        only_charge : bool
            Only run pdb2pqr to get charge
        only_bool : bool
            Only save threshld values (bool). Default is False
        calculate : bool
            Should calculate values, otherwise just look them up. Defautl True
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
            return self.get_charge_and_electrostatics_for_atom(atom,
                only_charge=only_charge, only_bool=only_bool, calculate=calculate)
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            return self.get_charge_and_electrostatics_for_residue(residue,
                only_charge=only_charge, only_bool=only_bool, calculate=calculate)
        else:
            raise RuntimeError("Input must be Atom or Residue: {}".format(type(atom_or_residue)))

    def get_charge_and_electrostatics_for_residue(self, atom_or_residue: Union[AtomType, ResidueType], only_charge: bool = False,
                                                  only_bool: bool = False, calculate: bool = True) -> pd.DataFrame:
        """Run pdb2par and APBS on an atom or residue to get charge and electrostatic information at the residue level

        Parameters
        ----------
        atom_or_residue : Bio.PDB.Atom or Bio.PDB.Residue
        only_charge : bool
            Only run pdb2pqr to get charge
        only_bool : bool
            Only save threshld values (bool). Default is False
        calculate : bool
            Should calculate values, otherwise just look them up. Defautl True
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

        if self.update_features is not None:
            if "electrostatic_potential" not in self.update_features and \
              "is_electropositive" not in self.update_features and \
              "is_electronegative" not in self.update_features:
                 only_charge = True
            if "get_charge_and_electrostatics" not in self.update_features:
                calculate = False

        atoms = pd.concat([self.get_charge_and_electrostatics_for_atom(
            self._remove_altloc(a)) for a in residue], axis=0)
        charge_value = atoms["charge"].sum()
        electrostatic_pot_value = atoms["electrostatic_potential"].sum()

        charge = [
            charge_value,
            default_features.check_threshold("neg_charge", charge_value, residue=True), #float(charge_value < 0)
            default_features.check_threshold("pos_charge", charge_value, residue=True) #float(charge_value > 0)
            ]
        cols = default_features.residue_features_by_category["get_charge_and_electrostatics"][:3]
        if not only_charge:
            charge += [
                electrostatic_pot_value,
                default_features.check_threshold("is_electronegative", electrostatic_pot_value, residue=True) #float(electrostatic_pot_value < 0)
                ]
            cols += default_features.residue_features_by_category["get_charge_and_electrostatics"][3:]

        if only_bool:
            charge = charge[1:3]
            cols = cols[1:3]
            if not only_charge:
                charge += charge[4:]
                cols += charge[4:]

        idx = residue.get_id()
        self.residue_features.loc[idx, cols] = charge
        return self.residue_features.loc[idx, cols]

    def get_charge_and_electrostatics_for_atom(self, atom: AtomType, only_charge: bool = False,
                                               only_bool: bool = False, calculat: bool = True) -> pd.DataFrame:
        """Run pdb2par and APBS on an atom to get carge and electrostatic information

        Parameters
        ----------
        atom: Bio.PDB.Atom
        only_charge : bool
            Only run pdb2pqr to get charge
        only_bool : bool
            Only save threshld values (bool). Default is False
        calculate : bool
            Should calculate values, otherwise just look them up. Defautl True
        """
        if not isinstance(atom, PDB.Atom.Atom):
            raise RuntimeError("Input must be Atom")

        if self.update_features is not None:
            if "electrostatic_potential" not in self.update_features and \
              "is_electropositive" not in self.update_features and \
              "is_electronegative" not in self.update_features:
                 only_charge = True
            if "get_charge_and_electrostatics" not in self.update_features:
                calculate = False

        if not hasattr(self, "_pqr"):
            self._pqr = {}
        if calculate and (len(self._pqr)==0 or (not only_charge and len(list(self._pqr.values())[0])==1)): #not hasattr(self, "_pqr")
            try:
                if only_charge:
                    pdb2pqr = Pdb2pqr(work_dir=self.work_dir, job=self.job)
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

        if calculate:
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
                default_features.check_threshold("neg_charge", charge_value), #float(charge_value < 0)
                default_features.check_threshold("pos_charge", charge_value) #float(charge_value > 0)
                ]

        cols = default_features.atom_features_by_category["get_charge_and_electrostatics"][:3]
        if not only_charge:
            if calculate:
                charge += [
                    electrostatic_pot_value,
                    check_threshold("is_electronegative", electrostatic_pot_value) #float(electrostatic_pot_value < 0)
                ]
            cols += default_features.atom_features_by_category["get_charge_and_electrostatics"][3:]

        if only_bool:
            if calculate:
                charge = charge[1:3]
            cols = cols[1:3]
            if not only_charge:
                if calculate:
                    charge += charge[4:]
                cols += charge[4:]

        idx = atom.serial_number
        if calculate:
            charge = np.array(charge)
            self.atom_features.loc[idx, cols] = charge

        return self.atom_features.loc[idx, cols]

    def get_concavity(self, atom_or_residue: Union[AtomType, ResidueType]) -> pd.DataFrame:
        """Get concavity of an atom or residue by running CX
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            atom = atom_or_residue
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            residue = atom_or_residue
            concavity_value = pd.concat([self.get_concavity(
                self._remove_altloc(a)) for a in residue], axis=0)["cx"].mean()

            concavity = np.array([
                concavity_value,
                default_features.check_threshold("is_concave", concavity_value, residue=True) #float(concavity_value <= 2)
                ])

            idx = residue.get_id()
            cols = default_features.residue_features_by_category["get_concavity"]
            self.residue_features.loc[idx, cols] = concavity

            return self.residue_features.loc[idx, cols]
        else:
            raise RuntimeError("Input must be Atom or Residue")

        if not hasattr(self, "_cx"):
            cx = CX(work_dir=self.work_dir, job=self.job)
            self._cx = cx.get_concavity(self.path)

        concavity_value = self._cx.get(atom.serial_number, np.NaN)


        concavity = np.array([
            concavity_value,
            default_features.check_threshold("is_concave", concavity_value) #float(concavity_value <= 2)
            ])

        idx = atom.serial_number
        cols = default_features.atom_features_by_category["get_concavity"]
        self.atom_features.loc[idx, cols] = concavity

        return self.atom_features.loc[idx, cols]

    def get_hydrophobicity(self, atom_or_residue: Union[AtomType, ResidueType]) -> pd.DataFrame:
        """Get the hydrophocity of an atom or residue using three different scales: 
            Kyte-Doolite (kd), Biological, and Octanal
        """
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
            default_features.check_threshold("is_hydrophobic", hydrophobicity, residue=not use_atom), #float(concavity_value <= 2)
            biological,
            octanal
            ])

        if use_atom:
            idx = atom.serial_number
            cols = default_features.atom_features_by_category["get_hydrophobicity"]
            self.atom_features.loc[idx, cols] = result
            return self.atom_features.loc[idx, cols]
        else:
            idx = residue.get_id()
            cols = default_features.residue_features_by_category["get_hydrophobicity"]
            self.residue_features.loc[idx, cols] = result
            return self.residue_features.loc[idx, cols]

    def get_accessible_surface_area(self, atom_or_residue: Union[AtomType, ResidueType], save: bool = True)  -> Union[pd.DataFrame, pd.Series]:
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

    def get_accessible_surface_area_atom(self, atom: AtomType, save: bool = True) -> Union[pd.DataFrame, pd.Series]:
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

    def get_accessible_surface_area_residue(self, atom_or_residue: Union[AtomType, ResidueType], acc_threshold: float = 0.2, 
                                            save: bool = True) -> Union[pd.DataFrame, pd.Series]:
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
            check_threshold("residue_buried", residue_rasa, residue=True)
        ])

        if is_atom:
            idx = atom.serial_number
            cols = default_features.atom_features_by_category["get_accessible_surface_area"][1:]
            self.atom_features.loc[idx, cols] = asa
            if save:
                return self.atom_features.loc[idx, cols]
            else:
                return pd.Series(asa, index=cols)
        else:
            idx = residue.get_id()
            cols = default_features.residue_features_by_category["get_accessible_surface_area"]
            self.residue_features.loc[idx, cols] = asa
            if save:
                return self.residue_features.loc[idx, cols]
            else:
                return pd.Series(asa, index=cols)

    def get_residue(self, atom_or_residue: Union[AtomType, ResidueType]) -> pd.DataFrame:
        """Get a one hote encoded represatnation the amino acid from atom or residue
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

        # residues = [0.]*(len(PDB.Polypeptide.aa1)+1) #one hot residue
        #
        # try:
        #     residues[PDB.Polypeptide.three_to_index(residue.get_resname())] = 1.
        # except (ValueError, KeyError) as e:
        #     residues[-1] = 1.
        # return np.array(residues)

        col = residue.get_resname() if residue.get_resname() in PDB.Polypeptide.aa3 \
            else "Unk_element"

        all_aas = list(PDB.Polypeptide.aa3)+["Unk_element"]

        if is_atom:
            self.atom_features.loc[atom.serial_number, all_aas] = 0.
            self.atom_features.loc[atom.serial_number, col] = 1.
            return self.atom_features.loc[atom.serial_number,
                default_features.atom_features_by_category["get_residue"]]
        else:
            self.residue_features.loc[residue.get_id(), all_aas] = 0.
            self.residue_features.loc[residue.get_id(), col] = 1.
            return self.residue_features.loc[residue.get_id(),
                default_features.residue_features_by_category["get_residue"]]

    def get_ss(self, atom_or_residue: Union[AtomType, ResidueType]) -> pd.Series:
        """Get 3- and 7- secondary strcutre codes from DSSP for an atom or residue as well as the phi and psi angles
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
            float(atom_ss in "GH"),
            float(atom_ss in "BE"),
            float(atom_ss not in "GHBE"),
            float(atom_ss == "H"),
            float(atom_ss == "B"),
            float(atom_ss == "E"),
            float(atom_ss == "G"),
            float(atom_ss == "I"),
            float(atom_ss == "T"),
            float(atom_ss == "S"),
            float(atom_ss in ["", "-", None, "None"])
        ])

        if is_atom:
            idx = atom.serial_number
            cols = default_features.atom_features_by_category["get_ss"]
            self.atom_features.loc[idx, cols] = ss
            return self.atom_features.loc[idx, cols]
        else:
            idx = residue.get_id()
            cols = default_features.residue_features_by_category["get_ss"]
            self.residue_features.loc[idx, cols] = ss
            return self.residue_features.loc[idx, cols]

    def get_deepsite_features(self, atom: AtomType, calc_charge: bool = True, calc_conservation: bool = True):
        """Use DeepSite rules for autodock atom types: 
            is_hydrophobic (C or A)
            is_aromatic (A)
            is_hbond_acceptor (NA, NS, OA, OS, or SA)
            is_hbond_donor (HS or HD)
            is_metal (MG, ZN, MN, CA, FA) (not relavent since all hetatms are likley removed)

        Parameters
        ----------
        atom : Bio.PDB.Atom
            Atom to apply features to
        calc_charge : bool
            Also calculate charge if not pre-computed. Defualt is True.
        calc_conservation : bool
            Also get conservation values
        """
        if not isinstance(atom, PDB.Atom.Atom):
            raise RuntimeError("Input must be Atom")

        if not hasattr(self, "_autodock"):
            prep = mgltools.PrepareReceptor(job=self.job, work_dir=self.work_dir)
            self._autodock = prep.get_autodock_atom_types(self.path, verify=True)

        try:
            atom_type, h_bond_donor = self._autodock[atom.serial_number]
        except KeyError:
            atom_type = "  "

        idx = atom.serial_number
        cols = default_features.atom_features_by_category["get_deepsite_features"][:5]

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
            charge = self.get_charge_and_electrostatics(atom, only_charge=True, only_bool=True, calculate=False)
            cols += charge.index.tolist()

        if calc_conservation:
            cons = self.get_evolutionary_conservation_score(atom, only_bool=True)
            cols += cons.index.tolist()

        return self.atom_features.loc[idx, cols]

    def get_evolutionary_conservation_score(self, atom_or_residue: Union[AtomType, ResidueType], eppic=True,
                                            only_bool=False, run_eppic_for_domain_on_failure=False) -> pd.Series:
        """Get the evolutionary conservation score measured by the EPPIC entropy value.

        Parameters
        ----------
        eppic : bool
            DEPRECATED. It is always using EPPIC now
        only_bool : bool
            Only save threshold values
        run_eppic_for_domain_on_failure : bool
            If protein does not exist in the EPPIC server, run local. Default False
        """
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            use_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
            cols = default_features.atom_features_by_category["get_evolutionary_conservation_score"]
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            use_atom = False
            residue = atom_or_residue
            cols = default_features.residue_features_by_category["get_evolutionary_conservation_score"]
        else:
            raise RuntimeError("Input must be Atom or Residue")

        idx = residue.get_id()

        if self.update_features is not None and \
          ("get_evolutionary_conservation_score" not in self.update_features or \
          "is_conserved" not in self.update_features or \
          "eppic_entropy" not in self.update_features):
            #No need to update
            return self.atom_features.loc[atom.serial_number, cols] if use_atom else \
                self.residue_features.loc[idx, cols]

        if not hasattr(self, "_eppic"):
            try:
                eppic_api = EPPICApi(self.pdb[:4], data_stores(self.job).eppic_store, data_stores(self.job).pdbe_store,
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
        result["is_conserved"] = default_features.check_threshold("is_conserved", result["eppic_entropy"], residue=not use_atom)

        if only_bool:
            return result["is_conserved"]

        if use_atom:
            idx = atom.serial_number
            self.atom_features.loc[idx, result.index] = result
        else:
            idx = residue.get_id()
            self.residue_features.loc[idx, result.index] = result

        return result

    def get_frustration(self, atom_or_residue: Union[AtomType, ResidueType]) -> None:
        """Calculate frustration for an atom or residue"""
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            use_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
            cols = default_features.atom_features_by_category["get_frustration"]
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            use_atom = False
            residue = atom_or_residue
            cols = default_features.residue_features_by_category["get_frustration"]
        else:
            raise RuntimeError("Input must be Atom or Residue")

        idx = residue.get_id()

        if self.update_features is not None and "get_frustration" not in self.update_features:
            #No need to update
            return self.atom_features.loc[atom.serial_number, cols] if use_atom else \
                self.residue_features.loc[idx, cols]

        if not hasattr(self, "_frustration_singleresidue"):
            frust = FrustratometeR(work_dir=self.work_dir)
            self._frustration_singleresidue = frust.run(pdb_file=self.path, mode="singleresidue")

        try:
            frust_values = self._frustration_singleresidue[
                (self._frustration_singleresidue.Res=="".join(map(str,idx[1:])).strip())&\
                (self._frustration_singleresidue.ChainRes==self.chain)].iloc[0]
        except IndexError:
            frust_values = {}

        result = pd.Series(np.empty(len(cols)), index=cols, dtype=np.float64)
        result["density_res"] = frust_values.get("DensityRes", np.nan)
        result["native_energy"] = frust_values.get("NativeEnergy", np.nan)
        result["decoy_energy"] = frust_values.get("DecoyEnergy", np.nan)
        result["sd_energy"] = frust_values.get("SDEnergy", np.nan)
        result["frustration_index"] = frust_values.get("FrstIndex", np.nan)

        result["is_highly_frustrated"] = default_features.check_threshold("is_highly_frustrated", frust_values.get("FrstIndex", np.nan), residue=not use_atom)
        result["is_minimally_frustrated"] = default_features.check_threshold("is_minimally_frustrated", frust_values.get("FrstIndex", np.nan), residue=not use_atom)
        result["has_nuetral_frustration"] = default_features.check_threshold("has_nuetral_frustration", frust_values.get("FrstIndex", np.nan), residue=not use_atom)
        

    def calculate_graph(self, d_cutoff: float = 100., edgelist: bool = False, write: bool = True) -> tuple[nx.Graph, str]:
        """Build a residue-residue network, linking residues if the distances between is less than given cutoff and 
        give attributes to each edge. Network is saved as a pd.DataFrame or edge list.
        """
        structure_graph = nx.Graph()
        for r1, r2 in self.calculate_neighbors(d_cutoff=d_cutoff):
            structure_graph.add_edge(r1.get_id(), r2.get_id(),
                attr_dict=self.get_edge_features(r1, r2))

        edge_file = os.path.join(self.work_dir, "{}.edges.gz".format(self.id))
        if write:
            nx.write_edgelist(structure_graph, edge_file)

        if edgelist:
            edges = [{"src":u, "dst":v, **dict(d)["attr_dict"]} for u, v, d in \
                structure_graph.edges(data=True)]
            structure_graph = pd.DataFrame(edges)

        return structure_graph, edge_file

    def get_edge_features(self, r1: ResidueType, r2: ResidueType) -> dict[str, float]:
        """Calculate edge features for a pair of residues.

        Features:
        distance,
        angle,
        omega,
        theta, (or dihedral)
        
        #Frustration:
        native_energy
        decoy_energy
        sd_energy
        frustration_index
        is_highly_frustrated
        is_minimally_frustrated
        has_nuetral_frustration
        is_water_mediated_welltype
        is_short_welltype
        is_long_welltype
        """
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

        #Calculate pairwise frustration
        #NativeEnergy DecoyEnergy SDEnergy FrstIndex Welltype FrstState

        r1_idx = r1.get_id()
        r2_idx = r2.get_id()

        if not hasattr(self, "_frustration_configutational"):
            frust = FrustratometeR(work_dir=self.work_dir)
            self._frustration_configutational = frust.run(pdb_file=self.path, mode="configurational")

        try:
            frust_values = self._frustration_configutational[
                (self._frustration_configutational.Res1=="".join(map(str, r1_idx[1:])).strip())&\
                (self._frustration_configutational.ChainRes1==self.chain)&\
                (self._frustration_configutational.Res2=="".join(map(str, r2_idx[1:])).strip())&\
                (self._frustration_configutational.ChainRes2==self.chain)].iloc[0]
        except IndexError:
            frust_values = {}

        return {
            "distance":distance,
            "angle":angle,
            "omega":omega,
            "theta":theta,
            #"chirality":chirality

            #Frustration:
            "native_energy": frust_values.get("NativeEnergy", np.nan),
            "decoy_energy": frust_values.get("DecoyEnergy", np.nan),
            "sd_energy": frust_values.get("SDEnergy", np.nan),
            "frustration_index": frust_values.get("FrstIndex", np.nan),
            "is_highly_frustrated": frust_values.get("FrstState","")=="highly",
            "is_minimally_frustrated": frust_values.get("FrstState","")=="minimally",
            "has_nuetral_frustration": frust_values.get("FrstState","")=="neutral",
            "is_water_mediated_welltype": frust_values.get("Welltype","")=="water-mediated",
            "is_short_welltype": frust_values.get("Welltype","")=="short",
            "is_long_welltype": frust_values.get("Welltype","")=="long",
        }

    def calculate_custom_features(self) -> None:
        if not hasattr(self, "custom_feature_values"):
            self.custom_feature_values = {}
            self.custom_feature_modules = {}
        
        for category in custom_features.features:
            for parser_name, features in groupby(list(category.values())[0], key=lambda feat: feat["parser"]):
                try:
                    parser = self.custom_feature_modules[category]
                except KeyError:
                    parser = import_module(parser_name)(work_dir=self.work_dir)
                    self.custom_feature_modules[category] = parser
                
                atom_results, residue_results = parser.calculate_prop3D(self.path, self.structure)

                if atom_results is None:
                    residue_results = pd.DataFrame(np.nan, index=self.residue_features.index)
                
                if residue_results is None:
                    residue_results = pd.DataFrame(np.nan, index=self.residue_features.index)#, columns=custom_features.atom_features_by_category[category])
                
                if len(residue_results.columns) != custom_features.residue_features_by_category[category] or \
                    len(atom_results.columns) != custom_features.atom_features_by_category[category]:
                    #Map features from atoms to residues and residues to atoms
                    resi_to_atoms = {resi.get_id():[a.serial for a in resi] for resi in self.get_residues()}

                    #Create boolean/thresholded features
                    update_atom_cols = False
                    for feat in custom_features.atom_feature_categories[category]:
                        if feat["name"] not in atom_results:
                            try:
                                source_feature = feat["from_feature"]
                            except KeyError:
                                raise KeyError(f"Feature '{feat['name']} has no source feature ('from_feature') attribute")
                            
                            new_feat = atom_results[source_feature].apply(lambda v: custom_features.check_threshold(feat["name"], v, residue=False))
                            atom_results = atom_results.assign(**{feat["name"]: new_feat})
                            update_atom_cols = True
                    if update_atom_cols:
                        atom_results = atom_results[custom_features.atom_features_by_category(category)]

                    update_resi_cols = False
                    for feat in custom_features.residue_feature_categories[category]:
                        if feat["name"] not in residue_results:
                            try:
                                source_feature = feat["from_feature"]
                            except KeyError:
                                raise KeyError(f"Feature '{feat['name']} has no source feature ('from_feature') attribute")
                            
                            new_feat = residue_results[source_feature].apply(lambda v: custom_features.check_threshold(feat["name"], v, residue=True))
                            residue_results = residue_results.assign(**{feat["name"]: new_feat})
                            update_resi_cols = True
                    if update_resi_cols:
                        residue_results = residue_results[custom_features.residue_features_by_category(category)]

                    #Map atom features to residue feature through aggegation rules
                    aggregate_rules = {"sum": np.sum, "min":np.min, "max":np.max, "mean":np.mean, "avg":np.mean}
                    for feat in custom_features.residue_feature_categories[category]:
                        if feat["name"] not in residue_results and feat.get("from_feature", None) is None:
                            try:
                                agg = aggregate_rules[feat["aggregate"]]
                            except KeyError:
                                raise KeyError(f"Aggregate rule must be: {aggregate_rules.keys()}")
                            
                            residue_results = residue_results.assign(**{feat: np.nan})
                            
                            for resi in self.get_residues():
                                idx = resi_to_atoms[resi.get_id()]
                                residue_results.loc[resi.get_id(), feat] = agg(atom_results.loc[idx, feat])

                    #Map residue features to atoms by copying the features to all atoms in the residue
                    for feat in custom_features.atom_feature_categories[category]:
                        if feat["name"] not in atom_results:
                            atom_results = atom_results.assign(**{feat: np.nan})

                            for resi, atoms in resi_to_atoms.items():
                                atom_results.loc[atoms, feat] = residue_results.loc[resi.get_id(), feat]




