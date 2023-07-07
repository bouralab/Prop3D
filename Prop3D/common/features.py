import os
from typing import Any, TypeVar
import yaml
from collections import OrderedDict

from pathlib import Path
import pandas as pd
import numpy as np

"""Functions to load in defualt features to pandas DataFrames
"""

non_geom_features_names = ["get_atom_type", "get_charge_and_electrostatics",
        "get_charge_and_electrostatics", "get_hydrophobicity", "get_residue",
        "get_deepsite_features", "get_evolutionary_conservation_score"]

FeaturesType = TypeVar('FeaturesType', bound='Features')
class Features(object):
    """Read in feature YAML file for eas parsing

    Parameters
    ----------
    custom : bool
        Read in custom features file
    """
    def __init__(self, custom: bool = False) -> None:
        base = Path(__file__).parent
        if not custom:
            features_file = base / "features.yaml"
        else:
            features_file = base.parent / "custom_featurizers" / "custom_features.yaml"

        with features_file.open() as fh:
            self.features = yaml.safe_load(fh)

        self.atom_features = [feature["name"] for category in self.features for feature in \
            list(category.values())[0]]

        self.residue_features = [feature["name"] for category in self.features for feature in \
            list(category.values())[0] if feature["residue"]]

        self.default_atom_features = pd.Series(
            {feature["name"]:feature["default"] for category in self.features for \
                feature in list(category.values())[0]},
            index=self.atom_features,
            dtype=np.float64)

        self.default_residue_features = pd.Series(
            {feature["name"]:feature["default"] for category in self.features for \
                feature in list(category.values())[0] if feature["residue"]},
            index=self.residue_features,
            dtype=np.float64)

        self.atom_feature_thresholds = {feature["name"]:(feature["threshold"], feature["equality"]) \
            for category in self.features for feature in list(category.values())[0]\
            if "threshold" in feature}

        self.residue_feature_thresholds = {feature["name"]:(feature["threshold"], feature["equality"]) \
            for category in self.features for feature in list(category.values())[0] \
            if feature["residue"] and "threshold" in feature}

        self.atom_bool_features = [feature["name"] for category in self.features for feature in \
            list(category.values())[0] if feature["bool"]]

        self.residue_bool_features = [feature["name"] for category in self.features for feature in \
            list(category.values())[0] if feature["residue"] and feature["bool"]]

        self.atom_feature_categories = OrderedDict([(list(category.keys())[0],
            list(category.values())[0]) for category in self.features])

        self.atom_features_by_category = OrderedDict([(list(category.keys())[0],
            [f["name"] for f in list(category.values())[0]]) \
            for category in self.features])

        self.residue_feature_categories = OrderedDict([(list(category.keys())[0],
            [f for f in list(category.values())[0] if f["residue"]]) \
            for category in self.features])

        self.residue_features_by_category = [(list(category.keys())[0],
            [f["name"] for f in list(category.values())[0] if f["residue"]]) \
            for category in self.features]
        self.residue_features_by_category = OrderedDict([rfeat for rfeat in \
            self.residue_features_by_category if len(rfeat[1])>0])

        self.atom_feature_aggregegation = OrderedDict([(feature["name"],feature["aggregate"]) \
            for category in self.features for feature in list(category.values())[0]])

    def default_atom_feature_df(self, natoms: int) -> pd.DataFrame:
        """Get an atom feature DataFrame with default values for every feature. Columns are features names
        
        Parameters
        ----------
        natoms: int
            Number of atoms to include in new data frame
        """
        return pd.concat([self.default_atom_features] * natoms, axis=1).T

    def default_residue_feature_df(self, nres: int) -> pd.DataFrame:
        """Get a residue feature DataFrame with default values for every feature. Columns are features names
        
        Parameters
        ----------
        nres: int
            Number of residues to include in new data frame
        """
        return pd.concat([self.default_residue_features] * nres, axis=1).T

    def default_atom_feature_np(self, natoms: int) -> np.recarray:
        """Get an atom feature record array with default values for every feature. Columns are features names
        
        Parameters
        ----------
        natoms: int
            Number of atoms to include in new data frame
        """
        df = self.default_atom_feature_df(natoms)
        ds_dt = np.dtype({'names':df.columns,'formats':df.dtypes.values()})
        return np.rec.fromarrays(df.values , dtype=ds_dt)

    def default_residue_feature_np(self, nres: int) -> np.recarray:
        """Get a residue feature record array with default values for every feature. Columns are features names
        
        Parameters
        ----------
        nres: int
            Number of residues to include in new data frame
        """
        df = self.default_residue_feature_df(nres)
        ds_dt = np.dtype({'names':df.columns,'formats':df.dtypes.values()})
        return np.rec.fromarrays(df.values , dtype=ds_dt)

    def check_threshold(self, feature_name: str, raw_value: float, residue: bool = False) -> float:
        """Create new boolean features by thresholding continuous value features

        Parameters
        ----------
        feature_name : str
            Name of new feature
        raw_value : flaot
            Value of continuous value feature
        residue : bool
            Feature calculated at the residue level
        """
        if raw_value == np.nan:
            return np.nan

        if not residue:
            threshold, equality = self.atom_feature_thresholds[feature_name]
        else:
            threshold, equality = self.residue_feature_thresholds[feature_name]

        if isinstance(threshold, (list, tuple)):
            if not isinstance(equality, (list, tuple)) or len(threshold)!=2:
                raise RuntimeError("If using two inequality values, they must both be lists or tuples of length 2")
            
        else:
            threshold, equality = [threshold], [equality]

        val = True

        for e, t in zip(equality, threshold):
            if e == ">":
                val_ = raw_value>t
            elif e == "<":
                val_ = raw_value<t
            elif e == "<=":
                val_ = raw_value<=t
            elif e == ">=":
                val_ = raw_value>=t
            elif e == "!=":
                val_ = raw_value!=t
            else:
                raise RuntimeError("Unknown equality")

            val = val and val_

        return float(val)

    def number_of_features(self, only_aa: bool = False, only_atom: bool = False, non_geom_features: bool = False,
      use_deepsite_features: bool = False, coarse_grained: bool = False) -> int:
        """Get the number of features used
        """
        if self.coarse_grained:
            if non_geom_features:
                return sum(len(self.residue_feature_categories[n]) for n in non_geom_features_names \
                    if n in self.residue_feature_names)
            elif only_aa:
                return len(self.residue_feature_categories["get_residue"])
            else:
                return len(self.residue_features)
        else:
            if non_geom_features:
                return sum(len(self.atom_feature_categories[n]) for n in non_geom_features_names)

            elif only_atom:
                return len(self.atom_feature_categories["get_atom_type"])

            elif only_aa and use_deepsite_features:
                return len(self.atom_feature_categories["get_residue"]) + \
                    len(self.atom_feature_categories["get_deepsite_features"])
            elif only_aa:
                return len(self.atom_feature_categories["get_residue"])
            elif use_deepsite_features:
                return len(self.atom_feature_categories["get_deepsite_features"])
            else:
                return len(self.atom_features)

    def get_feature_names(self,  only_aa: bool = False, only_atom: bool = False, use_deepsite_features: bool = False, coarse_grained: bool = False) -> Any:
        if coarse_grained:
            return self.residue_features
        if only_atom:
            return self.atom_features_by_category["get_atom_type"]
        elif only_aa and use_deepsite_features:
            return self.atom_features_by_category["get_residue"]+self.atom_features_by_category["get_deepsite_features"]
        if only_aa:
            return self.atom_features_by_category["get_residue"]
        if use_deepsite_features:
            return self.atom_features_by_category["get_deepsite_features"]

        return self.atom_features
    
class AllFeatures(object):
    """Combine all sets of features together
    
    Parameters
    ----------
    *features : list of Features objects
    """
    def __init__(self, *features: FeaturesType) -> None:
        self.features = features

    def __getattr__(self, __name: str) -> Any:
        """Combine features from many features"""
        result = None
        for feats in self.features:
            if hasattr(feats, __name):
                attr = getattr(feats, __name)
                if result is None:
                    result = attr
                elif isinstance(attr, type(result)):
                    if isinstance(attr, (list, int)):
                        result += attr
                    elif isinstance(attr, dict):
                        result.update(attr)
                    elif isinstance(attr, (pd.DataFrame, pd.Series)):
                        result = pd.concat((result, attr))
                    else:
                        raise AttributeError(f"{__name} has an invalid type")
                else:
                    raise AttributeError(f"{__name} has an invalid type")
            else:
                break
                raise KeyError(f"{__name} not in features")
        if result is None:
            return super().__getattr__(__name)
            raise RuntimeError("Unable to combine features")
        return result

    
default_features = Features()
custom_features = Features(custom=True)
all_features = AllFeatures(default_features, custom_features)
