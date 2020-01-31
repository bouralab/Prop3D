import os
import yaml
from collections import OrderedDict
import pandas as pd

with open(os.path.join(os.path.dirname(__file__), "features.yaml")) as fh:
    features = yaml.safe_load(fh)

atom_features = [feature["name"] for category in features for feature in \
    list(category.values())[0]]

residue_features = [feature["name"] for category in features for feature in \
    list(category.values())[0] if feature["residue"]]

default_atom_features = pd.Series(
    {feature["name"]:feature["default"] for category in features for \
        feature in list(category.values())[0]},
    index=atom_features,
    dtype=pd.np.float64)

default_residue_features = pd.Series(
    {feature["name"]:feature["default"] for category in features for \
        feature in list(category.values())[0] if feature["residue"]},
    index=residue_features,
    dtype=pd.np.float64)

atom_bool_features = [feature["name"] for category in features for feature in \
    list(category.values())[0] if feature["bool"]]

residue_bool_features = [feature["name"] for category in features for feature in \
    list(category.values())[0] if feature["residue"] and feature["bool"]]

atom_feature_categories = OrderedDict([(list(category.keys())[0],
    list(category.values())[0]) for category in features])

atom_features_by_category = OrderedDict([(list(category.keys())[0],
    [f["name"] for f in list(category.values())[0]]) \
    for category in features])

residue_feature_categories = OrderedDict([(list(category.keys())[0],
    [f for f in list(category.values())[0] if f["residue"]]) \
    for category in features])

residue_features_by_category = OrderedDict([(list(category.keys())[0],
    [f["name"] for f in list(category.values())[0] if f["residue"]]) \
    for category in features])

def default_atom_feature_df(natoms):
    return pd.concat([default_atom_features] * natoms, axis=1).T

def default_residue_feature_df(nres):
    return pd.concat([default_residue_features] * nres, axis=1).T

non_geom_features_names = ["get_atom_type", "get_charge_and_electrostatics",
    "get_charge_and_electrostatics", "get_hydrophobicity", "get_residue",
    "get_deepsite_features", "get_evolutionary_conservation_score"]

def number_of_features(only_aa=False, only_atom=False, non_geom_features=False,
  use_deepsite_features=False, coarse_grained=False):
    if coarse_grained:
        if non_geom_features:
            return sum(len(residue_feature_categories[n]) for n in non_geom_features_names \
                if n in residue_feature_names)
        elif only_aa:
            return len(residue_feature_categories["get_residue"])
        else:
            return len(residue_features)
    else:
        if non_geom_features:
            return sum(len(atom_feature_categories[n]) for n in non_geom_features_names)

        elif only_atom:
            return len(atom_feature_categories["get_atom_type"])

        elif only_aa and use_deepsite_features:
            return len(atom_feature_categories["get_residue"]) + \
                len(atom_feature_categories["get_deepsite_features"])
        elif only_aa:
            return len(atom_feature_categories["get_residue"])
        elif use_deepsite_features:
            return len(atom_feature_categories["get_deepsite_features"])
        else:
            return len(atom_features)

def get_feature_names(only_aa=False, only_atom=False, use_deepsite_features=False, coarse_grained=False):
    if coarse_grained:
        return residue_features
    if only_atom:
        return atom_features_by_category["get_atom_type"]
    elif only_aa and use_deepsite_features:
        return atom_features_by_category["get_residue"]+atom_features_by_category["get_deepsite_features"]
    if only_aa:
        return atom_features_by_category["get_residue"]
    if use_deepsite_features:
        return atom_features_by_category["get_deepsite_features"]

    return atom_features

# atom_feature_names = OrderedDict([
#     ("get_atom_type", [
#         "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F",
#         "MG", "P", "SA", "S", "CL", "CA", "MN", "FE", "ZN", "BR", "I", "Unk_atom"]),
#     ("get_element_type", [
#         "C_elem", "N_elem", "O_elem", "S_elem", "H_elem", "F_elem",
#         "MG_elem", "P_elem", "CL_elem", "CA_elem", "MN_elem",
#         "FE_elem", "ZN_elem", "BR_elem", "I_elem", "Unk_elem"]),
#     ("get_vdw", [
#         "vdw_volume"]),
#     ("get_charge_and_electrostatics", [
#         "charge", "neg_charge", "pos_charge", "electrostatic_potential",
#         "is_electronegative"]),
#     ("get_concavity", [
#         "cx", "is_concave"]),
#     ("get_hydrophobicity", [
#         "hydrophobicity", "is_hydrophobic"]),
#     ("get_accessible_surface_area", [
#         "atom_asa", "residue_rasa", "residue_buried"]),
#     ("get_residue",
#         aa3+["Unk_residue"]),
#     ("get_ss",
#         ["is_helix", "is_sheet", "Unk_SS"]),
#     ("get_deepsite_features", [
#         "hydrophobic_atom",
#         "aromatic_atom",
#         "hbond_acceptor",
#         "hbond_donor",
#         "metal"]),
#     ("get_evolutionary_conservation_score", [
#         "eppic_entropy",
#         "is_conserved"])
# ])
#
# atom_feature_names_flat = [col_name for _, col_names in \
#     atom_feature_names.items() for col_name in col_names]
#
# residue_feature_names = OrderedDict([
#     ("get_charge_and_electrostatics", [
#         "charge", "neg_charge", "pos_charge", "electrostatic_potential",
#         "is_electronegative"]),
#     ("get_concavity", [
#         "cx", "is_concave"]),
#     ("get_hydrophobicity", [
#         "hydrophobicity","is_hydrophobic"]),
#     ("get_accessible_surface_area", [
#         "residue_rasa", "residue_buried"]),
#     ("get_residue",
#         aa3+["Unk_residue"]),
#     ("get_ss",
#         ["is_helix", "is_sheet", "Unk_SS"]),
#     ("get_evolutionary_conservation_score", [
#         "eppic_conservation", "is_conserved"])
# ])
#
# residue_feature_names_flat = [col_name for _, col_names in \
#     residue_feature_names.items() for col_name in col_names]
