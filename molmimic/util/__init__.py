import os
import sys
import re
from contextlib import contextmanager

import pandas as pd

data_path_prefix = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

structures_path_prefix = os.path.join(data_path_prefix, "structures")

features_path_prefix = os.path.join(data_path_prefix, "features")

interfaces_path_prefix = os.path.join(data_path_prefix, "interfaces")

def get_structures_path(dataset_name):
    return os.path.join(structures_path_prefix, dataset_name)

def get_features_path(dataset_name):
    features_path = os.path.join(features_path_prefix, dataset_name)
    if not os.path.isdir(features_path):
        features_path = os.path.join(features_path_prefix, "features_{}".format(dataset_name))
    return features_path

def get_interfaces_path(dataset_name):
    return os.path.join(interfaces_path_prefix, dataset_name)

def initialize_data(dataset_name):
    os.environ["MOLMIMIC_STRUCTURES"] = get_structures_path(dataset_name)
    os.environ["MOLMIMIC_FEATURES"] = get_features_path(dataset_name)
    os.environ["MOLMIMIC_INTERFACES"] = get_interfaces_path(dataset_name)

def iter_cdd(use_label=True, use_id=False, label=None, id=None):
    if use_label and not use_id:
        col = 1
    elif not use_label and use_id:
        col = 2

    CDD = pd.read_hdf(os.path.join(data_path_prefix, "MMDB.h5"), "Superfamilies")
    CDD = CDD[["label", "sfam_id"]].drop_duplicates().dropna()

    if label is not None:
        CDD = CDD[CDD["label"]==label]
    elif id is not None:
        CDD = CDD[CDD["sfam_id"]==id]

    CDD["label"] = CDD["label"].apply(lambda cdd: cdd.replace("/", "").replace("'", "\'") if isinstance(cdd, str) else cdd)
    CDD.sort_values("label", inplace=True)

    if use_label and use_id:
        for cdd in CDD.itertuples():
            yield cdd[1], cdd[2]
    else:
        for cdd in CDD.itertuples():
            yield cdd[col]

class InvalidPDB(RuntimeError):
    pass

def download_pdb(id):
    pdbl = PDB.PDBList()
    try:
        fname = pdbl.retrieve_pdb_file(id.upper(), file_format="mmCif")
        if not os.path.isfile(fname):
            raise InvalidPDB(id)
        return fname, "mmcif"
    except IOError:
        raise InvalidPDB(id)

def get_first_chain(pdb_file):
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    return line[21]
    except IOError:
        pass

    return None

def get_all_chains(pdb_file):
    chains = set()
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chains.add(line[21])
    except IOError:
        pass

    return chains

def atof(text, use_int=False):
    converter = int if use_int else float
    try:
        retval = converter(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text, use_int=False):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c, use_int=use_int) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', str(text)) ]

def to_int(s):
    if isinstance(s, int):
        return s
    elif isinstance(s, str):
        return int("".join([d for d in s if d.isdigit()]))
    else:
        raise RuntimeError("Must be an an in or string")

def number_of_lines(path):
    numlines = 0
    with open(path) as f:
        numlines = sum([1 for line in f])
    return numlines

@contextmanager
def silence_stdout():
    new_target = open(os.devnull, "w")
    old_target, sys.stdout = sys.stdout, new_target
    try:
        yield new_target
    finally:
        sys.stdout = old_target

@contextmanager
def silence_stderr():
    new_target = open(os.devnull, "w")
    old_target, sys.stderr = sys.stderr, new_target
    try:
        yield new_target
    finally:
        sys.stderr = old_target
