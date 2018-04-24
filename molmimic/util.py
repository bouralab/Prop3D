import os
import sys
from contextlib import contextmanager

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
