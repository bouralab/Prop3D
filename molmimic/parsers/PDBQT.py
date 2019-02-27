import os, sys
import subprocess
import shutil

from joblib import Memory

from molmimic.util import silence_stdout, silence_stderr
from molmimic.parsers.psize import Psize

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None

def run_open_babel(in_format, in_file, out_format, out_file, work_dir=None, docker=True, job=None):
    """Run APBS. Calculates correct size using Psize and defualt from Chimera
    """
    if work_dir is None:
        work_dir = os.getcwd()

    if docker and apiDockerCall is not None and job is not None:
        parameters = ["-i{}".format(in_format), os.path.basename(in_file),
            "-o{}".format(out_format), "-O", os.path.basename(out_file)]

        if not os.path.abspath(os.path.dirname(in_file)) == os.path.abspath(work_dir):
            shutil.copy(in_file, os.path.join(work_dir, in_file))

        try:
            apiDockerCall(job,
                          image='edraizen/openbabel:latest',
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise

    else:
        raise RuntimeError("Openbabel needs to be run in docker")

    out_file = os.path.join(work_dir, os.path.basename(out_file))
    assert os.path.isfile(out_file), "Outfile not found: {}".format(os.listdir(work_dir))
    return out_file

def get_autodock_features_from_pdbqt(pdbqt_file):
    autodock_features = {}
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_serial = int(line[6:11])
                autodock_type = line[77:79].strip()
                autodock_features[atom_serial] = autodock_type
    return autodock_features

def run_pdbqt_python(pdb_path):
    import pybel
    mol = next(pybel.readfile("pdb", pdb_path))

    mol.addh()

    pdbqt = mol.write("pdbqt")
    autodock_features = {}
    for atom_index in range(mol.OBMol.NumAtoms()):
        a = mol.OBMol.GetAtom(atom_index + 1)

        if a.IsCarbon() and a.IsAromatic():
            element = 'A'
        elif a.IsOxygen():
            element = 'OA'
        elif a.IsNitrogen() and a.IsHbondAcceptor():
            element = 'NA'
        elif a.IsSulfur() and a.IsHbondAcceptor():
            element ='SA'
        else:
            element = "".join([c for c in a.GetType() if c.isalnum()])

        autodock_features[a.GetIdx()] = (element, a.IsHbondDonor())
    return autodock_features

def get_autodock_features(pdb_path, work_dir=None, job=None):
    pdbqt_file = run_open_babel("pdb", pdb_path, "pdbqt", pdb_path+"qt", work_dir=work_dir, job=job)
    autodock_features = get_autodock_features_from_pdbqt(pdbqt_file)
    try:
        os.remove(pdb_path+"qt")
    except OSError:
        pass
    return autodock_features
