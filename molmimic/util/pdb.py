import os
import sys
import re
import shutil
import subprocess
import glob

import numpy as np
import pandas as pd
from botocore.exceptions import ClientError

from molmimic.util import SubprocessChain, natural_keys
from molmimic.util.iostore import IOStore

#Auto-scaling on AWS with toil has trouble finding modules? Heres the workaround
PDB_TOOLS = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdb_tools")

class InvalidPDB(RuntimeError):
    pass

def PDBTools(commands, output):
    cmds = [[sys.executable, "-m", "pdb-tools.pdb_{}".format(cmd[0])]+cmd[1:] \
        for cmd in commands]
    SubprocessChain(cmd, output)

def download_pdb(id):
    from Bio import PDB
    pdbl = PDB.PDBList()
    try:
        fname = pdbl.retrieve_pdb_file(id.upper(), file_format="mmCif")
        if not os.path.isfile(fname):
            raise InvalidPDB(id)
        return fname, "mmcif"
    except IOError:
        raise InvalidPDB(id)

def get_atom_lines(pdb_file):
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    yield line
    except IOError:
        pass

def get_first_chain(pdb_file):
    for line in get_atom_lines(pdb_file):
        return line[21]
    return None

def get_all_chains(pdb_file):
    chains = set()
    for line in get_atom_lines(pdb_file):
        chains.add(line[21])
    return chains

def is_ca_model(pdb_file):
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM") and line[13:16] != 'CA ':
                    return False
        return True
    except IOError:
        return False

def get_pdb_residues(pdb_file):
    try:
        with open(pdb_file) as f:
            prev_res = None
            for line in f:
                if line.startswith("ATOM"):
                    res = natural_keys(line[22:27], use_int=True) #inlcude icode
                    if not res == prev_res:
                        yield res
                    prev_res = res
    except IOError:
        pass

def read_pdb(file):
    return np.array([(float(line[30:38]), float(line[38:46]), float(line[46:54])) \
        for line in get_atom_lines(file)])

def replace_chains(pdb_file, new_file, **chains):
    """Modified from pdbotools"""
    coord_re = re.compile('^(ATOM|HETATM)')

    with open(pdb_file) as f, open(new_file, "w") as new:
        for line in f:
            if coord_re.match(line) and line[21] in chains:
                print(line[:21] + chains[line[21]] + line[22:], file=new)
            else:
                print(line, file=new)

    print("NEW FILE", new_file)
    assert os.path.isfile(new_file), new_file
    return new_file

def extract_chains(pdb_file, chains, rename=None, new_file=None):
    """Modified from pdbotools"""
    coord_re = re.compile('^(ATOM|HETATM)')

    if new_file is None:
        name, ext = os.path.splitext(pdb_file)
        new_file = "{}.{}.pdb".format(name, chains)

    if isinstance(rename, (str, list, tuple)):
        assert len(chains) == len(rename), "'{}', '{}'".format(chains, rename)
        replace = dict(list(zip(chains, rename)))
        get_line = lambda l: line[:21] + replace[line[21]] + line[22:]
    else:
        get_line = lambda l: l

    with open(pdb_file) as f, open(new_file, "w") as new:
        for line in f:
            if coord_re.match(line) and line[21] in chains:
                new.write(get_line(line))

    return new_file

def update_xyz(old_pdb, new_pdb, updated_pdb=None, process_new_lines=None):
    if updated_pdb is None:
        updated_pdb = "{}.rottrans.pdb".format(os.path.splitext(old_pdb)[0])

    if process_new_lines is None:
        def process_new_lines(f):
            for line in get_atom_lines(f):
                yield line[30:54]

    with open(updated_pdb, "w") as updated:
        for old_line, new_line in zip(
          get_atom_lines(old_pdb),
          process_new_lines(new_pdb)):
            updated.write(old_line[:30]+new_line+old_line[54:])

    return updated_pdb

def rottrans(moving_pdb, matrix_file, new_file=None):
    coords = read_pdb(moving_pdb)
    m = pd.read_table(matrix_file, skiprows=1, nrows=3, delim_whitespace=True, index_col=0)
    M = m.iloc[:, 1:].values
    t = m.iloc[:, 1].values
    coords = np.dot(coords, M)+t
    def get_coords(mat):
        for x, y, z in mat:
            print("LINE:{:8.3f}{:8.3f}{:8.3f}".format(x, y, z))
            yield "{:8.3f}{:8.3f}{:8.3f}".format(x, y, z)
    return update_xyz(moving_pdb, coords, updated_pdb=new_file, process_new_lines=get_coords)

def rottrans_from_matrix(moving_pdb, M, t, new_file=None):
    coords = read_pdb(moving_pdb)
    coords = np.dot(coords, M)+t
    def get_coords(mat):
        for xyz in mat:
            #print("LINE:{:8.3f}{:8.3f}{:8.3f}".format(x, y, z))
            #yield "{:8.3f}{:8.3f}{:8.3f}".format(x, y, z)
            print("".join(("{:<8.3f}".format(i)[:8] for i  in xyz)))
            yield "".join(("{:<8.3f}".format(i)[:8] for i  in xyz))
    return update_xyz(moving_pdb, coords, updated_pdb=new_file, process_new_lines=get_coords)

def tidy(pdb_file, replace=False, new_file=None):
    _new_file = pdb_file+".tidy.pdb" if new_file is None else new_file
    with open(_new_file, "w") as f:
        subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), pdb_file], stdout=f)

    if replace and new_file is None:
        try:
            os.remove(pdb_file)
        except OSEror:
            pass
        shutil.move(_new_file, pdb_file)
        return pdb_file
    else:
        return _new_file

def delocc_pdb(pdb_file, updated_pdb=None):
    if updated_pdb is None:
        updated_pdb = "{}.delocc.pdb".format(os.path.splitext(pdb_file)[0])

    with open(updated_pdb, "w") as f:
        subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py"), pdb_file], stdout=f)

    with open(updated_pdb) as f:
        print("UPDATED PDB", f.read())
    return updated_pdb

def remove_ter_lines(pdb_file, updated_pdb=None):
    if updated_pdb is None:
        updated_pdb = "{}.untidy.pdb".format(os.path.splitext(pdb_file)[0])

    with open(pdb_file) as f, open(updated_pdb, "w") as updated:
        for line in f:
            if not line.startswith("TER"):
                updated.write(line)

    return updated_pdb

def build_atom_unique_id(atom_line):
    """Returns a unique identifying tuple from an ATOM line. From pdb_tools"""

    # unique_id: (name, altloc, resi, insert, chain, segid)
    unique_id = (atom_line[12:16],
                 atom_line[16],
                 int(atom_line[22:26]),
                 atom_line[26], atom_line[21],
                 atom_line[72:76].strip())

    return unique_id

def replace_occ_b(reference, target, occ=1.0, bfactor=20.0, updated_pdb=None):
    assert occ<999999
    assert bfactor<999999

    if updated_pdb is None:
        updated_pdb = "{}.occ_b.pdb".format(os.path.splitext(target)[0])

    reference_lines = {build_atom_unique_id(line):line for line in \
        get_atom_lines(reference)}

    with open(updated_pdb, "w") as fh:
        for line in get_atom_lines(target):
            atom_id = build_atom_unique_id(line)
            if atom_id in reference_lines:
                print(line[:54]+reference_lines[atom_id][54:].rstrip(), file=fh)
            else:
                print(line[:54]+"{:6.2f}".format(occ).ljust(6)[:6]+\
                      "{:6.2f}".format(bfactor).ljust(6)[:6], file=fh)

    return updated_pdb

def normalize_b(b):
    pdb_dir = os.path.dirname(pdb_file)
    bfactor_norm_file = os.path.join(pdb_dir, "bfactor.norm")
    if os.path.isfile(bfactor_norm_file):
        with open(bfactor_norm_file) as f:
            mean, std = map(float, f.read().split())
    else:
        mean, std = get_bfactor_norm_mean(pdb_dir)

    return ((b-mean)/std)

def get_bfactor_norm_mean(pdb_dir, force=False):
    from scipy.stats import median_absolute_deviation
    bfactor_norm_file = os.path.join(pdb_dir, "bfactor.norm")

    if not force and os.path.isfile(bfactor_norm_file):
        return

    bfactors_file = os.path.join(pdb_dir, "bfactors.txt")

    #Combine all B-facotrs from every PDB
    with open(bfactors_file, "w") as fh:
        for pdb_file in glob.glob(os.path.join(pdb_dir, "*.pdb")):
            pdb_bfactors = "\n".join(line[60:66] for line in \
                get_atom_lines(pdb_file))
            print(pdb_bfactors, file=fh)

    bfactors = np.loadtxt(bfactors_file)
    raw_mean, raw_std = bfactors.mean(), bfactors.std()
    mad = median_absolute_deviation(bfactors)

    #From Chung, Wang, Bourne. "Exploiting sequence and structure homologs to identify
    #protein–protein binding sites." Proteins: Structure, Function, Bioinformatics. 2005.
    #https://doi.org/10.1002/prot.20741
    M = 0.6745*((bfactors-mean)/mad)
    bfactors = bfactors[M<=3.5]
    mean, std = bfactors.mean(), bfactors.std()

    with open(bfactor_norm_file, "w") as fh:
        print("{} {}".format(mean, std), file=fh)

    return mean, std

def s3_download_pdb(pdb, work_dir=None, remove=False):
    if work_dir is None:
        work_dir = os.getcwd()

    store = IOStore.get("aws:us-east-1:molmimic-pdb")

    pdb_file_base = os.path.join(pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    pdb_file = os.path.join(work_dir, "pdb{}.ent.gz".format(pdb.lower()))

    format = "pdb"

    for file_format, obs in (("pdb", False), ("mmCif", False), ("pdb", True), ("mmCif", True)):
        _pdb_file_base = pdb_file_base.replace(".ent.", ".mmcif.") if file_format == "mmcif" else pdb_file_base
        _pdb_file = pdb_file.replace(".ent.", ".mmcif.") if file_format == "mmcif" else pdb_file
        _pdb_file_base = "obsolete/"+_pdb_file_base if obs else _pdb_file_base

        try:
            store.read_input_file(_pdb_file_base, _pdb_file)
        except ClientError:
            continue

        if os.path.isfile(_pdb_file):
            pdb_file_base = _pdb_file_base
            pdb_file = _pdb_file
            break
    else:
        from Bio.PDB.PDBList import PDBList
        import gzip

        # obsolete = False
        # pdb_file = PDBList().retrieve_pdb_file(pdb, pdir=work_dir, file_format="pdb")
        # if not os.path.isfile(pdb_file):
        #     obsolete = True
        #     pdb_file = PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=work_dir, file_format="pdb")
        #     if not os.path.isfile(pdb_file):
        #         raise IOError("{} not found".format(r))

        for file_format, obs in (("pdb", False), ("mmCif", False), ("pdb", True), ("mmCif", True)):
            pdb_file = PDBList().retrieve_pdb_file(pdb, obsolete=obs, pdir=work_dir, file_format=file_format)
            if os.path.isfile(pdb_file):
                obsolete = obs
                format = file_format
                break
        else:
            raise IOError("{} not found".format(pdb))

        if file_format == "mmcif":
            pdb_file_base = "{}.cif.gz".format(pdb_file[:-7])
            pdb_file = "{}.cif.gz".format(pdb_file[:-7])

        with open(pdb_file, 'rb') as f_in, gzip.open(pdb_file+'.gz', 'wb') as f_out:
            f_out.writelines(f_in)

        store.write_output_file(pdb_file+".gz", "{}{}".format("obsolete/" if obsolete else "", pdb_file_base))

        try:
            os.remove(pdb_file+".gz")
        except OSError:
            pass

    return pdb_file, format
