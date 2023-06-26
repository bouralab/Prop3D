import os
import sys
import re
import shutil
import subprocess
import glob
import functools

import numpy as np
import pandas as pd
import requests
from botocore.exceptions import ClientError

from Prop3D.util import SubprocessChain, natural_keys
from Prop3D.util.iostore import IOStore

#Auto-scaling on AWS with toil has trouble finding modules? Heres the workaround
PDB_TOOLS = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdb_tools")

class InvalidPDB(RuntimeError):
    pass

def PDBTools(commands, output):
    cmds = [[sys.executable, "-m", "pdb-tools.pdb_{}".format(cmd[0])]+cmd[1:] \
        for cmd in commands]
    SubprocessChain(cmd, output)

class PDBTools(object):
    def __init__(self, pdb_file):
        self.pdb = pdb_file

    def __get__(self, command, args):
        available_commands = [os.path.basename(f)[:-3] for f in \
            os.listdir(PDB_TOOLS) if f.endswith(".py")]
        if command in available_commands or "pdb_"+command in available_commands:
            return self.run(cmd, *args)
        return super().__get__(command, args)

    def run(self, cmd, *args):
        cmd = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_{}.py".format(cmd))]]
        with open("{}_{}.out".format(pdb_file, cmd)) as f:
            SubprocessChain(cmd, f)

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

def get_pdb_residues(pdb_file, include_resn=False, use_raw_resi=False):
    try:
        with open(pdb_file) as f:
            prev_resi = None
            prev_resn = None
            for line in f:
                if line.startswith("ATOM"):
                    resi_raw = line[22:27] #inlcude icode
                    resi = natural_keys(resi_raw, use_int=True)
                    resn = line[17:20]
                    if not resi == prev_resi:
                        if include_resn:
                            yield resn, resi_raw.strip() if use_raw_resi else resi
                        else:
                            yield resi_raw.strip() if use_raw_resi else resi
                    prev_resi = resi
                    prev_resn = resn
    except IOError:
        pass

def get_b(pdb_file):
    for line in get_atom_lines(pdb_file):
        yield float(line[60:66].strip())

def read_pdb(file):
    return np.array([(float(line[30:38]), float(line[38:46]), float(line[46:54])) \
        for line in get_atom_lines(file)])

def replace_chains(pdb_file, new_file=None, new_chain=None, **chains):
    """Modified from pdbotools"""
    assert new_chain is not None or len(chains)>0
    coord_re = re.compile('^(ATOM|HETATM)')
    if new_file is None:
        name, ext = os.path.splitext(pdb_file)
        new_file = "{}.{}.pdb".format(name, new_chain if new_chain is not None \
            else "".join(chains.keys()))

    with open(pdb_file) as f, open(new_file, "w") as new:
        for line in f:
            line = line.rstrip()
            if coord_re.match(line) and line[21] in chains:
                print(line[:21] + chains[line[21]] + line[22:], file=new)
            elif coord_re.match(line) and new_chain is not None and isinstance(new_chain, str):
                print(line[:21] + new_chain + line[22:], file=new)
            else:
                print(line, file=new)

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

def rottrans_from_symop(moving_pdb, symmetry_operation, new_file=None):
    if symmetry_operation in [None, "X,Y,Z"]:
        return moving_pdb

    sym_ops = [dim.split("+") for dim in symmetry_operation.split(",")]
    rot, trans = zip(*sym_ops)
    trans = map(float, trans)

    rot_mat = np.eye(3)
    for i, (r,d) in enumerate(zip(rot, ("X","Y","Z"))):
        rot_mat[i,i] = float(r.replace(d, "1"))

    return rottrans_from_matrix(moving_pdb, rot_mat, trans, new_file=new_file)

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

    return updated_pdb

def reres_pdb(pdb_file, resid=1, chain=None, updated_pdb=None):
    if updated_pdb is None:
        updated_pdb = "{}.reres.pdb".format(os.path.splitext(pdb_file)[0])

    cmd = [sys.executable, os.path.join(PDB_TOOLS, "pdb_reres.py")] #, pdb_file]

    if chain is not None:
        if chain == "all":
            cmd += ["-chain"]
        else:
            cmd += ["-chain", str(chain)]

    cmd += ["-resid", str(resid)]

    with open(pdb_file) as p, open(updated_pdb, "w") as f:
        out, err = subprocess.Popen(cmd, stdin=p, stdout=f, stderr=subprocess.PIPE).communicate()
        if err:
            raise RuntimeError(err)

    return updated_pdb

def remove_ter_lines(pdb_file, updated_pdb=None):
    if updated_pdb is None:
        updated_pdb = "{}.untidy.pdb".format(os.path.splitext(pdb_file)[0])

    with open(pdb_file) as f, open(updated_pdb, "w") as updated:
        for line in f:
            if not line.startswith("TER"):
                updated.write(line)

    return updated_pdb

def split_xyz(pdb_file, updated_pdb=None):
    if updated_pdb is None:
        updated_pdb = "{}.xyz.pdb".format(os.path.splitext(pdb_file)[0])

    with open(pdb_file) as f, open(updated_pdb, "w") as updated:
        for line in f:
            if line.startswith("ATOM"):
                new = line[:30]+" "+line[30:38]+" "+line[38:46]+" "+line[46:53]+" "+line[53:]
                updated.write(new)
                print(new)
            else:
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
    if os.path.isfile("bfactor.norm"):
        with open("bfactor.norm") as f:
            mean, std = map(float, f.read().split())
    else:
        raise RuntimeError("Must call get_bfactor_norm_mean before this")

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
            pdb_bfactors = "\n".join(line[60:66].strip() for line in \
                get_atom_lines(pdb_file))
            print(pdb_bfactors, file=fh)

    bfactors = np.loadtxt(bfactors_file)
    raw_mean, raw_std = bfactors.mean(), bfactors.std()
    mad = median_absolute_deviation(bfactors)

    #From Chung, Wang, Bourne. "Exploiting sequence and structure homologs to identify
    #protein–protein binding sites." Proteins: Structure, Function, Bioinformatics. 2005.
    #https://doi.org/10.1002/prot.20741
    M = 0.6745*((bfactors-raw_mean)/mad)
    bfactors = bfactors[M<=3.5]
    mean, std = bfactors.mean(), bfactors.std()

    with open(bfactor_norm_file, "w") as fh:
        print("{} {}".format(mean, std), file=fh)

    return mean, std

def s3_download_pdb(pdb, work_dir=None, remove=False, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    from Prop3D.generate_data.data_stores import data_stores
    from Prop3D.parsers.pdbe import PDBEApi

    store = data_stores(job).raw_pdb_store

    if len(pdb) >= 6 and pdb[4] == ".":
        pdb, chain = pdb.split(".")
    if len(pdb) != 4:
        #UniProt -> get alphafold or best pdb
        uniprot_id = pdb.upper()
        uniprot_file = os.path.join(work_dir, "{}.pdb".format(uniprot_id))
        alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.pdb"
        s3_file = "AlphaFold/"+uniprot_file
        chain = "A" #Default chain for all AlphaFOld files

        try:
            store.read_input_file(s3_file, uniprot_file)
            if not os.path.isfile(uniprot_file):
                raise RuntimeError
            return uniprot_file, chain, "pdb"
        except (ClientError, RuntimeError):
            try:
                with requests.get(alphafold_url, stream=True) as r:
                    r.raw.read = functools.partial(r.raw.read, decode_content=True)
                    with open(uniprot_file, 'wb') as f:
                        shutil.copyfileobj(r.raw, f)
                store.write_output_file(uniprot_file, s3_file)
                return uniprot_file, "pdb"
            except requests.RequestException:
                #Cannot download, check pdb
                pdbe = PDBEApi(work_dir=work_dir)
                try:
                    mappings = pdbe.get_uniprot_mapping(uniprot_id)
                    if uniprot_id not in mappings or "PDB" not in mappings[uniprot_id]:
                        raise KeyError
                    pdb = min(mappings[uniprot_id]["PDB"].keys())
                    chain = mappings[uniprot_id]["PDB"][pdb]["chain_id"]
                except KeyError:
                    if len(pdb)==5:
                        pdb, chain = pdb[:4], pdb[5]
                    else:
                        raise RuntimeError(f"Cannot find PDB file or AlphaFoldDB file for {pdb}")

    pdb_file_base = os.path.join(pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    pdb_file = os.path.join(work_dir, "pdb{}.ent.gz".format(pdb.lower()))

    format = "pdb"

    for file_format, obs in (("pdb", False), ("mmCif", False), ("pdb", True), ("mmCif", True)):
        _pdb_file_base = pdb_file_base.replace(".ent.", ".mmcif.") if file_format == "mmcif" else pdb_file_base
        _pdb_file = pdb_file.replace(".ent.", ".mmcif.") if file_format == "mmcif" else pdb_file
        _pdb_file_base = "obsolete/"+_pdb_file_base if obs else _pdb_file_base

        try:
            store.read_input_file(_pdb_file_base, _pdb_file)
        except (ClientError, RuntimeError):
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

    if chain is None:
        chain = get_first_chain(pdb_file)

    return pdb_file, chain, format
