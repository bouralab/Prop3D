from __future__ import print_function
import os
import sys
import re
import subprocess
from contextlib import contextmanager
from itertools import izip

import pandas as pd
import numpy as np
import toil.fileStore

#Auto-scaling on AWS with toil has trouble finding modules? Heres the workaround
PDB_TOOLS = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdb_tools")

def SubprocessChain(commands, output):
    if len(commands) > 2:
        prev_proc = subprocess.Popen(
            commands[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=os.environ)
        for cmd in commands[1:-1]:
            proc = subprocess.Popen(
                cmd,
                stdin=prev_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=os.environ)
            prev_proc = proc
        final_proc = subprocess.Popen(
            commands[-1],
            stdin=prev_proc.stdout,
            stdout=output,
            stderr=subprocess.PIPE,
            env=os.environ)
        return final_proc.communicate()
    elif len(commands) == 2:
        prev_proc = subprocess.Popen(
            commands[0],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=os.environ)
        final_proc = subprocess.Popen(
            commands[1],
            stdin=prev_proc.stdout,
            stdout=output,
            stderr=subprocess.PIPE,
            env=os.environ)
        print("Ran 2 commands", commands[1])
        return final_proc.communicate()
    elif len(commands) == 1:
        final_proc = subprocess.Popen(
            commands[0],
            stdout=output,
            stderr=subprocess.PIPE,
            env=os.environ)
    else:
        raise RuntimeError
    return final_proc.communicate()

def get_file(job, prefix, path_or_fileStoreID, work_dir=None):
    if isinstance(path_or_fileStoreID, str) and os.path.isfile(path_or_fileStoreID):
        return path_or_fileStoreID
    else:
        work_dir = work_dir or job.fileStore.getLocalTempDir()
        new_file = os.path.join(work_dir, prefix)

        if isinstance(path_or_fileStoreID, (toil.fileStore.FileID, str)):
            with job.fileStore.readGlobalFileStream(path_or_fileStoreID) as fs, open(new_file, "w") as nf:
                for line in fs:
                    nf.write(line)
        elif hasattr(path_or_fileStoreID, "read_input_file"):
            #Might be file store itself
            path_or_fileStoreID.read_input_file(prefix, new_file)

        else:
            raise RuntimeError("Invalid path_or_fileStoreID {} {}".format(type(path_or_fileStoreID), path_or_fileStoreID))

        return new_file

def filter_hdf(hdf_path, dataset, column=None, value=None, columns=None, **query):
    assert len(query) > 0 or (column is not None and value is not None)
    if len(query) == 0 and (column is not None and value is not None):
        if not isinstance(column, (list, tuple)) or not isinstance(value, (list, tuple)):
            where = "{}={}".format(column, value)
        elif len(column) == len(value):
            where = ["{}={}".format(c, v) for c, v in zip(column, value)]
        else:
            raise RuntimeError("Cols and values must match")
    else:
        where = ["{}={}".format(c,v) for c, v in query.iteritems()]
    try:
        df = pd.read_hdf(unicode(hdf_path), dataset, where=where, columns=columns)
        if df.shape[0] == 0: raise KeyError
    except (KeyError, ValueError, SyntaxError):
        df = filter_hdf_chunks(hdf_path, dataset, column=column, value=value, columns=columns, **query)
    return df

def filter_hdf_chunks(hdf_path, dataset, column=None, value=None, columns=None, chunksize=500, **query):
    df = None
    for _df in pd.read_hdf(unicode(hdf_path), dataset, chunksize=chunksize):
        if len(query) > 0:
            try:
                filtered_df = _df.query("("+") & (".join(["{}=={}".format(c, v) for c, v in query.iteritems()])+")")
            except SyntaxError:
                expression = None
                for c, v in query.iteritems():
                    _exp = _df[c]==v
                    if expression is None:
                        expression = _exp
                    else:
                        expression &= _exp
                filtered_df = [expression]
        elif column is not None and value is not None:
            if not isinstance(column, (list, tuple)) and not isinstance(value, (list, tuple)):
                filtered_df = _df[_df[column]==value].copy()
            elif len(column) == len(value):
                filtered_df = _df.query(" & ".join(["{}=={}".format(c, v) for c, v in zip(column, value)]))
            else:
                raise RuntimeError("Cols and values must match")
        else:
            filtered_df = _df.copy()
        if filtered_df.shape[0]>0:
            if columns is not None:
                filtered_df = filtered_df[columns]
            df = pd.concat((df, filtered_df), axis=0) if df is not None else filtered_df
            del _df
            _df = None
    if df is None:
        raise TypeError("Unable to parse HDF")
    return df

def PDBTools(commands, output):
    cmds = [[sys.executable, "-m", "pdb-tools.pdb_{}".format(cmd[0])]+cmd[1:] \
        for cmd in commands]
    SubprocessChain(cmd, output)

def get_jobstore_name(job, name="raw-files"):
    return "{}-{}".format(job._fileStore.jobStore.locator, name)

def get_jobstore(job, name="raw-files"):
    return Toil.getJobStore(get_jobstore_name(job, None))

def iter_cdd(use_label=True, use_id=False, label=None, id=None, group_superfam=False, all_superfam=False):
    if use_label and not use_id:
        col = 1
    elif not use_label and use_id:
        col = 2

    CDD = pd.read_hdf(unicode(os.path.join(data_path_prefix, "MMDB.h5")), "Superfamilies")
    CDD = CDD[["label", "sfam_id"]].drop_duplicates().dropna()

    if label is not None:
        CDD = CDD[CDD["label"]==label]
    elif id is not None:
        CDD = CDD[CDD["sfam_id"]==id]

    CDD["label"] = CDD["label"].apply(lambda cdd: cdd.replace("/", "").replace("'", "\'") if isinstance(cdd, str) else cdd)
    CDD.sort_values("label", inplace=True)

    if group_superfam:
        groups = CDD.groupby("sfam_id")
        if all_superfam:
            if use_label and use_id:
                for superfam_id, families in groups:
                    yield families.itertuples(), superfam_id
            else:
                for superfam_id, families in groups:
                    yield families.itertuples()
        else:
            if use_label and use_id:
                for superfam_id, families in groups:
                    yield next(families.itertuples())[1], superfam_id
            else:
                for superfam_id, families in groups:
                    yield next(families.itertuples())[1]
    else:
        if use_label and use_id:
            for cdd in CDD.itertuples():
                yield cdd[1], cdd[2]
        else:
            for cdd in CDD.itertuples():
                yield cdd[col]

def iter_unique_superfams():
    CDD = pd.read_hdf(unicode(os.path.join(data_path_prefix, "MMDB.h5")), "Superfamilies")
    CDD = CDD["sfam_id"].drop_duplicates().dropna()
    for cdd in CDD:
        yield cdd
    del CDD

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
                if line.startswith("ATOM") and line[13:15] != "CA":
                    return False
        return True
    except IOError:
        return False

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

def get_pdb_residues(pdb_file):
    try:
        with open(pdb_file) as f:
            prev_res = None
            for line in f:
                res = natural_keys(line[22:27], use_int=True) #inlcude icode
                if line.startswith("ATOM") and not res == prev_res:
                    yield res
    except IOError:
        pass

def izip_missing(iterA, iterB, **kwds):
    """Iterate through two iterables, while making sure they are in the same
    order. If there are missing values, you can skip the value entirely or
    return only the iterator with the value and a special fill value.
    Parameters:
    ___________
    iterA : the first iterator
    iterB : the second iterator
    key : function that returns items to compare. Must return strings, ints, or
        an object with the __lt__, __gt__, and __eq__ methods. Optional.
    fillvalue : The value to return if the item is missing. Optional.
    Returns:
    ________
    A : item from first iterator, or fillValue
    B : item from second iterator, or fillValue
    """
    #Get the comparison functions
    key = kwds.get("key", lambda x: x)
    keyA = kwds.get("keyA", key)
    keyB = kwds.get("keyB", key)

    useMissing = "fillvalue" in kwds
    fillvalue = kwds.get("fillvalue")

    verbose = kwds.get("verbose", False)

    #Start both iterators
    A = next(iterA)
    B = next(iterB)
    try:
        while True:
            if keyA(A) == keyB(B):
                if verbose:
                    print(keyA(A), "==", keyB(B))
                yield A, B
                A = next(iterA)
                B = next(iterB)
            elif keyA(A) < keyB(B):
                if verbose:
                    print(keyA(A), "<", keyB(B))
                if useMissing:
                    yield(A, fillvalue)
                A = next(iterA)
            elif keyA(A) > keyB(B):
                if verbose:
                    print(keyA(A), ">", keyB(B))
                if useMissing:
                    yield fillvalue, B
                B = next(iterB)
            else:
                raise RuntimeError("Invalid comparator")
    except StopIteration:
        pass

def make_h5_tables(files, iostore):
    for f in files:
        iostore.read_input_file(f, f)
        store = pd.HDFStore(unicode(f))
        for key in store.keys():
            df = store.get(key)
            df.to_hdf(unicode(f+".new"), "table", format="table",
                table=True, complevel=9, complib="bzip2", min_itemsize=1024)
        iostore.write_output_file(f+".new", f)

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
        replace = dict(izip(chains, rename))
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
        for old_line, new_line in izip(
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
