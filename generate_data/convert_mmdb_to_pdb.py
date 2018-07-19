import os, sys
import shutil
from itertools import izip

import numpy as np
import pandas as pd
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

from map_residues import mmdb_to_pdb_resi, InvalidSIFTS

data_path_prefix = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

NUM_WORKERS = 4

def convert_row(job, row):
    resi = [str(row["from"]), str(row["to"])]
    try:
        frm, to = list(mmdb_to_pdb_resi(row["pdbId"], row["chnLett"], resi, replace_nulls=True))
    except (IOError, TypeError, ValueError, InvalidSIFTS) as e:
        job.log("Error mapping mmdb for {} '{}' {} {}".format(row["pdbId"], row["chnLett"], resi, e))
        frm, to = np.nan, np.nan

    return pd.Series({"from":frm, "to":to})

def mmdb2pdb(job, pdb_group=None, dask=False, cores=NUM_WORKERS):
    mmdb_path = os.path.join(data_path_prefix, "MMDB.h5")
    cols = ["from", "to"]
    sdoms = pd.read_hdf(mmdb_path, "StrucutralDomains")

    if pdb_group is not None:
        tmp_path = os.path.join(data_path_prefix, "__mmdb_to_pdb", "{}.h5".format(pdb_group))
        pdb_groups = sdoms.groupby(lambda x: sdoms["pdbId"].loc[x][1:3])
        try:
            sdoms = pdb_groups.get_group(pdb_group.upper()).copy()
        except KeyError:
            return
            del pdb_groups
    else:
        tmp_path = os.path.join(data_path_prefix, "PDB.h5")

    if dask:
        ddf = dd.from_pandas(sdoms, npartitions=NUM_WORKERS)
        meta = pd.DataFrame({c:[1] for c in cols})

        with ProgressBar():
            pdb_map = ddf.map_partitions(lambda _df: _df.apply(\
                lambda row: convert_row(job, row), axis=1), meta=meta)\
                    .compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    else:
        pdb_map = sdoms.apply(lambda row: convert_row(job, row), axis=1)

    sdoms.loc[:, ["from", "to"]] = pdb_map
    sdoms = sdoms.dropna()
    sdoms["from"] = sdoms["from"].astype(int)
    sdoms["to"] = sdoms["to"].astype(int)
    sdoms.to_hdf(unicode(tmp_path), "table", complevel=9, complib="bzip2", min_itemsize=756)

    print "Finished StructuralDomains"

def merge(job, cores=4):
    pdb_dir = os.path.join(data_path_prefix, "pdb", "pdb")
    output = os.path.join(data_path_prefix, "PDB.h5")
    tmp_path = os.path.join(data_path_prefix, "__mmdb_to_pdb")
    for pdb_group in next(os.walk(pdb_dir))[1]:
        job.log(pdb_group)
        pdb_group = os.path.basename(pdb_group)
        pdb_group_path = os.path.join(tmp_path, "{}.h5".format(pdb_group.lower()))
        try:
            df = pd.read_hdf(pdb_group_path, "table")
        except IOError:
            continue
        df[['from', 'to']] = df[['from', 'to']].astype(int)
        df.to_hdf(unicode(output+".tmp"), "StructuralDomains", table=True, mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)
        del df

    sfams = pd.read_hdf(os.path.join(data_path_prefix, "MMDB.h5"), "Superfamilies")
    sfams = df[['sdsf_id', 'sdi', 'sfam_id', 'label', 'whole_chn', 'mmdb_id', 'mol_id', 'pssm']]
    sfams.to_hdf(unicode(output+".tmp"), "Superfamilies", complevel=9, complib="bzip2", min_itemsize=768)

    sdoms = pd.read_hdf(output+".tmp", "StructuralDomains")
    merged = pd.merge(sdoms, sfams, how="left", on="sdi").dropna()
    merged.to_hdf(unicode(output+".tmp"), "merged", complevel=9, complib="bzip2", min_itemsize=768)

    shutil.move(output+".tmp", output)
    shutil.rmtree(tmp_path)

# def submit_jobs(job_name="mmdb2pdb"):
#     pdb_dir = os.path.join(data_path_prefix, "pdb", "pdb")
#     tmp_path = os.path.join(data_path_prefix, "__mmdb_to_pdb")
#     if os.path.isdir(tmp_path):
#         shutil.rmtree(tmp_path)
#     os.makedirs(tmp_path)
#     job = SlurmJob(job_name, walltime="3-00:00:00")
#     for pdb_group in next(os.walk(pdb_dir))[1]:
#         job += "{} {} {}\n".format(sys.executable, __file__, os.path.basename(pdb_group))
#     jid = job.run()
#
#     job2 = SlurmJob(job_name+"_merge", walltime="3:00:00", individual=True)
#     job2 += "{} {} merge\n".format(sys.executable, __file__)
#     jid2 = job2.submit_individual(dependency="afterany:"+jid)
#     return jid2

def start_toil(job, name="mdb2pdb"):
    """Initiate Toil Jobs for converting MMDB to PDB

    Paramters
    ---------
    job : toil.Job

    Return
    ------
    mmdb2pdb : toil.Job
    """
    if os.path.isfile(os.path.join(data_path_prefix, "PDB.h5")):
        return

    pdb_dir = os.path.join(data_path_prefix, "pdb", "pdb")
    tmp_path = os.path.join(data_path_prefix, "__mmdb_to_pdb")
    if os.path.isdir(tmp_path):
        shutil.rmtree(tmp_path)
    os.makedirs(tmp_path)

    for pdb_group in next(os.walk(pdb_dir))[1]:
        pdb_group = os.path.basename(pdb_group)
        job.addChildJobFn(mmdb2pdb, pdb_group=pdb_group)

    job.addFollowOnJobFn(merge)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        submit_jobs()
    elif len(sys.argv) == 2 and sys.argv[1] == "merge":
        merge()
    elif len(sys.argv) == 2:
        convert_pdb_group(sys.argv[1])
