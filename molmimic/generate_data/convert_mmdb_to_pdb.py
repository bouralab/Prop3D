import os, sys
import shutil
import string
from itertools import izip, product

import numpy as np
import pandas as pd
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

from molmimic.generate_data.map_residues import mmdb_to_pdb_resi, InvalidSIFTS

data_path_prefix = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

NUM_WORKERS = 4

def convert_row(job, row):
    resi = [str(row["from"]), str(row["to"])]
    try:
        frm, to = list(mmdb_to_pdb_resi(row["pdbId"], row["chnLett"], resi, replace_nulls=True))
    except (IOError, TypeError, ValueError, InvalidSIFTS) as e:
        job.log("Error mapping mmdb for {} '{}' {} {}".format(row["pdbId"], row["chnLett"], resi, e))
        return None

    return pd.Series({"from":frm, "to":to})

def mmdb2pdb(job, mmdb_file=None, pdb_group=None, dask=True, cores=NUM_WORKERS, preemtable=True):
    work_dir = job.fileStore.getLocalTempDir()

    if mmdb_file is None:
        prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
        jobStore = IOStore.get("{}:molmimic-ibis".format(prefix))
        mmdb_path = os.path.join(work_dir, "MMDB.h5")
        jobStore.read_input_file("MMDB.h5", mmdb_path)

    sdoms = pd.read_hdf(unicode(mmdb_path), "StrucutralDomains")

    if pdb_group is not None:
        tmp_path = os.path.join(work_dir, "{}.h5".format(pdb_group))
        pdb_groups = sdoms.groupby(lambda x: sdoms["pdbId"].loc[x][1:3])
        try:
            sdoms = pdb_groups.get_group(pdb_group.upper()).copy()
        except KeyError:
            return
            del pdb_groups
    else:
        tmp_path = os.path.join(work_dir, "PDB.h5")

    if dask:
        ddf = dd.from_pandas(sdoms, npartitions=NUM_WORKERS)
        meta = pd.DataFrame({"from":[1], "to":[1]})

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

    return job.fileStore.writeGlobalFile(tmp_path)

def merge(job, jobStoreIDs, cores=4):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    jobStore = IOStore.get("{}:molmimic-ibis".format(prefix))

    output = os.path.join(work_dir, "PDB.h5")

    for jobStoreID in jobStoreIDs:
        if jobStoreID is None: continue

        pdb_group_path = job.fileStore.readGlobalFile(jobStoreID)

        try:
            df = pd.read_hdf(unicode(pdb_group_path), "table")
        except IOError:
            continue

        df[['from', 'to']] = df[['from', 'to']].astype(int)
        df.to_hdf(unicode(output), "StructuralDomains", table=True, mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)
        del df

        job.fileStore.deleteGlobalFile(jobStoreID)

    sfams = pd.read_hdf(unicode(os.path.join(data_path_prefix, "MMDB.h5")), "Superfamilies")
    sfams = df[['sdsf_id', 'sdi', 'sfam_id', 'label', 'whole_chn', 'mmdb_id', 'mol_id', 'pssm']]
    sfams.to_hdf(unicode(output), "Superfamilies", format="table", table=True, complevel=9, complib="bzip2", min_itemsize=768)

    sdoms = pd.read_hdf(unicode(output+".tmp"), "StructuralDomains")
    merged = pd.merge(sdoms, sfams, how="left", on="sdi").dropna()
    merged.to_hdf(unicode(output), "merged", format="table", table=True, complevel=9, complib="bzip2", min_itemsize=768)

    resolu = pd.read_table("ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx",
        header=None, names=["pdbId", "resolution"], skiprows=6, sep="\t;\t")
    resolu.to_hdf(unicode(output+".tmp"), "resolu", format="table", table=True, complevel=9, complib="bzip2", min_itemsize=768)

    jobStore.write_output_file(output, "PDB.h5")

def start_toil(job, name="mmdb2pdb"):
    """Initiate Toil Jobs for converting MMDB to PDB

    Paramters
    ---------
    job : toil.Job

    Return
    ------
    mmdb2pdb : toil.Job
    """
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    jobStore = IOStore.get("{}:molmimic-ibis".format(prefix))

    if jobStore.exists("PDB.h5"):
        return

    try:
        mmdb_file = jobStore.write_input_file("MMDB.h5")
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise RuntimeError("download_data must be run before mmdb2pdb")

    group_files = [job.addChildJobFn(mmdb2pdb, mmdb_file=mmdb_file, pdb_group=group, cores=
        job.fileStore.jobStore.config.maxCores).rv() for group in \
        product(string.digits+string.ascii_lowercase, repeat=2)]
    job.addFollowOnJobFn(merge, group_files)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"

    print "Running"

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
