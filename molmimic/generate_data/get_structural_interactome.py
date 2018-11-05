import os, sys
import glob
from multiprocessing.pool import ThreadPool
import shutil
import itertools as it

import numpy as np
import pandas as pd
import dask
import dask.dataframe as dd
from joblib import Parallel, delayed

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(*args, **kwds):
        return args

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import iter_unique_superfams, get_file
from molmimic.generate_data.job_utils import map_job, map_job_rv, map_job_rv_list
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS

dask.config.set(scheduler='multiprocessing', num_workers=4)
dask.config.set(pool=ThreadPool(4))

def process_observed_interaction(job, int_id, ibisFileStoreID, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    mmdb_path = get_file(job, "PDB.h5", pdbFileStoreID)

    ibis_path = get_file(job, "IBIS_obs.h5", ibisFileStoreID)
    row = pd.read_hdf(ibis_path, "ObsInt", where="obs_int_id={}".format(int_id))

    try:
        #Read in face1 residues
        face1 = pd.read_hdf(unicode(ibis_path), "MolResFace", where="obs_int_id=row['obs_int_id']")
        face1.columns = ["obs_int_id", "mol_res"]
        print face1.shape
        #Keep entries from current CDD
        row = pd.merge(row, face1, how="left", on="obs_int_id")
        del face1

        #Read in face2 residues and convert gzipped asn1 into res numbers
        face2 = pd.read_hdf(unicode(ibis_path), "IntResFace", where="obs_int_id=row['obs_int_id']")
        face2.columns = ["obs_int_id", "int_res"]
        print face2.shape

        #Keep entries from current CDD
        row = pd.merge(row, face2, how="left", on="obs_int_id")
        print row.shape
        del face2

        st_domain_mol = pd.read_hdf(unicode(mmdb_path), "StructuralDomains", where="sdi=row['mol_sdi_id']")
        st_domain_mol.columns = ['mol_sdi_id', 'mol_domNo', 'mol_gi', 'mol_pdb', 'mol_chain', 'mol_sdi_from', 'mol_sdi_to']
        row = pd.merge(row, st_domain_mol, how="left", on="mol_sdi_id")
        del st_domain_mol

        st_domain_int = pd.read_hdf(unicode(mmdb_path), "StructuralDomains", where="sdi=row['int_sdi_id']")
        st_domain_int.columns = ['int_sdi_id', 'int_domNo', 'int_gi', 'int_pdb', 'int_chain', 'int_sdi_from', 'int_sdi_to']
        row = pd.merge(row, st_domain_int, how="left", on="int_sdi_id")
        del st_domain_int

        updated_resi = {"mol_res":[], "int_res": []}

        for resi in row.itertuples():
            try:
                updated_resi["mol_res"].append(decode_residues(resi.mol_pdb, resi.mol_chain, resi.mol_res, resi))
                updated_resi["int_res"].append(decode_residues(resi.int_pdb, resi.int_chain, resi.int_res, resi))
            except InvalidSIFTS:
                #This Row has failed return None
                return None
        if len(updated_resi["mol_res"]) > 0:
            row = row.assign(**updated_resi)

        path = job.fileStore.getLocalTempFileName()
        row.to_hdf(path, "table", table=True, format="table", complib="bzip2", complevel=9, min_itemsize=1024)
        rowFileStore = job.fileStore.writeGlobalFile(path)
        print "Done row", int_id
        return rowFileStore
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        raise
        return None

def merge_interactome_rows(job, sfam_id, converted_residues, inferred_table=None):
    print "Start merge", sfam_id
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    status =  "observed"
    new_cols = ["mol_res", "int_res"]
    resi_prefix = "{}.observed_interactome".format(int(sfam_id))
    data_cols = ["obs_int_id", "mol_sdi", "int_sdi"]

    resi_path = os.path.join(work_dir, resi_prefix)
    job.log("CON RESI PROMISES: {}".format(converted_residues))

    if len(converted_residues) == 0:
        job.log("FAILED {} no converted_residues".format(resi_path))
        print "FAILED {} no converted_residues".format(resi_path)
        return

    #Combine residues into dataframe
    for conv_store in it.chain.from_iterable(converted_residues):
        job.log("Running {} {}".format(conv_store, type(conv_store)))
        if not conv_store: continue
        conv_file = job.fileStore.readGlobalFile(conv_store)
        df = pd.read_hdf(conv_file, "table")
        try:
            df.to_hdf(unicode(resi_path), "table", table=True, format="table", mode="a",
                data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)
            job.fileStore.deleteGlobalFile(conv_store)
        except TypeError:
            job.log("Failed writing {}: {}".format(sfam_id, resi_path))
            try:
                os.remove(resi_path)
            except OSError:
                pass
            raise

    if os.path.isfile(resi_path):
        #Upload to S3
        out_store.write_output_file(resi_path, os.path.join(str(int(sfam_id)), resi_prefix))

        #Cleanup
        os.remove(resi_path)
        print "End merge", sfam_id
    else:
        job.log("Failed merging: {}".format(resi_path))
        print "Failed merging: {}".format(resi_path)

def get_observed_structural_interactome(job, sfam_id, pdbFileStoreID, ibisObsFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    ibis_obs_path = get_file(job, "IBIS_obs.h5", ibisObsFileStoreID)
    df = pd.read_hdf(unicode(ibis_obs_path), "ObsInt", where=["mol_superfam_id=={}".format(float(sfam_id)), "mol_superfam_id=={}".format(sfam_id)])
    int_ids = df["obs_int_id"].drop_duplicates().dropna()
    if len(int_ids) == 0:
        job.log("FAILED OBS SFAM {}".format(sfam_id))
        print "FAILED OBS SFAM {}".format(sfam_id)
        return

    #Add jobs for each interaction
    rows = map_job_rv(process_interaction, int_ids, ibisObsFileStoreID, pdbFileStoreID)
    #rows = [job.addChildJobFn(process_interaction, int_id, ibisObsFileStoreID, pdbFileStoreID).rv() for int_id in int_ids]
    job.log("{}".format(rows))
    #Merge converted residues
    job.addFollowOnJobFn(merge_interactome_rows, sfam_id, rows)

def start_toil_observed(job, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    ibis_obs_prefix = "IBIS_observed.h5"
    ibis_obs_path = os.path.join(work_dir, ibis_obs_prefix)
    in_store.read_input_file(ibis_obs_prefix, ibis_obs_path)

    #Add ibis info into local job store
    ibisObsFileStoreID = job.fileStore.writeGlobalFile(ibis_obs_path)

    print "start obs 1"
    pdb_path = job.fileStore.readGlobalFile(pdbFileStoreID)
    skip_sfam = set([float(f.key.split("/", 1)[0]) for f in out_store.list_input_directory(None) \
        if f.key.endswith(".observed_interactome")])
    job.log("start obs 2")
    pdb = pd.read_hdf(unicode(pdb_path), "Superfamilies", columns=["sfam_id"])
    sfams = pdb[~pdb["sfam_id"].isin(skip_sfam)]["sfam_id"].drop_duplicates().dropna().astype(int)
    print "SFAMS", sfams
    map_job(job, get_observed_structural_interactome, sfams, pdbFileStoreID, ibisObsFileStoreID)
    print "Finished adding jobs"
    os.remove(pdb_path)
    print "FINISHED obs"

def start_toil(job):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    #Download PDB info
    pdb_file = os.path.join(work_dir, "PDB.h5")
    in_store.read_input_file("PDB.h5", pdb_file)

    #Add pdb info into local job store
    pdbFileStoreID = job.fileStore.writeGlobalFile(pdb_file)

    obsjob = job.addChildJobFn(start_toil_observed, pdbFileStoreID)
    #infjob = obsjob.addFollowOnJobFn(start_toil_inferred, pdbFileStoreID)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    print "Running"

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
