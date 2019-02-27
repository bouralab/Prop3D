import os, sys
import glob
from multiprocessing.pool import ThreadPool
import shutil
import itertools as it

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(*args, **kwds):
        return args

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import iter_unique_superfams, get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job, map_job_rv, map_job_rv_list
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS

from toil.realtimeLogger import RealtimeLogger

def process_observed_interaction(job, int_id, sfam_id, ibisFileStoreID, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    mmdb_path = get_file(job, "PDB.h5", pdbFileStoreID)

    ibis_path = get_file(job, "IBIS_obs.h5", ibisFileStoreID)
    row = filter_hdf(ibis_path, "ObsInt", "obs_int_id", int_id)

    try:
        #Read in face1 residues
        face1 = filter_hdf(str(ibis_path), "MolResFace", "obs_int_id", int_id)
        face1.columns = ["obs_int_id", "mol_res"]
        #Keep entries from current CDD
        row = pd.merge(row, face1, how="left", on="obs_int_id")
        del face1

        #Read in face2 residues and convert gzipped asn1 into res numbers
        face2 = filter_hdf(str(ibis_path), "IntResFace", "obs_int_id", int_id)
        face2.columns = ["obs_int_id", "int_res"]

        #Keep entries from current CDD
        row = pd.merge(row, face2, how="left", on="obs_int_id")
        del face2

        try:
            st_domain_mol = filter_hdf(str(mmdb_path), "StructuralDomains", "sdi", row.iloc[0]['mol_sdi_id'])
            st_domain_mol.columns = ['mol_sdi_id', 'mol_domNo', 'mol_gi', 'mol_pdb', 'mol_chain', 'mol_sdi_from', 'mol_sdi_to']
            row = pd.merge(row, st_domain_mol, how="left", on="mol_sdi_id")
            del st_domain_mol

            st_domain_int = filter_hdf(str(mmdb_path), "StructuralDomains", "sdi", row.iloc[0]['int_sdi_id'])
            st_domain_int.columns = ['int_sdi_id', 'int_domNo', 'int_gi', 'int_pdb', 'int_chain', 'int_sdi_from', 'int_sdi_to']
            row = pd.merge(row, st_domain_int, how="left", on="int_sdi_id")
            del st_domain_int
        except TypeError:
            #SDI's doesn't exists, must be obsolete
            RealtimeLogger.info("Done row", int_id)
            return

        updated_resi = {"mol_res":[], "int_res": []}

        for resi in row.itertuples():
            try:
                updated_resi["mol_res"].append(decode_residues(job, resi.mol_pdb, resi.mol_chain, resi.mol_res, resi))
                updated_resi["int_res"].append(decode_residues(job, resi.int_pdb, resi.int_chain, resi.int_res, resi))
            except InvalidSIFTS:
                #This Row has failed converting binary to string skip iti
                updated_resi["mol_res"].append(np.NaN)
                updated_resi["int_res"].append(np.NaN)
                continue

        if len(updated_resi["mol_res"].dropna()) > 0:
            row = row.assign(**updated_resi)
            row.dropna()
        else:
            #This entire interation failed; returned None
            return None

        path = "{}.h5".format(int_id)
        row.to_hdf(path, "table", table=True, format="table", complib="bzip2", complevel=9, min_itemsize=1024)
        out_store.write_output_file(path, "{}/_obsrows/{}".format(int(sfam_id), path))
        RealtimeLogger.info("Done row", int_id)
        try:
            os.remove(fail_file)
        except OSError:
            pass
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        job.log("FAILED {} {}".format(int_id, e, tb))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write(str(e))
            f.write(str(tb))
        out_store.write_output_file(fail_file, "{}/_obsrows/{}.failed".format(int(sfam_id), int_id))
        try:
            os.remove(fail_file)
        except OSError:
            pass

def merge_interactome_rows(job, sfam_id):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    status =  "observed"
    new_cols = ["mol_res", "int_res"]
    resi_prefix = "{}.observed_interactome".format(int(sfam_id))
    data_cols = ["obs_int_id", "mol_sdi", "int_sdi"]

    resi_path = os.path.join(work_dir, resi_prefix)

    #Combine residues into dataframe
    possible_errors = []
    nrows = None
    for nrows, row_prefix in enumerate(out_store.list_input_directory("{}/_obsrows/".format(int(sfam_id)))):
        if row_prefix.endswith("failed"): continue
        job.log("Running {} {}".format(sfam_id, row_prefix))

        row_file = os.path.join(work_dir, os.path.basename(row_prefix))
        out_store.read_input_file(row_prefix, row_file)

        df = pd.read_hdf(row_file, "table")
        try:
            df.to_hdf(str(resi_path), "table", table=True, format="table", append=True, mode="a",
                data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            import traceback
            tb = traceback.format_exc()
            job.log("Failed writing {}: {} {}".format(sfam_id, resi_path, tb))
            possible_errors.append(tb)
            continue

        try:
            os.remove(row_file)
        except OSError:
            pass

        out_store.remove_file(row_prefix)

    if os.path.isfile(resi_path):
        #Upload to S3
        out_store.write_output_file(resi_path, os.path.join(str(int(sfam_id)), resi_prefix))

        #Cleanup
        os.remove(resi_path)
    elif nrows is not None:
        RealtimeLogger.info("Failed merging: {}".format(resi_path))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write("No rows?")
            for e in possible_errors:
                f.write(e)
                f.write("\n")
        out_store.write_output_file(fail_file, "{}/{}.failed".format(int(sfam_id), resi_prefix))
        try:
            os.remove(fail_file)
        except OSError:
            pass

def get_observed_structural_interactome(job, sfam_id, pdbFileStoreID, ibisObsFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    ibis_obs_path = get_file(job, "IBIS_obs.h5", ibisObsFileStoreID)
    try:
        df = filter_hdf(ibis_obs_path, "ObsInt", "mol_superfam_id", float(sfam_id), columns=["obs_int_id"])
        int_ids = df["obs_int_id"].drop_duplicates().dropna().astype(int)
        if len(int_ids) == 0:
            RealtimeLogger.info("EMPTY OBS SFAM {}".format(sfam_id))
            return
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        RealtimeLogger.info("FAILED OBS SFAM {} {}".format(sfam_id, e))
        return

    current_rows = set(int(os.path.basename(key)[:-3]) for key in out_store.list_input_directory("{}/_obsrows".format(int(sfam_id))) if not key.endswith("failed"))
    int_ids = list(set(int_ids)-current_rows)
    RealtimeLogger.info("Will run {} ids: {}".format(len(int_ids), int_ids))

    if len(int_ids) > 0:
        #Add jobs for each interaction
        map_job(job, process_observed_interaction, int_ids, sfam_id, ibisObsFileStoreID, pdbFileStoreID)

    #Merge converted residues
    job.addFollowOnJobFn(merge_interactome_rows, sfam_id)

def cleanup(job):
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-interfaces".format(prefix))
    keys = list(in_store.list_input_directory())
    finished = set(key.split("/")[0] for key in keys if key.endswith("observed_interactome"))
    failed = set(key.split("/")[0] for key in keys if "_obsrows" in key and key.endswith("failed"))

    for key in failed-finished:
        in_store.remove_file(key)

def start_toil(job, pdbFileStoreID=None):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    if pdbFileStoreID is None:
        #Download PDB info
        pdb_path = os.path.join(work_dir, "PDB.h5")
        in_store.read_input_file("PDB.h5", pdb_path)

        #Add pdb info into local job store
        pdbFileStoreID = job.fileStore.writeGlobalFile(pdb_path)
    else:
        pdb_path = job.fileStore.readGlobalFile(pdbFileStoreID)

    ibis_obs_prefix = "IBIS_observed.h5"
    ibis_obs_path = os.path.join(work_dir, ibis_obs_prefix)
    in_store.read_input_file(ibis_obs_prefix, ibis_obs_path)

    #Add ibis info into local job store
    ibisObsFileStoreID = job.fileStore.writeGlobalFile(ibis_obs_path)

    #Choose which superfamilies to run, skip those already present
    skip_sfam = set([float(f.split("/", 1)[0]) for f in out_store.list_input_directory() \
        if f.endswith(".observed_interactome")])
    pdb = filter_hdf_chunks(str(ibis_obs_path), "ObsInt", columns=["mol_superfam_id"]).drop_duplicates()
    sfams = pdb[~pdb["mol_superfam_id"].isin(skip_sfam)]["mol_superfam_id"].drop_duplicates().dropna().astype(int)
    pRealtimeLogger.info("Will run a total of {} SFAMS".format(len(sfams)))

    #Run all superfamilies
    map_job(job, get_observed_structural_interactome, sfams, pdbFileStoreID, ibisObsFileStoreID)

    #Cleanup
    job.addFollowOnJobFn(cleanup)
    os.remove(ibis_obs_path)
    os.remove(pdb_path)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
