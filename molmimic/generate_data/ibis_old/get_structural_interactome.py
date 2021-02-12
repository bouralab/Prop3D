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

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError

def process_observed_interaction(job, int_id, sfam_id, ibisFileStoreID, pdbFileStoreID, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")

    sfam = int(sfam_id) if sfam_id is not None else sfam_id

    mmdb_path = get_file(job, "PDB.h5", pdbFileStoreID)

    ibis_path = get_file(job, "IBIS_obs.h5", ibisFileStoreID)
    row = filter_hdf(ibis_path, "ObsInt", "obs_int_id", int_id)

    #Remove GI since they are
    #row = row.drop(columns={"mol_gi"})

    row["mol_superfam_acc"] = row["mol_superfam_acc"].fillna("").str.encode("ascii")
    row["int_superfam_acc"] = row["int_superfam_acc"].fillna("").str.encode("ascii")

    try:
        #Read in face1 residues
        face1 = filter_hdf(str(ibis_path), "MolResFace", "obs_int_id", int_id)
        face1 = face1.rename(columns={"seqloc": "mol_res"})

        #Keep entries from current SDI
        row = pd.merge(row, face1, how="left", on="obs_int_id")
        del face1

        #Read in face2 residues and convert gzipped asn1 into res numbers
        face2 = filter_hdf(str(ibis_path), "IntResFace", "obs_int_id", int_id)
        face2 = face2.rename(columns={"seqloc": "int_res"})

        #Keep entries from current SDI
        row = pd.merge(row, face2, how="left", on="obs_int_id")
        del face2

        try:
            RealtimeLogger.info("Get SDI {}".formet(row))
            st_domain_mol = filter_hdf(str(mmdb_path), "StructuralDomains", "sdi", row.iloc[0]['mol_sdi_id'])
            RealtimeLogger.info("MOL ROW {}".format(st_domain_mol))
            st_domain_mol = st_domain_mol.rename(columns={"sdi":"mol_sdi_id",
                "domNo": "mol_domNo", "gi": "mol_gi", "pdbId": "mol_pdb",
                "chnLett": "mol_chain", "from": "mol_sdi_from", "to": "mol_sdi_to",
                "num_residues": "mol_nres", "whole_chain": "mol_whole_chain"})
            st_domain_mol = st_domain_mol.drop(columns=["mol_sdi_from", "mol_sdi_to"])
            st_domain_mol["mol_gi"] = st_domain_mol["mol_gi"].astype(float)
            st_domain_mol = st_domain_mol.drop_duplicates()
            RealtimeLogger.info("MOL ROW {}".format(st_domain_mol))
            row = pd.merge(row, st_domain_mol, how="left", on=["mol_sdi_id"])
            del st_domain_mol
            RealtimeLogger.info("MOL ROW {}".format(row))

            st_domain_int = filter_hdf(str(mmdb_path), "StructuralDomains", "sdi", row.iloc[0]['int_sdi_id'])
            st_domain_int = st_domain_int.rename(columns={"sdi":"int_sdi_id",
                "domNo": "int_domNo", "gi": "int_gi", "pdbId": "int_pdb",
                "chnLett": "int_chain", "from": "int_sdi_from", "to": "int_sdi_to",
                "num_residues": "int_nres", "whole_chain": "int_whole_chain"})
            st_domain_int = st_domain_int.drop(columns=["int_sdi_from", "int_sdi_to"])
            st_domain_int = st_domain_int.drop_duplicates()
            row = pd.merge(row, st_domain_int, how="left", on=["int_sdi_id", "int_gi"])
            del st_domain_int
            RealtimeLogger.info("MOL_INT ROW {}".format(row.iloc[0]))
        except TypeError as e:
            #SDI's doesn't exists, must be obsolete
            RealtimeLogger.info("Done row {}".format(e))
            return

        updated_resi = {"mol_res":[], "int_res": []}

        for resi in row.itertuples():
            try:
                mres = decode_residues(job, resi.mol_gi, resi.mol_pdb, resi.mol_chain, resi.mol_res, work_dir, resi)
                updated_resi["mol_res"].append(mres)
            except InvalidSIFTS as e:
                RealtimeLogger.info("Failed converting mol res names {} {}".format(updated_resi, e))
                updated_resi["mol_res"].append("")

            try:
                ires = decode_residues(job, resi.int_gi, resi.int_pdb, resi.int_chain, resi.int_res, work_dir, resi)
                updated_resi["int_res"].append(ires)
            except InvalidSIFTS as e:
                #This Row has failed converting binary to string skip it
                RealtimeLogger.info("Failed converting int res names {} {}".format(updated_resi, e))
                updated_resi["int_res"].append("")

        updated_resi = pd.DataFrame(updated_resi)

        row = row.drop(columns=["mol_res", "int_res"])

        if len(updated_resi.dropna()) > 0:
            row = pd.concat((row, updated_resi), axis=1)
            row.dropna()
        else:
            #This entire interation failed; returned None
            RealtimeLogger.info("Failed converting res names {}".format(updated_resi))
            return None

        RealtimeLogger.info("ROW {}".format(row.iloc[0]))

        row["mol_res"] = row["mol_res"].fillna("").astype(str)
        row["int_res"] = row["int_res"].fillna("").astype(str)
        row["int_domNo"] = row["int_domNo"].fillna(0).astype(int)
        row["int_nres"] = row["int_nres"].fillna(0).astype(int)
        row["int_sdi_id"] = row["int_sdi_id"].fillna(-1).astype(int)
        row["int_whole_chain"] = row["int_whole_chain"].fillna(0).astype(int)
        row["mol_sdi_id"] = row["mol_sdi_id"].fillna(-1).astype(int)

        RealtimeLogger.info("ROW {}".format(row.iloc[0]))

        path = "{}.h5".format(int_id)
        row.to_hdf(path, "table", table=True, format="table", complib="bzip2", complevel=9, min_itemsize=1024)
        out_store.write_output_file(path, "{}/_obsrows/{}".format(sfam_id, path))
        RealtimeLogger.info("Done row {}".format(int_id))

    except (SystemExit, KeyboardInterrupt):
        for f in (mmdb_path, ibis_path):
            try:
                os.remove(f)
            except (OSError, FileNotFoundError):
                pass
        raise
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        RealtimeLogger.info("FAILED {} {}".format(int_id, e, tb))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            f.write(str(e))
            f.write(str(tb))
        out_store.write_output_file(fail_file, "{}/_obsrows/{}.failed".format(sfam_id, int_id))
        try:
            os.remove(fail_file)
        except OSError:
            pass

    for f in (mmdb_path, ibis_path, path):
        try:
            os.remove(f)
        except (OSError, FileNotFoundError):
            pass


def merge_interactome_rows(job, sfam_id, pdbFileStoreID, ibisObsFileStoreID, update=False, cores=4, preemptable=True): #8
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")

    sfam_id = int(sfam_id) if sfam_id is not None else sfam_id

    status =  "observed"
    new_cols = ["mol_res", "int_res"]
    resi_prefix = "{}.observed_interactome".format(sfam_id)
    data_cols = ["obs_int_id", "mol_sdi", "int_sdi"]

    resi_path = os.path.join(work_dir, resi_prefix)

    store.download_input_directory("{}/_obsrows/".format(sfam_id), work_dir,
        postfix=".h5")

    import dask
    import dask.dataframe as dd
    from multiprocessing.pool import Pool

    dask.config.set(scheduler="processes")
    dask.config.set(pool=Pool(cores-1))

    import glob
    #
    # RealtimeLogger.info("WORK_DIR {} {}".format(os.path.join(work_dir),
    #     os.listdir(os.path.join(work_dir))))
    #
    all_rows = glob.glob(str(os.path.join(work_dir, "*.h5")))
    ddf = dd.read_hdf(all_rows, "table")
    ddf = ddf.repartition(npartitions=cores-1)

    ddf["int_domNo"] = ddf["int_domNo"].fillna(0).astype(int)
    ddf["int_nres"] = ddf["int_nres"].fillna(0).astype(int)
    ddf["int_sdi_id"] = ddf["int_sdi_id"].fillna(-1).astype(int)
    ddf["int_whole_chain"] = ddf["int_whole_chain"].fillna(0).astype(int)
    ddf["mol_sdi_id"] = ddf["mol_sdi_id"].fillna(-1).astype(int)


    RealtimeLogger.info("Starting merge for {}".format(sfam_id))
    try:
        os.remove(resi_path)
    except (OSError, FileNotFoundError):
        pass
    try:
        ddf.to_hdf(resi_path, "table", format="table", table=True, complevel=9,
            data_columns=data_cols, complib="bzip2", min_itemsize=1024)
    except ValueError:
        try:
            os.remove(resi_path)
        except (OSError, FileNotFoundError):
            pass
        old_df = None
        old_dtypes = None
        for i, f in enumerate(all_rows):
            df = pd.read_hdf(f, "table", mode="r")
            RealtimeLogger.info("DF is {}".format(df.iloc[0]))
            df["int_domNo"] = df["int_domNo"].fillna(0).astype(int)
            df["int_nres"] = df["int_nres"].fillna(0).astype(int)
            df["int_sdi_id"] = df["int_sdi_id"].fillna(-1).astype(int)
            df["int_whole_chain"] = df["int_whole_chain"].fillna(0).astype(int)
            df["mol_sdi_id"] = df["mol_sdi_id"].fillna(-1).astype(int)
            RealtimeLogger.info("DF2 is {}".format(df.iloc[0]))
            try:
                if i==0:
                    df.to_hdf(resi_path, "table", format="table", table=True,
                        data_columns=data_cols, complevel=9, complib="bzip2",
                        min_itemsize=1024)
                else:
                    df.to_hdf(resi_path, "table", mode="a", append=True,
                        format="table", table=True, data_columns=data_cols,
                        complevel=9, complib="bzip2", min_itemsize=1024)
            except ValueError as e:
                import traceback
                tb = traceback.format_exc()
                RealtimeLogger.info("FAILED merging rows for {} {}".format(sfam_id, tb))
                RealtimeLogger.info("PREV COLS {}".format(old_dtypes))
                RealtimeLogger.info("NEW  COLS {}".format(df[["int_domNo","int_nres","int_sdi_id","int_whole_chain","mol_sdi_id"]].dtypes))
                RealtimeLogger.info("PREV failed rows {}".format(
                    old_df[["int_domNo","int_nres","int_sdi_id","int_whole_chain","mol_sdi_id"]]))
                RealtimeLogger.info("NEW failed rows {}".format(
                    df[["int_domNo","int_nres","int_sdi_id","int_whole_chain","mol_sdi_id"]]))
                RealtimeLogger.info("CONCAT {}".format(pd.concat((old_df,
                    df[["int_domNo","int_nres","int_sdi_id","int_whole_chain","mol_sdi_id"]]
                ), axis=0)))
                raise
            old_df = df[["int_domNo","int_nres","int_sdi_id","int_whole_chain","mol_sdi_id"]]
            old_dtypes = df[["int_domNo","int_nres","int_sdi_id","int_whole_chain","mol_sdi_id"]].dtypes
            #del df



    RealtimeLogger.info("Finished merge for {}".format(sfam_id))
    #Combine residues into dataframe
    # possible_errors = []
    # nrows = None
    # for nrows, row_prefix in enumerate(out_store.list_input_directory("{}/_obsrows/".format(sfam_id))):
    #     if row_prefix.endswith("failed"): continue
    #     job.log("Running {} {}".format(sfam_id, row_prefix))
    #
    #     row_file = os.path.join(work_dir, os.path.basename(row_prefix))
    #     out_store.read_input_file(row_prefix, row_file)
    #
    #     df = pd.read_hdf(row_file, "table")
    #     try:
    #         df.to_hdf(str(resi_path), "table", table=True, format="table", append=True, mode="a",
    #             data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except:
    #         import traceback
    #         tb = traceback.format_exc()
    #         job.log("Failed writing {}: {} {}".format(sfam_id, resi_path, tb))
    #         possible_errors.append(tb)
    #         continue
    #
    #     try:
    #         os.remove(row_file)
    #     except OSError:
    #         pass
    #
    #     out_store.remove_file(row_prefix)

    files_to_remove = all_rows

    if os.path.isfile(resi_path):
        #Upload to S3
        store.write_output_file(resi_path, os.path.join(str(sfam_id), resi_prefix))
        files_to_remove.append(resi_path)

        if update:
            done_file = os.path.join(work_dir, "DONE")
            with open(done_file, "w") as f:
                pass
            store.write_output_file(done_file, "{sfam_id}/DONE".format(sfam_id=str(sfam_id)))
            files_to_remove.append(done_file)

    else:
        RealtimeLogger.info("Failed merging: {}".format(resi_path))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            pass

        store.write_output_file(fail_file, "{}/{}.failed".format(sfam_id, resi_prefix))

        files_to_remove.append(fail_file)

    for f in files_to_remove:
        try:
            os.remove(f)
        except (OSError, FileNotFoundError):
            pass

    #job.addFollowOnJobFn(cleanup_sfam, sfam_id, pdbFileStoreID, ibisObsFileStoreID)

def get_observed_structural_interactome(job, sfam_id, pdbFileStoreID, ibisObsFileStoreID, further_parallelize=True, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")

    ibis_obs_path = get_file(job, "IBIS_obs.h5", ibisObsFileStoreID)

    try:
        if sfam_id is not None:
            df = filter_hdf(ibis_obs_path, "ObsInt", "mol_superfam_id", float(sfam_id), columns=["obs_int_id"])
            sfam_id = int(sfam_id)
        else:
            df = pd.read_hdf(ibis_obs_path, "ObsInt", columns=["obs_int_id", "mol_superfam_id"])
            df = df[df["mol_superfam_id"]==-1][["obs_int_id"]]
        int_ids = df["obs_int_id"].drop_duplicates().dropna().astype(int)
        if len(int_ids) == 0:
            RealtimeLogger.info("EMPTY OBS SFAM {}".format(sfam_id))
            return
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        RealtimeLogger.info("FAILED OBS SFAM {} {}".format(sfam_id, e))
        return

    RealtimeLogger.info("IDs {}".format(list(int_ids)))
    current_rows = set(int(os.path.basename(key)[:-3]) for key in out_store.list_input_directory("{}/_obsrows".format(sfam_id)) if not key.endswith("failed"))
    int_ids = list(set(int_ids)-current_rows)
    RealtimeLogger.info("Will run {} ids: {}".format(len(int_ids), int_ids))

    # int_ids = [608369]

    if len(int_ids) > 0:
        if further_parallelize:
            #Add jobs for each interaction
            map_job(job, process_observed_interaction, int_ids, sfam_id, ibisObsFileStoreID, pdbFileStoreID)
        else:
            for int_id in int_ids:
                process_observed_interaction(job, int_id, sfam_id, ibisObsFileStoreID, pdbFileStoreID)

    #Merge converted residues
    job.addFollowOnJobFn(merge_interactome_rows, sfam_id, pdbFileStoreID, ibisObsFileStoreID)

    try:
        os.remove(ibis_obs_path)
    except (OSError, FileNotFoundError):
        pass

def cleanup(job, preemptable=True):
    in_store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")
    keys = list(in_store.list_input_directory())
    finished = set(key.split("/")[0] for key in keys if key.endswith("observed_interactome"))
    failed = set(key.split("/")[0] for key in keys if "_obsrows" in key and key.endswith("failed"))

    for key in failed-finished:
        in_store.remove_file(key)

def cleanup_sfam(job, sfam_id, pdbFileStoreID, ibisObsFileStoreID, further_parallelize=True, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")

    sfam_id = int(sfam_id) if sfam_id is not None else sfam_id

    key = "{sfam_id}/{sfam_id}.observed_interactome".format(sfam_id=sfam_id)
    obs_int_file = "{sfam_id}.observed_interactome".format(sfam_id=sfam_id)

    try:
        store.read_input_file(key, obs_int_file)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        RealtimeLogger.info("Could not get observed_interactome file for {}".format(sfam_id))
        return

    try:
        obs_int = pd.read_hdf(str(obs_int_file), "table", mode="r")
    except (IOError, KeyError):
        obs_int = None

    if obs_int is not None and len(obs_int)>0:
        try:
            obs_int = obs_int.drop(columns=["mol_sdi_from", "mol_sdi_to", "int_sdi_from", "int_sdi_to"])
        except KeyError:
            pass

        obs_int = obs_int.drop_duplicates()

        obs_int = obs_int[obs_int["mol_superfam_id"]==sfam_id]

        data_cols = ["obs_int_id", "mol_sdi", "int_sdi"]
        obs_int.to_hdf(str(obs_int_file)+".fixed", "table", table=True, format="table",
            data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)

        new_key = "{sfam_id}/_obsrows/to_update.h5".format(sfam_id=sfam_id)
        store.write_output_file(str(obs_int_file)+".fixed", new_key)

        skip_ids_file = set(obs_int["obs_int_id"].astype(int))
    else:
        skip_ids_file = set()

    def only_int(s):
        try:
            return int(s)
        except ValueError:
            return None

    done_rows = set(only_int(os.path.splitext(os.path.basename(k))[0]) for k in store.list_input_directory(
        "{sfam_id}/_obsrows/".format(sfam_id=sfam_id)) if k.endswith("h5") and \
        not k.endswith("to_update.h5"))

    skip_ids = skip_ids_file.union(done_rows)

    RealtimeLogger.info("RUNNING sfam {}: {}+{}={}".format(sfam_id, len(skip_ids_file),
        len(done_rows), len(skip_ids)))
    #RealtimeLogger.info("SKIP ROWS: {}".format(skip_ids))

    ibis_path = get_file(job, "IBIS_obs.h5", ibisObsFileStoreID)
    all_int_ids = set(filter_hdf(ibis_path, "ObsInt", "mol_superfam_id", sfam_id,
        columns=["obs_int_id"])["obs_int_id"])
    rerun_int_ids = list(all_int_ids-skip_ids)

    RealtimeLogger.info("RUNNING sfam {}: {}/{} interactions".format(sfam_id,
        len(rerun_int_ids), len(all_int_ids)))
    RealtimeLogger.info("RUNNING ROWS: {}".format(rerun_int_ids))

    if len(rerun_int_ids) > 0:
        if further_parallelize:
            map_job(job, process_observed_interaction, rerun_int_ids, sfam_id, ibisObsFileStoreID, pdbFileStoreID)
        else:
            for _id in rerun_int_ids:
                try:
                    process_observed_interaction(job, _id, sfam_id, ibisObsFileStoreID, pdbFileStoreID)
                except (SystemExit, KeyboardInterrupt):
                    raise
                except:
                    pass
    else:
        done_file = os.path.join(work_dir, "DONE")
        with open(done_file, "w") as f:
            pass
        store.write_output_file(done_file, "{sfam_id}/DONE".format(sfam_id=sfam_id))
        try:
            os.remove(done_file)
        except OSError:
            pass

    if len(skip_ids_file)+len(done_rows)>=len(all_int_ids):
        job.addFollowOnJobFn(merge_interactome_rows, sfam_id, pdbFileStoreID,
            ibisObsFileStoreID, True)

    for f in (obs_int_file, str(obs_int_file)+".fixed", "IBIS_obs.h5", 'DONE'):
        try:
            os.remove(f)
        except OSError:
            pass

def merge_full_interactome(job, cores=8, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")

    resi_prefix = "observed_interactome.h5"
    data_cols = ["obs_int_id", "mol_sdi", "int_sdi"]

    resi_path = os.path.join(work_dir, resi_prefix)

    store.download_input_directory("", work_dir, postfix=".observed_interactome")

    import dask
    import dask.dataframe as dd
    from multiprocessing.pool import Pool

    dask.config.set(scheduler="processes")
    dask.config.set(pool=Pool(cores-1))

    import glob
    #
    # RealtimeLogger.info("WORK_DIR {} {}".format(os.path.join(work_dir),
    #     os.listdir(os.path.join(work_dir))))
    #
    all_sfams = glob.glob(str(os.path.join(work_dir, "*.observed_interactome")))
    ddf = dd.read_hdf(all_sfams, "table")
    ddf = ddf.repartition(npartitions=cores-1)
    ddf.to_hdf(resi_path, "table", format="table", table=True, complevel=9,
        data_columns=data_cols, complib="bzip2", min_itemsize=1024)

    if os.path.isfile(resi_path):
        #Upload to S3
        out_store.write_output_file(resi_path, resi_prefix)

        #Cleanup
        for f in all_sfams+[resi_path]:
            try:
                os.remove(resi_path)
            except (OSError, FileNotFoundError):
                pass

    else:
        RealtimeLogger.info("Failed merging: {}".format(resi_path))
        fail_file = os.path.join(work_dir, "fail_file")
        with open(fail_file, "w") as f:
            pass
        out_store.write_output_file(fail_file, "full_merge.failed")
        for f in all_rows+[fail_file]:
            try:
                os.remove(resi_path)
            except (OSError, FileNotFoundError):
                pass

def start_toil(job, update=False, pdbFileStoreID=None, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    out_store = IOStore.get("aws:us-east-1:molmimic-interfaces-2")

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

    #Read in superfamilies
    #pdb = filter_hdf_chunks(str(ibis_obs_path), "ObsInt", columns=["mol_superfam_id"]).drop_duplicates()

    if not update:
        #Choose which superfamilies to run, skip those already present
        skip_sfam = set([float(f.split("/", 1)[0]) for f in out_store.list_input_directory() \
            if f.endswith(".observed_interactome")])

        sfams = pdb[~pdb["mol_superfam_id"].isin(skip_sfam)]["mol_superfam_id"].drop_duplicates().dropna().astype(int)
        #sfams.append(None) #Process those with no superfamily info
        RealtimeLogger.info("Will run a total of {} SFAMS: {}".format(len(sfams), sfams))

        #Run all superfamilies
        map_job(job, get_observed_structural_interactome, sfams, pdbFileStoreID, ibisObsFileStoreID)
    else:
        RealtimeLogger.info("Running UPDATE")
        sfams = [279669] #pdb["mol_superfam_id"]
        map_job(job, cleanup_sfam, sfams, pdbFileStoreID, ibisObsFileStoreID)

    #Cleanup
    job.addFollowOnJobFn(cleanup)
    job.addFollowOnJobFn(merge_full_interactome)
    os.remove(ibis_obs_path)
    os.remove(pdb_path)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    os.environ['TOIL_CUSTOM_DOCKER_INIT_COMMAND'] ='pip install awscli && eval $(aws ecr get-login --no-include-email --region us-east-1)'

    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--update", default=False, action="store_true")
    options = parser.parse_args()
    #options.logLevel = "OFF"
    options.clean = "always"
    options.targetTime = 1


    job = Job.wrapJobFn(start_toil, options.update)
    with Toil(options) as workflow:
        workflow.start(job)
