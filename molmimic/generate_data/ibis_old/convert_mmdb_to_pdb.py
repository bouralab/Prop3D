import os, sys
import shutil
import string
import logging
from itertools import product
#
import numpy as np
import pandas as pd
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

from molmimic.generate_data.map_residues import mmdb_to_pdb_resi, mmdb_to_pdb_resi_obsolete, InvalidSIFTS, ChangePDB
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.util import get_file

from botocore.exceptions import ClientError

from toil.realtimeLogger import RealtimeLogger

def convert_row(job, row, work_dir=None):
    resi = [str(row["from"]), str(row["to"])]
    # try:
    #     frm, to = list(mmdb_to_pdb_resi(row["pdbId"], row["chnLett"], resi, replace_nulls=True, job=job))
    # except (IOError, TypeError, ValueError, InvalidSIFTS) as e:
    try:
        frm, to, num_residues = list(mmdb_to_pdb_resi_obsolete(row["gi"], row["pdbId"],
            row["chnLett"], resi, shrink_or_expand=True, return_gi_len=True,
            replace_nulls=False, job=job))
    except ChangePDB as e:
        if None in e.new_res:
            return pd.Series({
                "pdbId": row["pdbId"],
                "chnLett": row["chnLett"],
                "num_residues": np.NaN,
                "whole_chain": np.NaN,
                "from": np.NaN,
                "to":   np.NaN})
        else:
            whole_chain = e.num_residues-row["to"]+row["from"]<12
            return pd.Series({
                "pdbId": e.new_pdb,
                "chnLett": e.new_chain,
                "num_residues": e.num_residues,
                "whole_chain": whole_chain,
                "from": "".join(map(str, e.new_res[0][1:])).strip(),
                "to":   "".join(map(str, e.new_res[1][1:])).strip()})
    except (SystemExit, KeyboardInterrupt):
        raise
    except InvalidSIFTS as e:
        import traceback
        RealtimeLogger.info("Error mapping mmdb for {} '{}' {} {} {} {} {}".format(
            row["pdbId"], row["chnLett"], row["sdi"], resi, e.__class__.__name__, e, traceback.format_exc()))
        return pd.Series({
            "pdbId": row["pdbId"],
            "chnLett": row["chnLett"],
            "num_residues": np.NaN,
            "whole_chain": np.NaN,
            "from": np.NaN,
            "to":   np.NaN})
    except Exception as e:
        raise
        import traceback
        RealtimeLogger.info("Error mapping mmdb for {} '{}' {} {} {} {} {}".format(
            row["pdbId"], row["chnLett"], row["sdi"], resi, e.__class__.__name__, e, traceback.format_exc()))
        return pd.Series({
            "pdbId": row["pdbId"],
            "chnLett": row["chnLett"],
            "num_residues": np.NaN,
            "whole_chain": np.NaN,
            "from": np.NaN,
            "to":   np.NaN})

    if None in (frm, to):
        return pd.Series({
            "pdbId": row["pdbId"],
            "chnLett": row["chnLett"],
            "num_residues": np.NaN,
            "whole_chain": np.NaN,
            "from": np.NaN,
            "to":   np.NaN})

    whole_chain = num_residues-row["to"]+row["from"]<12
    return pd.Series({
        "pdbId": row["pdbId"],
        "chnLett": row["chnLett"],
        "num_residues": num_residues,
        "whole_chain": whole_chain,
        "from": "".join(map(str, frm[1:])).strip(),
        "to":   "".join(map(str, to[1:])).strip()})

def mmdb2pdb(job, pdb_group, mmdb_file, missing=False, use_sdi=False, use_chunk=True, dask=False, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()

    store = IOStore.get("aws:us-east-1:molmimic-ibis")
    mmdb_path = os.path.join(work_dir, "MMDB.h5")

    if mmdb_file is None:
        store.read_input_file("MMDB.h5", mmdb_path)
    else:
        get_file(job, mmdb_path, mmdb_file, work_dir=work_dir)

    if pdb_group is not None:
        tmp_path = os.path.join(work_dir, "{}{}.h5".format("sdi" if missing or use_sdi else "", pdb_group))
        if not missing and not use_sdi and not use_chunk:
            sdoms = pd.read_hdf(str(mmdb_path), "StrucutralDomains")
            if len(pdb_group) == 2:
                sdoms = sdoms[sdoms["pdbId"].str[:3]==pdb_group.upper()]
            elif len(pdb_group) == 3:
                sdoms = sdoms[sdoms["pdbId"].str.startswith(pdb_group)]
            elif len(pdb_group) == 4:
                sdoms = sdoms[sdoms["pdbId"]==pdb_group]
        elif use_sdi:
            sdoms = pd.read_hdf(str(mmdb_path), "StrucutralDomains")
            sdoms = sdoms[sdoms["sdi"]==int(pdb_group)]
        elif use_chunk:
            for i, sdoms in enumerate(pd.read_hdf(str(mmdb_path), "StrucutralDomains", chunksize=100)):
                if i == pdb_group:
                    break
            else:
                return
    else:
        tmp_path = os.path.join(work_dir, "PDB.h5")

    #sdoms = sdoms.assign(length=sdoms["to"]-sdoms["from"])
    #sdoms = sdoms[sdoms["length"]>4]
    #del sdoms["length"]

    # sdoms_l = sdoms.assign(length=sdoms["to"]-sdoms["from"])
    # short_sdi = sdoms_l[sdoms_l["length"]<15]["sdi"].drop_duplicates()
    # total_lens = sdoms_l[sdoms_l["sdi"].isin(short_sdi)].groupby("sdi")["length"].sum()
    # bad_sdi = total_lens[total_lens<15].index
    # sdoms = sdoms[~sdoms.isin(bad_sdi)]
    # del sdoms_l

    if dask:
        ddf = dd.from_pandas(sdoms, npartitions=NUM_WORKERS)
        meta = pd.DataFrame({"from":[str], "to":[str]})

        with ProgressBar():
            pdb_map = ddf.map_partitions(lambda _df: _df.apply(\
                lambda row: convert_row(job, row), axis=1), meta=meta)\
                    .compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    else:
        pdb_map = sdoms.apply(lambda row: convert_row(job, row, work_dir=work_dir), axis=1)

    sdoms.loc[:, ["pdbId", "chnLett", "from", "to"]] = pdb_map #sdoms = pd.concat((sdoms, pdb_map), axis=1) #
    sdoms = sdoms.assign(
        num_residues = pdb_map["num_residues"],
        whole_chain = pdb_map["whole_chain"])

    sdoms = sdoms.dropna()
    # sdoms["from"] = sdoms["from"].astype(str)
    # sdoms["to"] = sdoms["to"].astype(str)
    sdoms[["pdbId", "chnLett", "from", "to"]] = sdoms[["pdbId", "chnLett", "from", "to"]].astype(str)
    sdoms[["sdi", "domNo", "gi", "num_residues", "whole_chain"]] = \
        sdoms[["sdi", "domNo", "gi", "num_residues", "whole_chain"]].astype(int)


    if len(sdoms)>0:
        sdoms.to_hdf(str(tmp_path), "table", table=True, format="table", complevel=9, complib="bzip2", min_itemsize=756)
        store.write_output_file(tmp_path, "{}/{}".format("sdi" if missing else "chunks_", os.path.basename(tmp_path)))
    else:
        RealtimeLogger.info("Failed group {}".format(pdb_group))

# def dl(key, output_dir):
#     store_ = IOStore.get("aws:us-east-1:molmimic-ibis")
#     outfile = os.path.join(output_dir, os.path.basename(key))
#     store_.read_input_file(key, outfile)
#     del store_

def merge(job, mmdb_file, missing=False, cores=8, preemptable=False):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-ibis")

    merged_file = os.path.join(work_dir, "PDB.h5")
    prefix = "chunks_/"

    if missing:
        store.read_input_file("PDB.h5", merged_file)
        prefix = "sdi/"
        pdbStore = pd.HDFStore("PDB.h5")
        for k in ("StructuralDomains", "resolu"):
            try:
                del pdbStore[k]
            except KeyError:
                pass
        pdbStore.close()

    # from joblib import Parallel, delayed
    #
    # Parallel(n_jobs=-1)(delayed(dl)(k, work_dir) for k in store.list_input_directory(prefix))

    # for i, pdb_group_key in enumerate(store.list_input_directory(prefix)):
    #     if i >= 10: break
    #     pdb_group_file = os.path.join(work_dir, os.path.basename(pdb_group_key))
    #     store.read_input_file(pdb_group_key, pdb_group_file)

    store.download_input_directory(prefix, work_dir)

    output = os.path.join(work_dir, "PDB.h5")
    import dask
    from multiprocessing.pool import Pool

    dask.config.set(scheduler="processes")
    dask.config.set(pool=Pool(cores-1))

    import glob
    #
    # RealtimeLogger.info("WORK_DIR {} {}".format(os.path.join(work_dir),
    #     os.listdir(os.path.join(work_dir))))
    #
    ddf = dd.read_hdf(glob.glob(str(os.path.join(work_dir, "*.h5"))), "table")
    ddf = ddf.repartition(npartitions=cores-1)
    ddf.to_hdf(output, "StructuralDomains", format="table", table=True, complevel=9,
        complib="bzip2", min_itemsize=768)
    # # remove_keys = []
    # for pdb_group_key in store.list_input_directory(prefix):
    #     pdb_group_file = os.path.join(work_dir, os.path.basename(pdb_group_key))
    #     store.read_input_file(pdb_group_key, pdb_group_file)
    #
    #     try:
    #         df = pd.read_hdf(str(pdb_group_file), "table")
    #     except IOError as e:
    #         RealtimeLogger.info("Stuck {}".format(e))
    #         continue
    #
    #     RealtimeLogger.info("{} DF {} {}".format(pdb_group_key, df.head(), df.shape))
    #
    #     df[["sdi", "domNo", "gi"]].astype(int)
    #     df[['pdbId', 'chnLett', 'from', 'to']] = df[['pdbId', 'chnLett', 'from', 'to']].astype(str)
    #
    #     try:
    #         df.to_hdf(str(merged_file), "StructuralDomains", format="table", table=True,
    #             mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)
    #     except ValueError:
    #         pass
    #
    #     del df
    #
    #     try:
    #         os.remove(pdb_group_file)
    #     except OSError:
    #         pass

        #remove_keys.append(pdb_group_key)

    mmdb_path = os.path.join(work_dir, "MMDB.h5")
    if mmdb_file is None:
        store.read_input_file("MMDB.h5", mmdb_path)
    else:
        get_file(job, mmdb_path, mmdb_file, work_dir=work_dir)

    sfams = pd.read_hdf(str(mmdb_path), "Superfamilies", mode="r")
    sfams = sfams[['sdsf_id', 'sdi', 'sfam_id', 'label', 'whole_chn', 'mmdb_id', 'mol_id', 'pssm']]
    sfams.to_hdf(str(merged_file), "Superfamilies", format="table", table=True,
        mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)

    sdoms = pd.read_hdf(str(merged_file), "StructuralDomains")
    merged = pd.merge(sdoms, sfams, how="left", on="sdi").dropna()
    merged.to_hdf(str(merged_file), "StructuralDomainFamilies", format="table",
        table=True, mode='a', append=True, complevel=9, complib="bzip2",
        min_itemsize=768)

    resolu = pd.read_table("ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx",
        header=None, names=["pdbId", "resolution"], skiprows=6, sep="\t;\t")
    resolu.to_hdf(str(merged_file), "resolu", format="table", table=True,
        mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)

    store.write_output_file(merged_file, "PDB.h5")

    # for key in remove_keys:
    #     pass #store.remove_file(key)

    try:
        os.remove(merged_file)
    except OSError:
        pass


def start_toil(job, missing=False, use_sdi=False, use_chunk=True, force_merge=True, name="mmdb2pdb"): #, preemptable=True):
    """Initiate Toil Jobs for converting MMDB to PDB

    Paramters
    ---------
    job : toil.Job

    Return
    ------
    mmdb2pdb : toil.Job
    """
    work_dir = job.fileStore.getLocalTempDir()

    store = IOStore.get("aws:us-east-1:molmimic-ibis")

    # if jobStore.exists("PDB.h5"):
    #     return

    try:
        mmdb_file = os.path.join(work_dir, "MMDB.h5")
        store.read_input_file("MMDB.h5", mmdb_file)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        raise RuntimeError("download_data must be run before mmdb2pdb")

    mmdbFileStoreID = job.fileStore.writeGlobalFile(mmdb_file)

    if not force_merge:

        mmdb = pd.read_hdf(mmdb_file, "StrucutralDomains", mode="r")
        mmdb["sdi"] = mmdb["sdi"].astype(int)
        #mmdb = mmdb.assign(length=mmdb["to"]-mmdb["from"])
        # groups = list(range(int(mmdb.shape[0]/100)+1))
        # use_chunk = True


        if missing and store.exists("PDB.h5"):
            mmdb["sdi"] = mmdb["sdi"].astype(int)
            mmdb = mmdb.assign(length=mmdb["to"]-mmdb["from"])
            short_sdi = mmdb[mmdb["length"]<15]["sdi"].drop_duplicates()
            total_lens = mmdb[mmdb["sdi"].isin(short_sdi)].groupby("sdi")["length"].sum()
            bad_sdi = total_lens[total_lens<15].index
            mmdb = mmdb[~mmdb.isin(bad_sdi)].dropna()

            pdb_file = os.path.join(work_dir, "PDB.h5")
            RealtimeLogger.info("Download PDB file")
            store.read_input_file("PDB.h5", pdb_file)
            RealtimeLogger.info("Done Download PDB file")

            groups = set(mmdb["sdi"])
            del mmdb
            for pdb_sdi_group in pd.read_hdf(pdb_file, "StructuralDomains", columns=["sdi"], mode="r", chunksize=1000):
                groups -= set(pdb_sdi_group["sdi"].astype(int))

            #existing = set(int(k[7:-3]) for k in store.list_input_directory("sdi/"))
            #groups -= existing

            #range(int(mmdb.shape[0]/100)+1)

            groups = list(groups)[:1]
            use_sdi = True
            RealtimeLogger.info("USING SDIs Instead: {}".format(len(groups)))
        elif use_chunk:
            groups = list(range(int(mmdb.shape[0]/100)+1))
            use_chunk = True
            # skip_chunks = set(int(f[8:-3]) for f in store.list_input_directory("chunks_"))
            # RealtimeLogger.info("SKIP_CHUNKS={} {}".format(len(skip_chunks), list(skip_chunks)[:5]))
            # if len(skip_chunks) > 0:
            #     groups = list(set(groups)-skip_chunks)
            #     RealtimeLogger.info("RUN GRoup {}".format(len(groups)))

        else:
            existing = pd.Series([k[7:-3] for k in store.list_input_directory("groups/")]).astype(str)

            if use_sdi:
                groups = mmdb[~mmdb["sdi"].isin(existing[existing.str.startswith("sdi")].str[3:].astype(int))]["sdi"].drop_duplicates()
                RealtimeLogger.info("USING SDIs Instead: {}".format(len(groups)))
            else:
                groups = mmdb[(~mmdb["pdbId"].str[:3].isin(existing[existing.str.len()==3]))&\
                              (mmdb["length"]>4)]["pdbId"].str[:3].drop_duplicates()

                RealtimeLogger.info("TOTAL GROUPS: {}".format(len(groups)))

                if 10 <= groups.shape[0] < 30:
                    groups = mmdb[(~mmdb["pdbId"].str[:3].isin(existing[existing.str.len()==3]))&\
                                  (mmdb["length"]>4)]["pdbId"].drop_duplicates()
                    groups = groups[~groups.isin(existing[existing.str.len()==4])]
                    RealtimeLogger.info("USING PDBs Instead: {}".format(len(groups)))
                elif groups.shape[0] < 10:
                    groups = mmdb[(~mmdb["pdbId"].str[:3].isin(existing[existing.str.len()==3]))&\
                                  (mmdb["length"]>4)]
                    groups = groups[~groups["sdi"].isin(existing[existing.str.startswith("sdi")].str[3:].astype(int))]["sdi"].drop_duplicates()
                    use_sdi = True
                    RealtimeLogger.info("USING SDIs Instead: {}".format(len(groups)))

        map_job(job, mmdb2pdb, groups, mmdbFileStoreID, missing, use_sdi, use_chunk)

    job.addFollowOnJobFn(merge, mmdbFileStoreID, missing)

    try:
        os.remove(mmdb_file)
    except OSError:
        pass

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    os.environ["TOIL_CUSTOM_DOCKER_INIT_COMMAND"] = 'pip install awscli && eval $(aws ecr get-login --no-include-email --region us-east-1)'

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    #options.targetTime = 1

    #detector = logging.StreamHandler()
    #logging.getLogger().addHandler(detector)

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as workflow:
        from boto.utils import get_instance_metadata
        instanceMetadata = get_instance_metadata()["iam"]["security-credentials"]["toil_cluster_toil"]
        RealtimeLogger.info("CREDS={}".format(instanceMetadata))
        workflow.start(job)
