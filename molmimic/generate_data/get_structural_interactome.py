import os, sys
import glob
from multiprocessing.pool import ThreadPool

import binascii
import zlib
from pyasn1.codec.ber import decoder
import shutil

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
from molmimic.generate_data.map_residues import mmdb_to_pdb_resi

dask.config.set(scheduler='multiprocessing', num_workers=4)
dask.config.set(pool=ThreadPool(4))

def decode_residues(job, pdb, chain, res, row):
    residues = []

    if res.startswith("0x"):
        res = res[2:]
    try:
        res = binascii.unhexlify(res)
    except:
        pass

    try:
        code, rest = decoder.decode(zlib.decompress(res, 16 + zlib.MAX_WBITS))
    except Exception as e:
        if type(res, str) and "," in res:
            return res
        else:
            return np.NaN

    for i in xrange(len(code)):
        c = code[i]
        range_from, range_to, gi = tuple([c[j] for j in range(len(c))])
        for x in xrange(range_from, range_to + 1):
            residues.append(x)

    try:
        return ",".join(map(str, mmdb_to_pdb_resi(pdb, chain, residues, job=job)))
    except (IOError, TypeError) as error:
        print "Error mapping mmdb for", pdb, chain, error, row

        residues.insert(0, "mmdb")
        return ",".join(map(str,residues))

def convert_residues(job, row, only_mol=False):
    try:
        resi = {"mol_res": decode_residues(job, row["mol_pdb"], row["mol_chain"], row["mol_res"], row)}
        if not only_mol:
            resi["int_res"] = decode_residues(job, row["int_pdb"], row["int_chain"], row["int_res"], row)
        return pd.Series(resi)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        return None

def merge_interactome_resi(job, sfam_id, converted_residues, newIbisFileStoreID, inferred_table=None):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    if inferred_table is not None:
        status = "inferred (table {})".format(inferred_table)
        resi_prefix = "Intrac{}_{}.inferred_interactome".format(inferred_table, int(sfam_id))
        new_cols = ["mol_res"]
        data_cols = ["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"]
    else:
        status =  "observed"
        new_cols = ["mol_res", "int_res"]
        resi_prefix = "{}.observed_interactome".format(int(sfam_id))
        data_cols = ["obs_int_id", "mol_sdi", "int_sdi"]

    #Combine residues into dataframe
    conv_resi = pd.concat(converted_residues)

    if conv_resi.shape[0] == 0:
        job.log("Failed reading converted {} residues from {}".format(status, sfam_id))

    #Downlaod IBIS obs ints
    ibis_path = get_file(newIbisFileStoreID)

    #Merge resdides with inferred interactrome
    ints = pd.read_hdf(unicode(ibis_path), "table")
    ints.loc[:, new_cols] = conv_resi

    #Save to file
    resi_path = os.path.join(work_dir, resi_path)
    job.log("Writing {} interactome: {}".format(sfam_id, resi_path))
    try:
        obs_ints.to_hdf(unicode(path), "table", table=True, format="table",
            data_columns=data_cols, complib="bzip2", complevel=9, min_itemsize=1024)
    except TypeError:
        job.log("Failed writing {}: {}".format(sfam_id, resi_path))
        raise

    #Upload to S3
    out_store.write_output_file(resi_path, os.path.join(str(int(sfam_id)), resi_prefix))

    #Cleanup
    job.fileStore.deleteGlobalFile(newIbisFileStoreID)
    os.remove(resi_path)

def get_observed_structural_interactome(job, sfam_id, pdbFileStoreID, ibisObsFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    out_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    ibis_obs_path = get_file(job, "IBIS_obs.h5", ibisObsFileStoreID)
    mmdb_path = get_file(job, "PDB.h5", out_store)

    #Read in PPIs and save from the given CDD
    obs_ints = pd.read_hdf(unicode(ibis_obs_path), "ObsInt", where="mol_superfam_id=={}".format(sfam_id))
    #obs_ints = obs_ints[obs_ints["mol_superfam_id"]==sfam_id]
    print obs_ints.shape
    #Read in face1 residues
    face1 = pd.read_hdf(unicode(ibis_obs_path), "MolResFace", where="obs_int_id=obs_ints['obs_int_id']")
    face1.columns = ["obs_int_id", "mol_res"]
    print face1.shape
    #Keep entries from current CDD
    obs_ints = pd.merge(obs_ints, face1, how="left", on="obs_int_id")
    del face1

    #Read in face2 residues and convert gzipped asn1 into res numbers
    face2 = pd.read_hdf(unicode(ibis_obs_path), "IntResFace", where="obs_int_id=obs_ints['obs_int_id']")
    face2.columns = ["obs_int_id", "int_res"]
    print face2.shape
    #Keep entries from current CDD
    obs_ints = pd.merge(obs_ints, face2, how="left", on="obs_int_id")
    del face2
    print obs_ints.shape

    st_domain = pd.read_hdf(unicode(mmdb_path), "StructuralDomains")
    st_domain.columns = ['mol_sdi', 'mol_domNo', 'mol_gi', 'mol_pdb', 'mol_chain', 'mol_sdi_from', 'mol_sdi_to']

    obs_ints = pd.merge(obs_ints, st_domain, how="left", left_on="mol_sdi_id", right_on="mol_sdi")

    st_domain.columns = ['int_sdi', 'int_domNo', 'int_gi', 'int_pdb', 'int_chain', 'int_sdi_from', 'int_sdi_to']
    obs_ints = pd.merge(obs_ints, st_domain, how="left", left_on="int_sdi_id", right_on="int_sdi")
    del st_domain

    #sfams = pd.read_hdf(unicode(mmdb_path), "Superfamilies")[["sdi", "label", "sfam_id"]].dropna(axis=0, subset=["sfam_id"])
    #sfams = sfams[sfams["sfam_id"]==sfam_id].rename(columns={"sdi":"mol_sdi", "label":"int_superfam_label", "sfam_id":"int_superfam_id"})
    #obs_ints = pd.merge(obs_ints, sfams, how="left", on="mol_sdi").dropna(axis=0, subset=["mol_sdi"])
    #del sfams

    new_obs_path = os.path.join(work_dir, "{}.observed_interactome.tmp".format(int(sfam_id)))
    try:
        obs_ints.to_hdf(unicode(new_obs_path), "table", table=True,
            data_columns=["obs_int_id", "mol_sdi", "int_sdi"], format="table", complevel=9, complib="bzip2", min_itemsize=1024)
    except TypeError:
        job.log("Failed writing {}: {}".format(sfam_id, new_obs_path))
        raise

    #Add ibis info into local job store
    newIbisObsFileStoreID = job.fileStore.writeGlobalFile(new_obs_path)

    #Add jobs for each interaction
    save_cols = ["mol_pdb", "mol_chain", "mol_res" "int_pdb", "int_chain", "int_res"]
    resi = [job.addChildJobFn(convert_residues).rv() for obs_int in \
        obs_ints[save_cols].itertuples()]

    #Merge converted residues
    job.addFollowOnJobFn(merge_interactome_resi, sfam_id, resi, pdbFileStoreID, \
        newIbisObsFileStoreID)

def process_inferred_sfam_table(job, mol_sfam_id, table, groupInfFileStoreID, pdbFileStoreID, cores=2):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))
    out_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    df_key = "{}_{}.h5".format(mol_sfam_id, table)
    tmp_df_file = get_file(job, "tmp"+df_key, groupInfFileStoreID)

    sfams_file = get_file(job, "PDB.h5", pdbFileStoreID)
    sfams = pd.read_hdf(unicode(sfams_file), "Superfamilies", where="sfam_id=={}".format(float(mol_sfam_id)))
    sfams = sfams[["sdi", "label"]] #[sfams["sfam_id"]==float(mol_sfam_id)]

    convert_residues = lambda row: decode_residues(row["pdbId"], row["chnLett"], row["resi"], row)

    inferred_interfaces = dd.read_hdf([unicode(tmp_df_file)], "/table")
    inferred_interfaces = inferred_interfaces.repartition(npartitions=NUM_WORKERS)
    resi = inferred_interfaces.apply(convert_residues, axis=1, meta=pd.Series({"resi":[str]}))
    inferred_interfaces = inferred_interfaces.assign(resi=resi)

    obspath = get_file(job, "{}.observed_interactome".format(mol_sfam_id), out_store)
    try:
        observed_interactome = pd.read_hdf(unicode(obspath), "table")
    except TypeError:
        job.log("Failed reading {}".format(obspath))
        raise

    #Rename mol to nbr
    observed_interactome.rename(columns={
        'mol_sdi': 'nbr_sdi',
        'mol_domNo': 'nbr_domNo',
        'mol_gi':'nbr_gi',
        'mol_pdb': 'nbr_pdb',
        'mol_chain': 'nbr_chain'})

    #Add in neghibor information from observed interactome
    inferred_interfaces = inferred_interfaces.merge(observed_interactome,
        how="left", left_on="nbr_obs_int_id", right_on="obs_int_id",
        suffixes=["_inf", "_obs"])
    del observed_interactome

    #Select relevant columns
    if "int_superfam_id_inf" in inferred_interfaces.columns:
        int_superfam_id_col = "int_superfam_id_inf"
    elif "int_superfam_id_x" in inferred_interfaces.columns:
        #Suffix not working?
        int_superfam_id_col = "int_superfam_id_x"
    else:
        raise RuntimeError("Merge faled for obs and inf")
    try:
        inferred_interfaces = inferred_interfaces[["sdi", "nbr_sdi_id", "nbr_score", "nbr_taxid",
            "nbr_superfam_id", int_superfam_id_col, "nbr_obs_int_id", "resn", "resi",
            "domNo", "pdbId", "chnLett", "from", "to", "int_sdi", "int_taxid", "int_res",
            "int_domNo", "int_pdb", "int_chain", "int_sdi_from", "int_sdi_to",
            "nbr_sdi", "nbr_domNo", "nbr_gi", "nbr_pdb", "nbr_chain"]]
    except KeyError as e:
        job.log("Unable to filter df. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
        raise

    #Rename columns
    inferred_interfaces = inferred_interfaces.rename(columns={
        int_superfam_id_col:"int_superfam_id",
        "resn":"mol_resn",
        "sdi":"mol_sdi",
        "resi":"mol_resi",
        "domNo":"mol_domNo",
        "pdbId":"mol_pdb",
        "chnLett":"mol_chain",
        "from":"mol_sdi_from",
        "to":"mol_sdi_to"})

    #Add superfamily name
    #inferred_interfaces = inferred_interfaces.merge(sfams, how="inner", left_on="mol_sdi", right_on="sdi")
    #inferred_interfaces.rename(columns={"label":"mol_superfam_label"})

    #Add interacting superfamily name
    #inferred_interfaces = inferred_interfaces.merge(sfams, how="inner", left_on="int_sdi", right_on="sdi")
    #inferred_interfaces.rename(columns={"label":"int_superfam_label"})

    try:
        inferred_interfaces = inferred_interfaces[ \
            ~inferred_interfaces["mol_sdi"].isnull() & \
            ~inferred_interfaces["int_sdi"].isnull()]
    except KeyError as e:
        job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
        job.log("Sfam cols: {}".format(sfams.columns))
        raise

    try:
        inferred_interfaces = inferred_interfaces.compute(scheduler="multiprocessing", num_workers=cores)
    except Exception as e:
        print e
        job.log(e)
        raise

    df_file = os.path.join(work_dir, df_key)
    inferred_interfaces.to_hdf(unicode(df_file), "table", format="table",
        table=True, complevel=9, complib="bzip2", min_itemsize=1024,
        data_coumns=["nbr_obs_int_id", "nbr_sdi", "mol_sdi", "int_sdi"],
        #scheduler="multiprocessing", dask_kwargs={"num_workers":NUM_WORKERS}
        )
    job.log("Wrote "+df_file)

    #Add ibis info into local job store
    newIbisInfFileStoreID = job.fileStore.writeGlobalFile(df_file)

    save_cols = ["mol_pdb", "mol_chain", "mol_resi"]
    resi = [job.addChildJobFn(convert_residues, row, only_mol=True).rv() for \
        row in inferred_interfaces[save_cols].itertuples()]

    #Merge converted residues
    job.addFollowOnJobFn(merge_interactome_resi, sfam_id, resi, pdbFileStoreID, \
        newIbisInfFileStoreID)

def get_inferred_structural_interactome_by_table(job, table, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    pdb_file = get_file(job, "PDB.h5", pdbFileStoreID)
    struct_domains = pd.read_hdf(unicode(pdb_file), "StructuralDomains")

    #Read in H5 for entire table
    infpath = get_file(job, "IBIS_inferred_{}.h5".format(table), in_store)

    try:
        inferred_interfaces = pd.read_hdf(unicode(infpath), "Intrac{}".format(table))
    except KeyError as e:
        job.log("Unable to load inferred interfaces from {}: {}".format(infpath, e))
        return

    #Add PDB, chain, and sdi information. Inner join to only allow sdi's that are in databses, not obsolete
    inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="inner", left_on="mol_sdi_id", right_on="sdi")

    if inferred_interfaces.shape[0] == 0:
        return

    #Add column to group by
    inferred_interfaces["superfamily"] = inferred_interfaces["nbr_superfam_id"]

    submitted_jobs = 0
    for mol_sfam_id, group in inferred_interfaces.groupby("superfamily"):
        outpath = os.path.join(work_dir, "Intrac{}_{}.inferred_interactome".format(table, mol_sfam_id))
        if not os.path.isfile(outpath):
            group = group.copy()
            group.to_hdf(unicode(outpath), "table", format="table", table=True,
                complevel=9, complib="bzip2", min_itemsize=756)
            #Add group into local job store
            groupInfFileStoreID = job.fileStore.writeGlobalFile(outpath)
            job.addChildJobFn(process_inferred_sfam_table, mol_sfam_id, table, groupInfFileStoreID, pdbFileStoreID)
            submitted_jobs += 1
    job.log("Table {}: Submitted {} Superfamilies".format(table, submitted_jobs))

def merge_inferred_interactome(job, sfam_id):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    iostore = IOStore.get("{}:molmimic-ibis".format(prefix))

    sfam_prefix = "{}/Intrac".format(int(sfam_id))

    sfam_file = "{s}/{s}.inferred_interactome".formart(s=int(sfam_id))

    for table_prefix in iostore.list_input_directory(sfam_prefix):
        table_file = os.path.join(work_dir, table_prefix)
        iostore.read_input_file(table_prefix, table_file)
        try:
            for df in pd.read_hdf(unicode(table_file), "table", chunksize=1000):
                df[cols] = df[cols].astype(float)
                df.to_hdf(unicode(sfam_id), "table", mode="a", format="table",
                    table=True, complevel=9, complib="bzip2", min_itemsize=1024)
        except (IOError, ValueError):
            job.log("Failed to read {}".format(table_file))
        iostore.remove_file(table_prefix)

def start_toil_observed(job, pdbFileStoreID):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ibis".format(prefix))

    ibis_obs_prefix = "IBIS_observed.h5"
    ibis_obs_path = os.path.join(work_dir, ibis_obs_prefix)
    in_store.read_input_file(ibis_obs_prefix, ibis_obs_path)

    #Add ibis info into local job store
    ibisObsFileStoreID = job.fileStore.writeGlobalFile(ibis_obs_path)
    
    pdb_path = job.fileStore.readGlobalFile(pdbFileStoreID)
    sfams= pd.read_hdf(unicode(pdb_path), "Superfamilies")
    sfams = sfams["sfam_id"].drop_duplicates().dropna()
    os.remove(pdb_path)
    for sfam_id in sfams:
        job.log("Running {}".format(sfam_id))
        try:
            job.addChildJobFn(get_observed_structural_interactome, int(sfam_id), pdbFileStoreID, ibisObsFileStoreID)
        except ValueError:
            job.log("Cannot convert {} to string".format(sfam_id))
        break #Only test first sfam

def start_toil_inferred(job, pdbFileStoreID):
    for table in xrange(1, 87):
        job.addChildJobFn(get_inferred_structural_interactome_by_table, table, pdbFileStoreID)
    for sfam_id in iter_unique_superfams():
        try:
            job.addFollowJobFn(merge_inferred_interactome, int(sfam_id), pdbFileStoreID)
        except ValueError:
            job.log("Cannot convert {} to string".format(sfam_id))

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
    dataset_name = options.jobStore.split(":")[-1]

    print "Running"

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
