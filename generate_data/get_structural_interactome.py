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

from util import get_interfaces_path, data_path_prefix, iter_cdd
from map_residues import mmdb_to_pdb_resi

NUM_WORKERS = 20
dask.config.set(scheduler='multiprocessing', num_workers=NUM_WORKERS)
dask.config.set(pool=ThreadPool(20))

def decode_residues(pdb, chain, res, row):
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
        return ",".join(map(str, mmdb_to_pdb_resi(pdb, chain, residues)))
    except (IOError, TypeError) as error:
        print "Error mapping mmdb for", pdb, chain, error, row

        residues.insert(0, "mmdb")
        return ",".join(map(str,residues))

def get_observed_structural_interactome(job, dataset_name, sfam_id, cores=NUM_WORKERS):
    ibis_obs_path = os.path.join(data_path_prefix, "IBIS_observed.h5")
    mmdb_path = os.path.join(data_path_prefix, "PDB.h5")

    obs_ints = pd.read_hdf(unicode(ibis_obs_path), "ObsInt")

    #Save PPIs from the given CDD
    obs_ints = obs_ints[obs_ints["mol_superfam_id"]==sfam_id]

    #Read in face1 residues
    face1 = pd.read_hdf(unicode(ibis_obs_path), "MolResFace")
    face1.columns = ["obs_int_id", "mol_res"]

    #Keep entries from current CDD
    obs_ints = pd.merge(obs_ints, face1, how="left", on="obs_int_id")
    del face1

    #Read in face2 residues and convert gzipped asn1 into res numbers
    face2 = pd.read_hdf(unicode(ibis_obs_path), "IntResFace")
    face2.columns = ["obs_int_id", "int_res"]

    #Keep entries from current CDD
    obs_ints = pd.merge(obs_ints, face2, how="left", on="obs_int_id")
    del face2

    st_domain = pd.read_hdf(unicode(mmdb_path), "StructuralDomains")
    st_domain.columns = ['mol_sdi', 'mol_domNo', 'mol_gi', 'mol_pdb', 'mol_chain', 'mol_sdi_from', 'mol_sdi_to']

    obs_ints = pd.merge(obs_ints, st_domain, how="left", left_on="mol_sdi_id", right_on="mol_sdi")

    st_domain.columns = ['int_sdi', 'int_domNo', 'int_gi', 'int_pdb', 'int_chain', 'int_sdi_from', 'int_sdi_to']
    obs_ints = pd.merge(obs_ints, st_domain, how="left", left_on="int_sdi_id", right_on="int_sdi")
    del st_domain

    convert_residues = lambda row: pd.Series({
        "mol_res": decode_residues(row["mol_pdb"], row["mol_chain"], row["mol_res"], row),
        "int_res": decode_residues(row["int_pdb"], row["int_chain"], row["int_res"], row)})

    sfams = pd.read_hdf(unicode(mmdb_path), "Superfamilies")[["sdi", "label", "sfam_id"]].dropna(axis=0, subset=["sfam_id"])
    sfams = sfams[sfams["sfam_id"]==sfam_id].rename(columns={"sdi":"mol_sdi", "label":"int_superfam_label", "sfam_id":"int_superfam_id"})
    obs_ints = pd.merge(obs_ints, sfams, how="left", on="mol_sdi").dropna(axis=0, subset=["mol_sdi"])
    del sfams

    df = dd.from_pandas(obs_ints, name="obs_ints", npartitions=NUM_WORKERS)
    meta = pd.DataFrame({"mol_res":[str], "int_res":[str]})
    obs_resi = df.map_partitions(lambda _df: _df.apply(convert_residues, axis=1), \
        meta=meta).compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    obs_ints.loc[:, ["mol_res", "int_res"]] = obs_resi

    path = os.path.join(get_interfaces_path(dataset_name), "by_superfamily", str(int(sfam_id)), "{}.observed_interactome".format(int(sfam_id)))
    job.log("Writing {} ({}) interactome: {}".format(cdd, sfam_id, path))
    try:
        obs_ints.to_hdf(unicode(path), "table", complevel=9, complib="bzip2")
    except TypeError:
        job.log("Failed writing {} ({}): {}".format(sfam_id, path))
        raise

def process_inferred_cdd(job, mol_sfam_id, dataset_name, table, cores=NUM_WORKERS):
    table = "Intrac{}".format(table)
    path = get_interfaces_path(dataset_name)
    df_file = os.path.join(path, table, "{}.inferred_interactome".format(mol_sfam_id))

    if os.path.isfile(df_file):
        return

    sfams = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "Superfamilies")
    sfams = sfams[sfams[sfam_id]==float(mol_sfam_id)][["sdi", "label"]]

    convert_residues = lambda row: decode_residues(row["pdbId"], row["chnLett"], row["resi"], row)

    inferred_interfaces = dd.read_hdf([unicode(df_file+".tmp")], "/table")
    inferred_interfaces = inferred_interfaces.repartition(npartitions=NUM_WORKERS)
    resi = inferred_interfaces.apply(convert_residues, axis=1, meta=pd.Series({"resi":[str]}))#, scheduler="multiprocessing", num_workers=NUM_WORKERS)
    inferred_interfaces = inferred_interfaces.assign(resi=resi)


    obspath = os.path.join(path, mol_cdd.replace("/", ""), mol_cdd.replace("/", "")+".observed_interactome")
    try:
        observed_interactome = pd.read_hdf(unicode(obspath), "table")
    except TypeError:
        job.log("Failed reading {}".format(obspath))
        raise

    #Add in neghibor information from observed interactome
    #inferred_interfaces = pd.merge(df, observed_interactome, how="left", left_on="nbr_obs_int_id", right_on="obs_int_id", suffixes=["_inf", "_obs"])
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
            "int_domNo", "int_pdb", "int_chain", "int_sdi_from", "int_sdi_to"]]
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

    sfams = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "Superfamilies")
    sfams = sfams[sfams["sfam_id"]==float(mol_sfam_id)]["sdi", "label"]

    #Add superfamily name
    inferred_interfaces = inferred_interfaces.merge(sfams, how="inner", left_on="mol_sdi", right_on="sdi")
    inferred_interfaces.rename(columns={"label":"mol_superfam_label"})

    #Add interacting superfamily name
    inferred_interfaces = inferred_interfaces.merge(sfams, how="inner", left_on="int_sdi", right_on="sdi")
    inferred_interfaces.rename(columns={"label":"int_superfam_label"})

    try:
        #inferred_interfaces.dropna(axis=0, subset=["mol_sdi", "int_sdi"])
        inferred_interfaces = inferred_interfaces[ \
            ~inferred_interfaces["mol_sdi"].isnull() & \
            ~inferred_interfaces["int_sdi"].isnull()]
    except KeyError as e:
        job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
        job.log("Sfam cols: {}".format(sfams.columns))
        raise

    try:
        inferred_interfaces = inferred_interfaces.compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    except Exception as e:
        print e
        job.log(e)
        raise

    inferred_interfaces.to_hdf(unicode(df_file), "table", format="table",
        table=True, complevel=9, complib="bzip2", min_itemsize=1024,
        #scheduler="multiprocessing", dask_kwargs={"num_workers":NUM_WORKERS}
        )
    job.log("Wrote "+df_file)

    return

def get_inferred_structural_interactome_by_table(job, dataset_name, table):
    struct_domains = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "StructuralDomains")

    #Read in H5 for entire table
    infpath = os.path.join(data_path_prefix, "IBIS_inferred.h5")
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
        outpath = os.path.join(get_interfaces_path(dataset_name), "Intrac{}".format(table),
            "{}.inferred_interactome".format(mol_sfam_id))
        if not os.path.isfile(outpath):
            group = group.copy()
            group.to_hdf(unicode(outpath+".tmp"), "table", format="table", table=True, complevel=9, complib="bzip2", min_itemsize=756)
            job.addChildJobFn(process_inferred_cdd, mol_sfam_id, dataset_name, table) #, memory="240000M")
            submitted_jobs += 1
            job.log("FAILED SFAM:{}\t{}".format(mol_sfam_id, table))
    job.log("Table {}: Submitted {} Superfamilies".format(table, submitted_jobs))

def read_inferred_interactome(dataset_name, table, cdd_id):
    print dataset_name, table, cdd_id
    path = os.path.join(get_interfaces_path(dataset_name))
    try:
        df = pd.read_hdf(unicode(os.path.join(path, "Intrac{}".format(table), "{}.inferred_interactome".format(cdd_id))), "table")
    except IOError:
        return None
    cols = ['mol_cdd', 'mol_resn', 'mol_resi', 'mol_pdb', 'mol_chain', 'int_res', 'int_pdb', 'int_chain', 'int_cdd']
    df[cols] = df[cols].apply(lambda s: s.str.encode('utf-8'))
    df[["int_superfam_id", "int_taxid"]] = df[["int_superfam_id", "int_taxid"]].astype(float)
    return df

def merge_inferred_interactome(job, dataset_name, cdd_id, cores=NUM_WORKERS):
    path = get_interfaces_path(dataset_name)
    inferred_file = os.path.join(path, "by_superfamily", str(int(cdd_id)), "{}.inferred_interactome".format(int(cdd_id)))

    if os.path.isfile(inferred_file):
        job.log("Skipping {} bc it exists".format(cdd_id))
        return

    pattern = os.path.join(path, "Intrac*", "{}.inferred_interactome".format(cdd_id))
    files = [unicode(f) for f in glob.glob(pattern) if "IntracGrouped" not in f] # and \
        #os.stat(f).st_size > 5000]
    files = sorted(files, key=lambda f: os.stat(f).st_size)

    while True:
        if len(files)==0:
            job.log("Failed {} bc there are no subfiles".format(cdd_id))
            return
        try:
            df = dd.read_hdf(files, key='/table', mode="r")
            break
        except (IOError, ValueError):
            job.log("Failed to read at least 1 file, removing {} ()".format(files[0], os.stat(files[0]).st_size))
            del files[0]

    df = df.repartition(npartitions=NUM_WORKERS)
    cols = ["mol_sdi","int_superfam_id","mol_domNo","mol_sdi_from","mol_sdi_to","int_taxid"]
    df[cols] = df[cols].astype(float)
    try:
        df.to_hdf(unicode(inferred_file), "/table", complevel=9, complib="bzip2",
            min_itemsize=1024, scheduler="multiprocessing", dask_kwargs={"num_workers":NUM_WORKERS})
        job.log("Successfuly saved {}".format(cdd_id))
    except Exception as e:
        job.log("Failed {} saving due to {}".format(cdd_id, e))

def split_familes(job, dataset_name, cdd_id):
    path = get_interfaces_path(dataset_name)
    inferred_file = os.path.join(path, "IntracGrouped", "{}.inferred_interactome".format(int(cdd_id)))
    interfaces = pd.read_hdf(inferred_file, "table")
    superfamily_name = None

    for families, sfam_id in iter_cdd(use_id=True, group_superfam=True, all_superfam=True):
        for family in families.itertuples():
            cdd_label = family["label"]
            family_interfaces = interfaces[interfaces["mol_cdd"]==cdd_label]
            family_path = os.path.join(path, "by_family", cdd_label, "{}.inferred_interactome".format(cdd_label))
            if not os.path.isdir(os.path.dirname(family_path)):
                os.makedirs(os.path.dirname(family_path))
            family_interfaces.to_hdf(unicode(family_path), "table")

            if superfamily_name:
                #Save first family as superfamily
                superfamily_name = cdd_label

    interfaces["mol_cdd"] = superfamily_name
    interfaces = interfaces.drop_duplicates()
    family_path = os.path.join(path, "by_superfamily", superfamily_name, "{}.inferred_interactome".format(superfamily_name))
    if not os.path.isdir(os.path.dirname(family_path)):
        os.makedirs(os.path.dirname(family_path))
    interfaces.to_hdf(family_path, "table", complevel=9, complib="bzip2", min_itemsize=768,)

def start_toil_observed(job, dataset_name):
    path = os.path.join(get_interfaces_path(dataset_name), "by_superfamily")
    if not os.path.isdir(path):
        os.makedirs(path)
    for cdd, sfam_id in iter_cdd(use_id=True, group_superfam=True):
        sfam_path = os.path.join(path, str(int(sfam_id)), "{}.observed_interactome".format(int(sfam_id)))
        if os.path.isfile(sfam_path):
            continue
        if not os.path.isdir(os.path.dirname(sfam_path)):
            os.makedirs(os.path.dirname(sfam_path))
        job.log("Running {} ({})".format(cdd, sfam_id))
        try:
            job.addChildJobFn(get_observed_structural_interactome, dataset_name, int(sfam_id), memory="120000M")
        except ValueError:
            job.log("Error adding {} (Cannot convert {} to string)".format(cdd, sfam_id))

def start_toil_inferred_merge(job, dataset_name):
    path = os.path.join(get_interfaces_path(dataset_name), "IntracGrouped")
    if not os.path.isdir(path):
        os.makedirs(path)
    for _, sfam_id in iter_cdd(use_id=True, group_superfam=True):
        if os.path.isfile(os.path.join(path, "{}.inferred_interactome".format(int(sfam_id)))):
            continue
        try:
            job.addChildJobFn(merge_inferred_interactome, dataset_name, int(sfam_id), memory="240000M")
        except ValueError:
            job.log("Error adding {} (Cannot convert {} to string)".format(cdd, sfam_id))

def start_toil_inferred(job, dataset_name):
    # stopped_jobs = glob.glob(os.path.join(get_interfaces_path(dataset_name), "Intrac*",
    #     "*.inferred_interactome.tmp"))
    # job.log("Will Issue {} jobs".format(len(stopped_jobs)))
    # if len(stopped_jobs) > 0:
    #     for f in stopped_jobs:
    #         mol_sfam_id = os.path.basename(f).split(".", 1)[0]
    #         table = os.path.basename(os.path.dirname(f))[6:]
    #         if not os.path.isfile(f[:-4]):
    #             job.addChildJobFn(process_inferred_cdd, mol_sfam_id, dataset_name, table)
    #             job.log("Issuing {}".format(f))
    #     job.log("Issued {} jobs".format(len(stopped_jobs)))
    #     return

    path = get_interfaces_path(dataset_name)
    for table in xrange(1, 87):
        table_dir = os.path.join(path, "Intrac{}".format(table))
        if not os.path.isdir(table_dir):
            os.makedirs(table_dir)
        job.addChildJobFn(get_inferred_structural_interactome_by_table, dataset_name, table, memory="120000M")

def start_toil(job, dataset_name):
    if not os.path.isdir(get_interfaces_path(dataset_name)):
        os.makedirs(get_interfaces_path(dataset_name))

    obsjob = job.addChildJobFn(start_toil_observed, dataset_name)
    infjob = obsjob.addFollowOnJobFn(start_toil_inferred, dataset_name)
    infmergejob = infjob.addFollowOnJobFn(start_toil_inferred_merge, dataset_name, cores=20, memory="240000M")
    return

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    dataset_name = options.jobStore.split(":")[-1]

    print "Running"

    job = Job.wrapJobFn(start_toil, dataset_name)
    with Toil(options) as toil:
        toil.start(job)
