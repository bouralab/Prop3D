import os, sys
import glob

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

NUM_WORKERS = 16
dask.config.set(scheduler='multiprocessing', num_worker=NUM_WORKERS)

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
            raise

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

tbl_map = None
def get_tbl(sdi):
    if tbl_map is None:
        tbl_map = pd.read_sql("SELECT tbl_id, min_sdi, max_sdi FROM Intrac.dbo.TblMap", conn)
    return tbl_map[(tbl_map["min_sdi"]<=sdi)|(tbl_map["max_sdi"]>=sdi)].loc[0, "tbl_id"]

def get_observed_structural_interactome(job, dataset_name, cdd, sfam_id, cores=NUM_WORKERS):
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
    sfams = sfams[sfams["sfam_id"]==sfam_id].rename(columns={"sdi":"mol_sdi", "label":"int_cdd", "sfam_id":"int_superfam_id"})
    obs_ints = pd.merge(obs_ints, sfams, how="left", on="mol_sdi").dropna(axis=0, subset=["mol_sdi"])
    del sfams

    df = dd.from_pandas(obs_ints, name="obs_ints", npartitions=NUM_WORKERS)
    meta = pd.DataFrame({"mol_res":[str], "int_res":[str]})
    obs_resi = df.map_partitions(lambda _df: _df.apply(convert_residues, axis=1), \
        meta=meta).compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    obs_ints.loc[:, ["mol_res", "int_res"]] = obs_resi

    cdd = cdd.replace("/", "")
    path = os.path.join(get_interfaces_path(dataset_name), cdd, "{}.observed_interactome".format(cdd))
    job.log("Writing {} ({}) interactome: {}".format(cdd, sfam_id, path))
    try:
        obs_ints.to_hdf(unicode(path), "table", complevel=9, complib="bzip2")
    except TypeError:
        job.log("Failed writing {} ({}): {}".format(cdd, sfam_id, path))
        raise

def process_inferred_cdd(job, mol_sfam_id, df_file, dataset_name, table, cores=NUM_WORKERS):
    table = "Intrac{}".format(table)
    path = get_interfaces_path(dataset_name)

    outpath = os.path.join(path, table, "{}.inferred_interactome".format(mol_sfam_id))
    if os.path.isfile(outpath):
        return pd.Series({"Done":mol_sfam_id})

    mol_cdd = next(iter_cdd(id=float(mol_sfam_id)))

    convert_residues = lambda row: decode_residues(row["pdbId"], row["chnLett"], row["resi"], row)
    #convert_row = lambda df: df.apply(convert_residues, axis=1)

    inferred_interfaces = dd.read_hdf(df_file, "/table")
    inferred_interfaces = inferred_interfaces.repartition(npartitions=NUM_WORKERS)
    inferred_interfaces = inferred_interfaces.assign(mol_cdd=mol_cdd)
    print "Start"
    resi = inferred_interfaces.apply(convert_residues, axis=1, meta=pd.Series({"resi":[str]}))#, scheduler="multiprocessing", num_workers=NUM_WORKERS)
    print "End"
    inferred_interfaces = inferred_interfaces.assign(resi=resi)

    # resi = ddf.map_partitions(lambda _df: _df.apply(convert_residues, axis=1), \
    #     meta=meta).compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    # df.loc[:, ["resi"]] = resi

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
            "nbr_superfam_id", int_superfam_id_col, "mol_cdd", "nbr_obs_int_id", "resn", "resi",
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

    #Add interacting superfamily name
    sfams = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "Superfamilies")[["label", "sfam_id"]].drop_duplicates()
    sfams = sfams.rename(columns={"sfam_id":"int_superfam_id", "label":"int_cdd"})
    inferred_interfaces = inferred_interfaces.merge(sfams, how="left", on="int_superfam_id")

    try:
        #inferred_interfaces.dropna(axis=0, subset=["mol_sdi", "int_sdi"])
        inferred_interfaces = inferred_interfaces[~inferred_interfaces["mol_sdi"].isnull()&\
            ~inferred_interfaces["int_sdi"].isnull()]
    except KeyError as e:
        job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
        job.log("Sfam cols: {}".format(sfams.columns))
        raise

    print "Save File", inferred_interfaces.npartitions
    inferred_interfaces.to_hdf(unicode(outpath), "table", format="table",
        table=True, complevel=9, complib="bzip2", min_itemsize=756,
        )
    job.log("Wrote "+outpath)

    return pd.Series({"Done":mol_sfam_id})

def get_inferred_structural_interactome_by_table(job, dataset_name, table):
    struct_domains = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "StructuralDomains")

    #Read in H5 for entire table
    infpath = os.path.join(data_path_prefix, "IBIS_inferred.h5")
    try:
        inferred_interfaces = pd.read_hdf(unicode(infpath), "Intrac{}".format(table))
    except KeyError as e:
        job.log("Unable to load inferred interfaces from {}: {}".format(infpath, e))
        return

    #Add PDB, chain, and sdi information
    inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="left", left_on="mol_sdi_id", right_on="sdi")

    #Add column to group by
    inferred_interfaces["superfamily"] = inferred_interfaces["nbr_superfam_id"]

    #Change job queue
    #old_slurm_args = os.environ["TOIL_SLURM_ARGS"]
    #os.environ["TOIL_SLURM_ARGS"] = "-t 10:00:00 --partition=largemem -A muragroup"

    submitted_jobs = 0
    for mol_sfam_id, group in inferred_interfaces.groupby("superfamily"):
        outpath = os.path.join(get_interfaces_path(dataset_name), "Intrac{}".format(table),
            "{}.inferred_interactome".format(mol_sfam_id))
        if not os.path.isfile(outpath):
            group = group.copy()
            group.to_hdf(unicode(outpath+".tmp"), "table", format="table", table=True, complevel=9, complib="bzip2", min_itemsize=756)
            job.addChildJobFn(process_inferred_cdd, mol_sfam_id, outpath+".tmp", dataset_name, table) #, memory="240000M")
            submitted_jobs += 1
    job.log("Table {}: Submitted {} Superfamilies".format(table, submitted_jobs))

    #os.environ["TOIL_SLURM_ARGS"] = old_slurm_args

    # sfams = pd.read_hdf(unicode(os.path.join(data_path_prefix, "PDB.h5")), "Superfamilies")[["label", "sfam_id"]].drop_duplicates()
    # sfams = sfams.rename(columns={"sfam_id":"int_superfam_id", "label":"int_cdd"})
    #
    # convert_residues = lambda row: pd.Series({
    #     "resi": decode_residues(row["pdbId"], row["chnLett"], row["resi"], row)})
    # convert_row = lambda df: df.apply(convert_residues, axis=1)
    #
    # def apply_df(df):
    #     df = df.copy()
    #     try:
    #         mol_sfam_id = df["nbr_superfam_id"].iloc[0]
    #     except KeyError:
    #         print "FAILED", df.index
    #         return pd.Series({"Done":None})
    #
    #     outpath = os.path.join(path, table, "{}.inferred_interactome".format(mol_sfam_id))
    #     if os.path.isfile(outpath):
    #         return pd.Series({"Done":mol_sfam_id})
    #
    #     mol_cdd = next(iter_cdd(id=float(mol_sfam_id)))
    #     df.loc[:, "mol_cdd"] = mol_cdd
    #
    #     inf_resi = df.apply(convert_residues, axis=1)
    #     df.loc[:, "resi"] = inf_resi
    #
    #     obspath = os.path.join(path, mol_cdd.replace("/", ""), mol_cdd.replace("/", "")+".observed_interactome")
    #     try:
    #         observed_interactome = pd.read_hdf(unicode(obspath), "table")
    #     except TypeError:
    #         job.log("Failed reading {}".format(obspath))
    #         raise
    #
    #     #Add in neghibor information from observed interactome
    #     inferred_interfaces = pd.merge(df, observed_interactome, how="left", left_on="nbr_obs_int_id", right_on="obs_int_id", suffixes=["_inf", "_obs"])
    #
    #     #Select relevant columns
    #     if "int_superfam_id_inf" in inferred_interfaces.columns:
    #         int_superfam_id_col = "int_superfam_id_inf"
    #     elif "int_superfam_id_x" in inferred_interfaces.columns:
    #         #Suffix not working?
    #         int_superfam_id_col = "int_superfam_id_x"
    #     else:
    #         raise RuntimeError("Merge faled for obs and inf")
    #     try:
    #         inferred_interfaces = inferred_interfaces[["sdi", "nbr_sdi_id", "nbr_score", "nbr_taxid",
    #             "nbr_superfam_id", int_superfam_id_col, "mol_cdd", "nbr_obs_int_id", "resn", "resi",
    #             "domNo", "pdbId", "chnLett", "from", "to", "int_sdi", "int_taxid", "int_res",
    #             "int_domNo", "int_pdb", "int_chain", "int_sdi_from", "int_sdi_to"]]
    #     except KeyError as e:
    #         job.log("Unable to filter df. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
    #         raise
    #
    #     #Rename columns
    #     inferred_interfaces = inferred_interfaces.rename(columns={
    #         int_superfam_id_col:"int_superfam_id",
    #         "resn":"mol_resn",
    #         "sdi":"mol_sdi",
    #         "resi":"mol_resi",
    #         "domNo":"mol_domNo",
    #         "pdbId":"mol_pdb",
    #         "chnLett":"mol_chain",
    #         "from":"mol_sdi_from",
    #         "to":"mol_sdi_to"})
    #
    #     #Add interacting superfamily name
    #     inferred_interfaces = pd.merge(inferred_interfaces, sfams, how="left", on="int_superfam_id")
    #     try:
    #         inferred_interfaces.dropna(axis=0, subset=["mol_sdi", "int_sdi"])
    #     except KeyError as e:
    #         job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
    #         job.log("Sfam cols: {}".format(sfams.columns))
    #         raise
    #
    #     inferred_interfaces.to_hdf(unicode(outpath), "table", complevel=9, complib="bzip2", min_itemsize=756)
    #     job.log("Wrote "+outpath)
    #
    #     return pd.Series({"Done":mol_sfam_id})

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

def merge_inferred_interactome(job, dataset_name, cdd_name, cdd_id, cores=NUM_WORKERS):
    cdd_name = cdd_name.replace("/", "")
    path = get_interfaces_path(dataset_name)
    inferred_file = os.path.join(path, cdd_name, "{}.inferred_interactome".format(cdd_name))

    if os.path.isfile(inferred_file):
        return

    df = dd.read_hdf(os.path.join(path, "Intrac*", "{}.inferred_interactome".format(cdd_id)), key='/table')
    df.to_hdf(inferred_file, "/table", complevel=9, complib="bzip2", min_itemsize=768)

    # inferred_interactome = Parallel(n_jobs=NUM_WORKERS)(delayed(read_inferred_interactome)(dataset_name, table, cdd_id) \
    #     for table in xrange(1, 87))
    # inferred_interactome = pd.concat(filter(lambda df: df is not None, inferred_interactome))
    # inferred_interactome.to_hdf(unicode(tmp_file), "table", complevel=9, complib="bzip2", min_itemsize=768)

    #
    # for table in xrange(1, 87):
    #     print table
    #     input = os.path.join(path, "Intrac{}".format(table), "{}.inferred_interactome".format(cdd_id))
    #     try:
    #         df = pd.read_hdf(unicode(input), "table")
    #     except IOError:
    #         job.log("Can't open {}".format(os.path.join(path, "Intrac{}".format(table), "{}.inferred_interactome".format(cdd_id))))
    #         continue
    #     print "Done reading"
    #     #print df.dtypes
    #
    #     # cols = ['mol_cdd', 'mol_resn', 'mol_resi', 'mol_pdb', 'mol_chain', 'int_res', 'int_pdb', 'int_chain', 'int_cdd']
    #     # df[cols] = df[cols].apply(lambda s: s.str.encode('utf-8'))
    #     #
    #     # #for col in ['mol_cdd', 'mol_resn', 'mol_resi', 'mol_pdb', 'mol_chain', 'int_res', 'int_pdb', 'int_chain', 'int_cdd']:
    #     # #    df[col] = df[col].str.encode('utf-8')
    #     # df[["int_superfam_id", "int_taxid"]] = df[["int_superfam_id", "int_taxid"]].astype(float)
    #     #
    #     # print df.dtypes
    #
    #     df.to_hdf(unicode(input+".table"), "table", table=True, complevel=9, complib="bzip2", min_itemsize=768)
    #     del df

    # return
    #
    # df = pd.read_hdf(unicode(tmp_file), "table").reset_index()
    # df.to_hdf(unicode(tmp_file), "table", table=True, mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)

    # shutil.move(
    #     os.path.join(path, cdd_name, "{}.inferred_interactome.tmp".format(cdd_name)),
    #     os.path.join(path, cdd_name, "{}.inferred_interactome".format(cdd_name))
    # )

def start_toil_observed(job, dataset_name):
    for cdd, sfam_id in iter_cdd(use_id=True):
        path = os.path.join(get_interfaces_path(dataset_name), "{}.observed_interactome".format(cdd.replace("/", "")))
        if os.path.isfile(path):
            continue
        job.log("Running {} ({})".format(cdd, sfam_id))
        try:
            job.addChildJobFn(get_observed_structural_interactome, dataset_name, cdd, int(sfam_id), memory="120000M")
        except ValueError:
            job.log("Error adding {} (Cannot convert {} to string)".format(cdd, sfam_id))

def start_toil_inferred_merge(job, dataset_name):
    path = os.path.join(get_interfaces_path(dataset_name))
    for cdd, sfam_id in iter_cdd(use_id=True):
        cdd_name = cdd.replace("/", "")
        if os.path.isfile(os.path.join(path, cdd_name, "{}.inferred_interactome".format(cdd_name))):
            continue
        try:
            job.addChildJobFn(merge_inferred_interactome, dataset_name, cdd, int(sfam_id))
        except ValueError:
            job.log("Error adding {} (Cannot convert {} to string)".format(cdd, sfam_id))

def start_toil_inferred(job, dataset_name):
    stopped_jobs = glob.glob(os.path.join(get_interfaces_path(dataset_name), "Intrac*",
        "*.inferred_interactome.tmp"))
    if len(stopped_jobs) > 0:
        for f in stopped_jobs:
            mol_sfam_id = os.path.basename(f).split(".", 1)[0]
            table = os.path.basename(os.path.dirname(f))[6:]
            job.addChildJobFn(process_inferred_cdd, mol_sfam_id, outpath+".tmp", dataset_name, table)
        return
        
    path = get_interfaces_path(dataset_name)
    for table in xrange(1, 87):
        table_dir = os.path.join(path, "Intrac{}".format(table))
        if not os.path.isdir(table_dir):
            os.makedirs(table_dir)
        job.addChildJobFn(get_inferred_structural_interactome_by_table, dataset_name, table, memory="240000M")
    #job.addFollowOnJobFn(start_toil_inferred_merge, dataset_name, memory="120000M")

def start_toil(job, dataset_name):
    if not os.path.isdir(get_interfaces_path(dataset_name)):
        os.makedirs(get_interfaces_path(dataset_name))

    #obsjob = job.addChildJobFn(start_toil_observed, dataset_name)
    infjob = job.addChildJobFn(start_toil_inferred, dataset_name)
    infmergejob = infjob.addFollowOnJobFn(start_toil_inferred_merge, dataset_name, cores=20, memory="240000M")

def fix2(job, f):
    try:
        df = pd.read_hdf(unicode(f), "table")
    except IOError:
        job.log("Can't open {}".format(os.path.join(path, "Intrac{}".format(table), "{}.inferred_interactome".format(cdd_id))))
        return
    print "Done reading"
    df.to_hdf(unicode(f), "table", format="table", table=True, complevel=9, complib="bzip2", min_itemsize=768)

def fix(job, dataset_name):
    import glob
    fixing = 0
    for f in glob.glob(os.path.join(get_interfaces_path(dataset_name), "Intrac*", "*.inferred_interactome")):
        with pd.HDFStore(f, mode="a") as hdf:
            stops = []
            divisions = []
            storer = hdf.get_storer("table")
            if storer.format_type != 'table':
                job.addChildJobFn(fix2, f)
                fixing += 1
    job.log("WILL FIX {} HDFs".format(fixing))


if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    dataset_name = options.jobStore.split(":")[-1]

    job = Job.wrapJobFn(start_toil, dataset_name)
    with Toil(options) as toil:
        toil.start(job)

# def submit_jobs(dataset_name, job_type, job_name="interactome", dependency=None):
#     if not os.path.isdir(get_interfaces_path(dataset_name)):
#         os.makedirs(get_interfaces_path(dataset_name))
#
#     job = SlurmJob(job_name+"_"+job_type, cpus=56, walltime="8:00:00")
#     for cdd, sfam_id in iter_cdd(use_id=True):
#         try:
#             job += "{} {} {} {} {} {}\n".format(sys.executable, __file__, job_type, dataset_name, cdd, int(sfam_id))
#         except ValueError:
#             print "Error for", cdd
#     jid = job.run(dependency=dependency)
#     print jid
#     return jid
#
# def submit_jobs_observed(dataset_name, job_name="interactome", dependency=None):
#     return submit_jobs(dataset_name, "obs", job_name=job_name, dependency=dependency)
#
# def submit_jobs_inferred(dataset_name, job_name="interactome", dependency=None):
#     job = SlurmJob(job_name+"_inf", cpus=56, walltime="8:00:00")
#     path = get_interfaces_path(dataset_name)
#     for table in xrange(1, 87):
#         table_dir = os.path.join(path, "Intrac{}".format(table))
#         if not os.path.isdir(table_dir):
#             os.makedirs(table_dir)
#         job += "{} {} inf {} {}\n".format(sys.executable, __file__, dataset_name, table)
#     jid = job.run(dependency=dependency)
#     return submit_jobs(dataset_name, "merge_inf", dependency="afterany:"+str(jid))

# def old_interactome():
#     interfaces_file = os.path.join(get_interfaces_path(dataset_name), "{}.extended_sdi".format(cdd.replace("/", "")))
#     interfaces = pd.read_table(interfaces_file, dtype="str")
#     observed = interfaces[interfaces["observed"]=="1"]
#     pd.merge(observed, obs_ints, left_on=("sdi", "resi"), right_on=("mol_sdi_id", "mol_res"), how="outer", validate="one_to_one")
#     inferred = interfaces[interfaces["observed"]=="0"]
#     inferred["tbl"] = inferred["sdi"].apply(get_tbl)
#     inferred["db"] = inferred["tbl"].apply(lambda t: int(t/5))
#     groups = inferred.group_by(["db", "tbl"])
#     for group_name, group in groups:
#         intrac = pd.read_sql("SELECT inf_id, mol_sdi_id, nbr_sdi_id, nbr_obs_int_id, face_seq FROM Inrac{db}.dbo.Inferred{tbl} WHERE nbr_superfam_id={sfam}".format(
#             db=group_name[0], tbl=group_name[1], sfam=sfam_id), conn)
#         group = pd.merge(group, intrac, how="left", left_on=("sdi", "resn"), left_on=("mol_sdi_id", "face_seq"))
#         for _, row in group.iterrows():
#             inf_obs = observed[observed["obs_int_int"]==row["nbr_obs_int_id"]]
#     with open(interfaces_file) as f:
#         f.next()
#         for line in f:
#             try:
#                 pdb, chain, sdi1, domNo1, resi1, resn1, is_multimer, cdd, partner_cdd, pdb_evidence, observed = line.rstrip().split("\t")
#             except ValueError:
#                 continue
#             tbl = get_tbl(sdi1)

# def get_inferred_structural_interactome_by_cdd(dataset_name, cdd, sfam_id):
#     print cdd, sfam_id
#     path = os.path.join(get_interfaces_path(dataset_name), cdd.replace("/", ""))
#     struct_domains = pd.read_hdf(os.path.join(data_path_prefix, "PDB.h5"), "StructuralDomains")
#
#     observed_interactome = pd.read_hdf(path+".observed_interactome", "table")
#
#     convert_residues = lambda row: pd.Series({
#         "resi": decode_residues(row["pdbId"], row["chnLett"], row["resi"])})
#
#     for table in xrange(1, 87):
#         print table
#         #Read in H5 for entire table
#         inferred_interfaces = pd.read_hdf(os.path.join(data_path_prefix, "IBIS_inferred.h5"), "Intrac{}".format(table))
#
#         #Filter table by CDD
#         inferred_interfaces = inferred_interfaces[inferred_interfaces["nbr_superfam_id"]==sfam_id]
#
#         #Add PDB, chain, and sdi information
#         inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="left", left_on="mol_sdi_id", right_on="sdi")
#
#         #Convert resn to PDB format
#         df = dd.from_pandas(inferred_interfaces, name="inferred_interfaces", npartitions=56)
#         meta = pd.DataFrame({"resi":[str]})
#         inf_resi = df.map_partitions(lambda _df: _df.apply(convert_residues, axis=1), meta=meta).compute(get=get)
#         inferred_interfaces["resi"] = inf_resi
#
#         #Add in neghibor information from observed interactome
#         inferred_interfaces = pd.merge(inferred_interfaces, observed_interactome, how="left", left_on="nbr_obs_int_id", right_on="obs_int_id", suffixes=["_inf", "_obs"])
#
#         #Select relevant columns
#         inferred_interfaces = inferred_interfaces[["sdi", "nbr_sdi_id", "nbr_score", "nbr_taxid",
#             "nbr_superfam_id", "int_superfam_id_inf", "nbr_obs_int_id", "resn", "resi",
#             "domNo", "pdbId", "chnLett", "from", "to", "int_sdi", "int_taxid", "int_res",
#             "int_domNo", "int_pdb", "int_chain", "int_sdi_from", "int_sdi_to"]]
#
#         #Rename columns
#         inferred_interfaces = inferred_interfaces.rename(columns={
#             "int_superfam_id_inf":"int_superfam_id",
#             "resn":"mol_resn", "sdi":"mol_sdi",
#             "resi":"mol_resi",
#             "domNo":"mol_domNo",
#             "pdbId":"mol_pdb",
#             "chnLett":"mol_chain",
#             "from":"mol_sdi_from",
#             "to":"mol_sdi_to"})
#
#         inferred_interfaces.to_hdf(path+".inferred_interactome.tmp", "table", table=True, mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=756)
#
#
#     try:
#         shutil.move(path+".inferred_interactome.tmp", path+".inferred_interactome")
#     except IOError:
#         print "Maybe no entries?"
#
# if __name__ == "__main__":
#     if len(sys.argv) == 3 and sys.argv[1] == "obs":
#         submit_jobs_observed(sys.argv[2])
#     elif len(sys.argv) == 3 and sys.argv[1] == "inf":
#         submit_jobs_inferred(sys.argv[2])
#     elif len(sys.argv) in [4, 5] and sys.argv[1] in ["obs", "merge_inf"]:
#         dataset_name, cdd = sys.argv[2:4]
#         if len(sys.argv) == 4:
#             CDD = iter_cdd(use_id=True, label=cdd)
#             _, sfam_id = next(CDD)
#         else:
#             sfam_id = sys.argv[4]
#
#         if sys.argv[1] == "obs":
#             get_observed_structural_interactome(dataset_name, cdd, int(sfam_id))
#         else:
#             merge_inferred_interactome(dataset_name, cdd, int(sfam_id))
#     elif len(sys.argv) == 4 and sys.argv[1] == "inf":
#         get_inferred_structural_interactome_by_table(sys.argv[2], sys.argv[3])
