import os, sys
sys.path.append("/data/draizene/molmimic")

import binascii
import zlib
from pyasn1.codec.ber import decoder
import shutil

import pandas as pd
import dask.dataframe as dd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(*args, **kwds):
        return args


from util import get_interfaces_path, data_path_prefix, iter_cdd
from map_residues import mmdb_to_pdb_resi

NUM_WORKERS = 4

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

    obs_ints = pd.read_hdf(ibis_obs_path, "ObsInt")

    #Save PPIs from the given CDD
    obs_ints = obs_ints[obs_ints["mol_superfam_id"]==sfam_id]

    #Read in face1 residues
    face1 = pd.read_hdf(ibis_obs_path, "MolResFace")
    face1.columns = ["obs_int_id", "mol_res"]

    #Keep entries from current CDD
    obs_ints = pd.merge(obs_ints, face1, how="left", on="obs_int_id")
    del face1

    #Read in face2 residues and convert gzipped asn1 into res numbers
    face2 = pd.read_hdf(ibis_obs_path, "IntResFace")
    face2.columns = ["obs_int_id", "int_res"]

    #Keep entries from current CDD
    obs_ints = pd.merge(obs_ints, face2, how="left", on="obs_int_id")
    del face2

    st_domain = pd.read_hdf(mmdb_path, "StructuralDomains")
    st_domain.columns = ['mol_sdi', 'mol_domNo', 'mol_gi', 'mol_pdb', 'mol_chain', 'mol_sdi_from', 'mol_sdi_to']

    obs_ints = pd.merge(obs_ints, st_domain, how="left", left_on="mol_sdi_id", right_on="mol_sdi")

    st_domain.columns = ['int_sdi', 'int_domNo', 'int_gi', 'int_pdb', 'int_chain', 'int_sdi_from', 'int_sdi_to']
    obs_ints = pd.merge(obs_ints, st_domain, how="left", left_on="int_sdi_id", right_on="int_sdi")
    del st_domain

    convert_residues = lambda row: pd.Series({
        "mol_res": decode_residues(row["mol_pdb"], row["mol_chain"], row["mol_res"], row),
        "int_res": decode_residues(row["int_pdb"], row["int_chain"], row["int_res"], row)})

    sfams = pd.read_hdf(mmdb_path, "Superfamilies")[["sdi", "label", "sfam_id"]].dropna(axis=0, subset=["sfam_id"])
    sfams = sfams[sfams["sfam_id"]==sfam_id].rename(columns={"sdi":"mol_sdi", "label":"int_cdd", "sfam_id":"int_superfam_id"})
    obs_ints = pd.merge(obs_ints, sfams, how="left", on="mol_sdi").dropna(axis=0, subset=["mol_sdi"])
    del sfams

    df = dd.from_pandas(obs_ints, name="obs_ints", npartitions=NUM_WORKERS)
    meta = pd.DataFrame({"mol_res":[str], "int_res":[str]})
    obs_resi = df.map_partitions(lambda _df: _df.apply(convert_residues, axis=1), \
        meta=meta).compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    obs_ints.loc[:, ["mol_res", "int_res"]] = obs_resi

    path = os.path.join(get_interfaces_path(dataset_name), "{}.observed_interactome".format(cdd.replace("/", "")))
    job.log("Writing {} ({}) interactome: {}".format(cdd, sfam_id, path))
    try:
        obs_ints.to_hdf(unicode(path), "table", complevel=9, complib="bzip2")
    except TypeError:
        job.log("Failed writing {} ({}): {}".format(cdd, sfam_id, path))
        raise

def get_inferred_structural_interactome_by_table(job, dataset_name, table, cores=NUM_WORKERS):
    table = "Intrac{}".format(table)
    path = get_interfaces_path(dataset_name)

    struct_domains = pd.read_hdf(os.path.join(data_path_prefix, "PDB.h5"), "StructuralDomains")

    sfams = pd.read_hdf(os.path.join(data_path_prefix, "PDB.h5"), "Superfamilies")[["label", "sfam_id"]].drop_duplicates()
    sfams = sfams.rename(columns={"sfam_id":"int_superfam_id", "label":"int_cdd"})

    convert_residues = lambda row: pd.Series({
        "resi": decode_residues(row["pdbId"], row["chnLett"], row["resi"], row)})
    convert_row = lambda df: df.apply(convert_residues, axis=1)

    def apply_df(df):
        df = df.copy()
        try:
            mol_sfam_id = df["nbr_superfam_id"].iloc[0]
        except KeyError:
            print "FAILED", df.index
            return pd.Series({"Done":None})

        mol_cdd = next(iter_cdd(id=float(mol_sfam_id)))
        df.loc[:, "mol_cdd"] = mol_cdd

        inf_resi = df.apply(convert_residues, axis=1)
        df.loc[:, "resi"] = inf_resi

        obspath = os.path.join(path, mol_cdd.replace("/", "")+".observed_interactome")
        try:
            observed_interactome = pd.read_hdf(unicode(obspath), "table")
        except TypeError:
            job.log("Failed reading {}".format(obspath))
            raise

        #Add in neghibor information from observed interactome
        inferred_interfaces = pd.merge(df, observed_interactome, how="left", left_on="nbr_obs_int_id", right_on="obs_int_id", suffixes=["_inf", "_obs"])

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
        inferred_interfaces = pd.merge(inferred_interfaces, sfams, how="left", on="int_superfam_id")
        try:
            inferred_interfaces.dropna(axis=0, subset=["mol_sdi", "int_sdi"])
        except KeyError as e:
            job.log("Unable to drop na. Columns are: {}. Error is: {}".format(inferred_interfaces.columns, e))
            job.log("Sfam cols: {}".format(sfams.columns))
            raise

        outpath = os.path.join(path, table, "{}.inferred_interactome".format(mol_sfam_id))
        inferred_interfaces.to_hdf(unicode(outpath), "table", complevel=9, complib="bzip2", min_itemsize=756)
        job.log("Wrote "+outpath)

        return pd.Series({"Done":mol_sfam_id})

    #Read in H5 for entire table
    infpath = os.path.join(data_path_prefix, "IBIS_inferred.h5")
    try:
        inferred_interfaces = pd.read_hdf(unicode(infpath), table)
    except KeyError as e:
        job.log("Unable to load inferred interfaces from {}: {}".format(infpath, e))
        return

    #Add PDB, chain, and sdi information
    inferred_interfaces = pd.merge(inferred_interfaces, struct_domains, how="left", left_on="mol_sdi_id", right_on="sdi")

    #Add column to group by
    inferred_interfaces["superfamily"] = inferred_interfaces["nbr_superfam_id"]

    #Convert resn to PDB format
    ddf = dd.from_pandas(inferred_interfaces, name="inferred_interfaces", npartitions=NUM_WORKERS)

    #Concert MMDB to PDB and merge with observed, then save, but on groups in parallel
    ddf.groupby("superfamily").apply(apply_df, meta=pd.DataFrame({"Done":[str]}))\
        .compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)

def merge_inferred_interactome(job, dataset_name, cdd_name, cdd_id, cores=4):
    path = get_interfaces_path(dataset_name)
    if os.path.isfile(os.path.join(path, "{}.inferred_interactome".format(cdd_name))):
        return
    for table in xrange(1, 87):
        try:
            df = pd.read_hdf(os.path.join(path, "Intrac{}".format(table), "{}.inferred_interactome".format(cdd_id)), "table")
        except IOError:
            print "Can't open", os.path.join(path, "Intrac{}".format(table), "{}.inferred_interactome".format(cdd_id))
            continue
        df.to_hdf(unicode(os.path.join(path, "{}.inferred_interactome".format(cdd_name))), "table", table=True, mode='a', append=True, complevel=9, complib="bzip2", min_itemsize=768)
        del df

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
    for cdd, sfam_id in iter_cdd(use_id=True):
        try:
            job.addChildJobFn(merge_inferred_interactome, dataset_name, cdd, int(sfam_id))
        except ValueError:
            job.log("Error adding {} (Cannot convert {} to string)".format(cdd, sfam_id))

def start_toil_inferred(job, dataset_name):
    path = get_interfaces_path(dataset_name)
    for table in xrange(1, 87):
        table_dir = os.path.join(path, "Intrac{}".format(table))
        if not os.path.isdir(table_dir):
            os.makedirs(table_dir)
        job.addChildJobFn(get_inferred_structural_interactome_by_table, dataset_name, table, memory="120000M")
    job.addFollowOnJobFn(start_toil_inferred_merge, dataset_name, memory="120000M")

def start_toil(job, dataset_name):
    if not os.path.isdir(get_interfaces_path(dataset_name)):
        os.makedirs(get_interfaces_path(dataset_name))

    obsjob = job.addChildJobFn(start_toil_observed, dataset_name)
    infjob = obsjob.addFollowOnJobFn(start_toil_inferred, dataset_name)
    infmergejob = infjob.addFollowOnJobFn(start_toil_inferred, dataset_name)

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
