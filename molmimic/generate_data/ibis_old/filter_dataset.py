import os, sys
import shutil
import glob

import numpy as np
import pandas as pd
import dask.dataframe as dd

from molmimic.generate_data.util import natural_keys, iter_cdd
from molmimic.generate_data.map_residues import compare_to_pdb

NUM_WORKERS = 20

def check_interface(pdb, chain, domNo, sdi, resi):
    pdb, chain, domNo, sdi, = row["mol_pdb"], row["mol_chain"], row["mol_domNo"], row["mol_sdi"]
    pdb_id = "{}_{}_sdi{}_{}".format(pdb, chain, sdi, domNum)
    features_path = os.path.join(get_features_path(dataset_name), row["mol_cdd"],
        pdb[1:3].lower(), "{}.npy".format(pdb_id))
    structure_path = os.path.join(get_structures_path(dataset_name), row["mol_cdd"],
        pdb[1:3].lower(),  "{}.pdb".format(pdb_id))

    if os.path.isfile(feature_path) and os.path.isfile(structure_path):
        #Remove residues not in structure
        return ",".join(("".join(r).strip() for r in compare_to_pdb(structure_path, mol_resi)))
    else:
        return np.nan

def binding_site_align(binding_site, aln, seq_len):
    binding_site = iter(binding_site)
    aln = iter(aln)

    bs_res = next(binding_site)
    aln_col = next(aln)

    seq_count = 0

    new_binding_site = []
    while seq_count < seq_len:
        rep = ""
        while aln_col.isdigit():
            rep += aln_col
            next(aln_col)
        try:
            rep = int(rep)
        except ValueError:
            rep = 1

        if seq_count+1 <= natural_keys(bs_res)[1] < seq_count+rep+1 and aln_col=="M":
            new_binding_site.append(bs_res)

        seq_count += rep
        bs_res = next(binding_site)
        aln_col = next(aln)

    return set(new_binding_site)

def process_cluster(df):
    centroid_pdb = df.name
    rep_pdb, rep_chain, rep_sdi, rep_domNo, _ = centroid_pdb
    rep_sdi, rep_domNo = rep_sdi[3:], rep_domNo[1:]
    binding_sites = binding_sites.copy()
    resi = binding_sites.apply(lambda row: binding_site_align(row["resi"].split(","), row["alignment"], row["length"]))
    binding_site = set()
    for r in resi:
        binding_site &= binding_site

    centroid = df[(df["int_pdb"]==rep_pdb)&(df["int_chain"]==rep_chain)&(df["int_sdi"]==rep_sdi)&(df["int_domNo"]==rep_domNo)].iloc[0].copy()
    centroid["int_resi"] = binding_site
    return centroid

def cluster_interactome(job, df, cdd, sfam_id=None, table="observed"):
    cdd = cdd.replace("/", "")
    if df is None and sfam_id is not None:
        int_path = os.path.join(get_interfaces_path(dataset_name), "{}.h5".format(cdd))
        df = pd.read_hdf(str(int_path), table)
    df["pdb_id"] = df.apply(lambda row: "{p}_{c}_sdi{s}_d{d}_{c}".format(
        p=row["int_pdb"], c=row["int_chain"], s=row["int_sdi"], d=row["int_domNo"]), axis=1)

    pdb_path = os.path.join(data_path_prefix, "structures", cdd, "{}_clusters.h5".format(cdd))
    pdb_clusters = pd.read_hdf(str(pdb_path), "table")

    cluster_int = pd.merge(df, pdb_clusters, how="left", left_on="pdb_id", right_on="label_query")
    cluster_int = cluster_int.groupby("label_target").apply(process_cluster)
    return cluster_int

def process_sdi(df):
    df = df.copy().reset_index()
    row = df.iloc[0]
    pdb, chain, domNo, sdi = row["mol_pdb"], row["mol_chain"], row["mol_domNo"], row["mol_sdi"]
    resi = sorted(set(",".join(df["mol_resi"]).split(",")), key=natural_keys)
    return pd.Series({
        "pdb": pdb,
        "chain": chain,
        "domNo": domNo,
        "sdi": sdi,
        "resi": check_interface(resi)
    })

def combine(df):
    """Aggregate function to combine dataframes that have already been verfied

    Parameters
    ----------
    df : pd.DataFrame

    Return
    ------
    Series with same cols as df, but resi have been merged
    """
    df = df.copy().reset_index()
    row = df.iloc[0]
    resi = sorted(set(",".join(df["resi"]).split(",")), key=natural_keys)
    return pd.Series({
        "pdb": row["pdb"],
        "chain": row["chain"],
        "domNo": row["domNo"],
        "sdi": row["sdi"],
        "resi": resi
    })

def filter_interfaces(job, interface, status, ppi_type):
    try:
        df = pd.read_hdf(str("{}_bsa.h5".format(interface)), status)
    except (IOError, KeyError):
        return

    if ppi_type != "all":
        df = df[df["ppi_type"]==ppi_type]

    ddf = dd.from_pandas(df, name="obs", npartitions=NUM_WORKERS)

    #Find interfaces the have structures and features and the residues found in structure
    meta = pd.DataFrame({"pdb":[str], "chain":[str], "domNo":[int], "sdi":[int], "resi":[str]})
    filtered_interfaces = ddf.groupby("sdi", as_index=False).apply(process_sdi, \
        meta=meta, axis=1).compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)

    #Remove all interfaces that do not match criterion
    filtered_interfaces = filtered_interfaces.dropna()
    filtered_interfaces["resi"] = filtered_interfaces["resi"].astype(str)

    filtered_interfaces.to_hdf(str(interface+"{}.{}.h5".format(interface, status)), ppi_type, complevel=9, complib="bzip2")
    clust_int = cluster_interactome(filtered_interfaces, cdd)
    clust_int.to_hdf(str(interface+".{}.clustered.h5".format(status)), ppi_type, complevel=9, complib="bzip2")
    del filtered_interfaces
    del clust_int
    del ddf
    del meta
    del df

def filter_mixed_interfaces(job, interface, ppi_type):
    try:
        obs = dd.read_hdf(str(interface+".observed.h5"), ppi_type)
        inf = dd.read_hdf(str(interface+".inferred.h5"), ppi_type)
    except (IOError, KeyError):
        return

    mixed = obs.append(inf).groupby("sdi", as_index=False).apply(combine, meta=meta).\
        compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    mixed.to_hdf(str(interface+".mixed.h5"), ppi_type, complevel=9, complib="bzip2")
    clust_int = cluster_interactome(mixed, cdd)
    clust_int.to_hdf(str(interface+".mixed.clustered.h5"), ppi_type, complevel=9, complib="bzip2")
    del mixed
    del obs

def filter_data(job, dataset_name, cdd, cores=NUM_WORKERS):
    interface = os.path.join(get_interfaces_path(dataset_name), cdd, cdd)
    # try:
    #     obs_df = pd.read_hdf(unicode(interface+".observed_bsa"), "table")
    # except (IOError, KeyError):
    #     return
    #
    # obs_ddf = dd.from_pandas(obs_df, name="obs", npartitions=NUM_WORKERS)
    #
    # #Find interfaces the have structures and features and the residues found in structure
    # meta = pd.DataFrame({"pdb":[str], "chain":[str], "domNo":[int], "sdi":[int], "resi":[str]})
    # filtered_interfaces = obs_ddf.groupby("sdi", as_index=False).apply(process_sdi, \
    #     axis=1).compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)
    #
    # #Remove all interfaces that do not match criterion
    # filtered_interfaces = filtered_interfaces.dropna()
    # filtered_interfaces["resi"] = filtered_interfaces["resi"].astype(str)

    #Save interfaces
    #filtered_interfaces.to_hdf(unicode("{}.h5".format(interface)), "observed")
    #obs_clustered = cluster_interactome(obs_df)
    #obs_clustered

    statuses = ["observed", "inferred"]
    ppi_types = ["all", "strong", "weak_transient", "transient"]

    for status in statuses:
        for ppi_type in ppi_types:
            job.addChildJobFn(filter_interfaces, interface, status, ppi_type)

    for ppi_type in ppi_types:
        job.addFollowOnJobFn(filter_mixed_interfaces, interface, ppi_type)

def start_toil(job, dataset_name):
    for cdd in iter_cdd():
        cdd = cdd.replace("/", "")
        job.addChildJobFn(filter_data, dataset_name, cdd)

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

    # if len(sys.argv) == 3:
    #     run_filter(sys.argv[1], sys.argv[2])
    # elif len(sys.argv) == 4 and sys.argv[1] == "filter":
    #     filter_ibis(sys.argv[2], sys.argv[3])

# def run_filter(dataset_name, ibis_data, job_name="filter_ibis", dependency=None):
#     if os.path.isdir(ibis_data):
#         ibis_data_files = glob.glob(os.path.join(ibis_data, "*.tsv"))
#     else:
#         ibis_data_files = [ibis_data]
#
#     job = SwarmJob(job_name, mem="1", individual=len(ibis_data_files)==1)
#     for ibis_data in ibis_data_files:
#         job += "python {} filter {}\n".format(__file__, dataset_name, ibis_data)
#
#     if len(ibis_data_files)==1:
#         jid = job.submit_individual(dependency=dependency)
#     else:
#         jid = job.run(dependency=dependency)
#
#     features_path = os.path.join(get_features_path(dataset_name), "done")
#     endjob = SwarmJob(job_name, mem="1", individual=True)
#     endjob += "touch {}\n".format(features_path)
#     endjob.submit_individual(dependency="afterany:"+jid)
#
#     while not os.path.isfile(features_path):
#         time.sleep(800)
#
# def submit_filter(dataset_name, ibis_data, job_name="filter_ibis", dependency=None):
#     job = SwarmJob(job_name, walltime="2-00:00:00", mem="1", individual=True)
#     job += "python {} {} {}\n".format(__file__, dataset_name, ibis_data)
#     job_id = job.submit_individual(dependency=dependency)
#     return job_id
