import os
import re
import sys
import urllib.request
import multiprocessing
from numbers import Number
from collections import defaultdict

import pandas as pd
from joblib import Parallel, delayed
from toil.realtimeLogger import RealtimeLogger

from molmimic.util.toil import map_job

try:
    import h5pyd as h5py
    DISTRIBUTED = True
except ImportError:
    try:
        import h5py
        DISTRIBUTED = False
    except:
        raise ImportError("h5pyd or h5py must be installed")

def split_superfamily_at_level(job, cath_full_h5, superfamily, sfam_df, level_key, level_name,
  split_size={"train":0.8, "validation":0.1, "test":0.1}):
    if isinstance(split_size, (list, tuple)):
        if len(split_size)==1:
            assert split_size[0]<1
            other_size = (1-split_size[0])/2
            split_size = dict(zip(["train", "validation", "test"],
                sorted([split_size[0], other_size, other_size])))
        else:
            assert sum(split_size) == 1
            split_size = {f"split{i}":s for i, s in enumerate(split_size)}
    elif isinstance(split_size, Number) and split_size<1:
        other_size = (1-split_size)/2
        split_size = {"train":split_size, "validation":other_size, "test":other_size}
    elif isinstance(split_size, dict) and sum(split_size.values()) == 1:
        #Correct
        pass
    else:
        raise RuntimeError("Invalid split_size. Must be dict {split_name:split_pct}, a list of split sizes (names automatically assinged), or a single number")

    start = 0
    subsets = []

    RealtimeLogger.info(f"Start splits for {superfamily} at level {level_name}")

    split_sizes = sorted(split_size.items(), key=lambda x: x[1], reverse=True)

    clusters = sfam_df.groupby(level_key)
    sorted_cluster_indices = list(sorted(clusters.indices.keys()))
    split_index_start = 0
    last_size = [None, None]

    for split_num, (split_name, split_pct) in enumerate(split_sizes):
        if split_num < len(split_sizes)-1:
            ideal_set1_size = int(clusters.ngroups*split_pct)

            while True:
                set1_clusters = sorted_cluster_indices[split_index_start:split_index_start+ideal_set1_size]
                set1 = [idx for cluster in set1_clusters for idx in clusters.get_group(cluster).index] #.indices[cluster]]
                size_pct = len(set1)/(len(sfam_df))
                print("size", len(set1), len(sfam_df), size_pct, ideal_set1_size)
                if size_pct in last_size:
                    break
                if size_pct > split_pct+.01:
                    ideal_set1_size -= 1
                elif size_pct < split_pct-.01:
                    ideal_set1_size += 1
                else:
                    break

                last_size[0] = last_size[1]
                last_size[1] = size_pct

            subset_idx = list(sorted(set1))
            subset = sfam_df[sfam_df.index.isin(set1)]["cath_domain"]

            #Reset index for next iteration to skip current domains
            split_index_start += ideal_set1_size
        else:
            set1_clusters = sorted_cluster_indices[split_index_start:]
            set1 = [idx for cluster in set1_clusters for idx in clusters.get_group(cluster).index] #.indices[cluster]]
            size_pct = len(set1)/(len(sfam_df))
            subset = sfam_df[sfam_df.index.isin(set1)]["cath_domain"]

        with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
            store.require_group(f"{superfamily}/data_splits/{level_name}")
            group = store.require_group(f"{superfamily}/data_splits/{level_name}/{split_name}")
            group.attrs["percent"] = size_pct

            RealtimeLogger.info(f'subset {subset}')

            for domain in subset:
                RealtimeLogger.info(f'{superfamily}/domains/{domain}')
                group[domain] = store[f'{superfamily}/domains/{domain}']


def create_splits_for_superfamily_levels(job, sfam, cath_full_h5):
    superfamily, sfam_df = sfam
    RealtimeLogger.info(f"Start splits for {superfamily}")

    with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
        store.require_group(f"{superfamily}/data_splits")

    for level_key, level_name in [("S", "S35"), (list("SO"), "S60"), (list("SOL"), "S95"), (list("SOLI"), "S100")]:
        job.addChildJobFn(split_superfamily_at_level, cath_full_h5, superfamily, sfam_df, level_key, level_name)

def create_representatives_for_superfamily(job, sfam, cath_full_h5):
    from molmimic.parsers.cath import CATHApi
    cath = CATHApi()
    hierarchy = cath.list_children_in_heirarchy(sfam, 5)
    representatives = [child["example_domain_id"] for child in hierarchy["children"]]

    key = f"{sfam.replace('.', '/')}/representatives"

    missing_domains = []

    with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
        group = store.require_group(key)

        for domain in representatives:
            RealtimeLogger.info(f"Adding {domain} from {sfam.replace('.', '/')}/domains/{domain} to {key}/{domain}")
            try:
                group[domain] = store[f"{sfam.replace('.', '/')}/domains/{domain}'"]
            except KeyError:
                missing_domains.append(domain)

        if len(missing_domains) > 0:
            store[key].attrs["missing_domains"] = missing_domains
        store[key].attrs["total_domains"] = len(representatives)

def create_splits(job, cath_full_h5, all_superfamilies):
    RealtimeLogger.info(f"Start all splits {all_superfamilies}")
    sfams = [g for g in all_superfamilies.groupby("h5_key")]
    RealtimeLogger.info(f"sfam splits {sfams}")
    map_job(job, create_splits_for_superfamily_levels, sfams, cath_full_h5)
    map_job(job, create_representatives_for_superfamily, [s.iloc[0].cath_code for _, s in sfams], cath_full_h5)
    job.addFollowOnJobFn(finish_section, cath_full_h5, "completed_domain_splits")
    job.addFollowOnJobFn(finish_section, cath_full_h5, "completed_representatives")

def process_cath_domain_list_for_group(job, group, cath_full_h5):
    name, group_df = group
    RealtimeLogger.info(f"process_cath_domain_list_for_group {name}")

    with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
        for _, row in group_df.iterrows():
            group = store.require_group(f"{row.h5_key}/domains/{row.cath_domain}")
            group.domain_length = row.domain_length
            group.resolution = row.resolution

def process_cath_domain_list(job, cath_full_h5, cathcode=None, skip_cathcode=None, force=False, work_dir=None):
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    run_domain_names = True
    run_splits = True

    if isinstance(force, bool) or (is_num(force) and int(force)<3):
        try:
            with h5py.File(cath_full_h5, mode="r", use_cache=False) as store:
                run_domain_names = not store.attrs.get("completed_domain_list", False)
                run_splits = not store.attrs.get("completed_domain_splits", False)
        except IOError:
            raise

    if not run_splits:
        #Already exists do not run again
        return

    #Run splits needs same info, so build pandas frame for both:

    cath_domain_list_file = os.path.join(work_dir, "cath-domain-list.txt")
    if not os.path.isfile(cath_domain_list_file):
        urllib.request.urlretrieve(
            "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt",
            cath_domain_list_file)

    names = pd.read_csv(cath_domain_list_file, delim_whitespace=True, header=None, comment="#",
        names=["cath_domain", *list("CATHSOLID"), "domain_length", "resolution"])
    names = names.assign(cath_code=names["C"].astype(str)+"."+names["A"].astype(str)+"."+names["T"].astype(str)+"."+names["H"].astype(str))
    names = names.assign(h5_key="/"+names["cath_code"].str.replace(".","/"))
    names = names.assign(group=names["cath_code"].str.split(".", expand=True)[[0,1]].fillna("").agg('.'.join, axis=1))

    RealtimeLogger.info(names)

    if cathcode is not None:
        if not isinstance(cathcode, (list, tuple)):
            cathcode = [cathcode]

        use_sfams = tuple([".".join(sfam)+"." if isinstance(sfam, (list,tuple)) else sfam.replace("/", ".")+"." \
            for sfam in cathcode])

        names = names[(names["cath_code"]+".").str.startswith(use_sfams)]

    if skip_cathcode is not None:
        if not isinstance(skip_cathcode, (list, tuple)):
            skip_cathcode = [skip_cathcode]

        skip_sfams = tuple([".".join(sfam)+"." if isinstance(sfam, (list,tuple)) else sfam.replace("/", ".")+"." \
            for sfam in skip_cathcode])

        names = names[~(names["cath_code"]+".").str.startswith(skip_sfams)]

    if run_domain_names:
        groups = [g for g in names.groupby("group")]
        map_job(job, process_cath_domain_list_for_group, groups, cath_full_h5)

        job.addFollowOnJobFn(finish_section, cath_full_h5, "completed_domain_list")

    #Create splits
    all_superfamilies = names[names['cath_code'].str.split('.').agg(len)==4]
    job.addFollowOnJobFn(create_splits, cath_full_h5, all_superfamilies)

def process_cath_names_for_group(job, group, cath_full_h5):
    name, group_df = group
    RealtimeLogger.info(f"process_cath_names_for_group {name}")

    with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
        for _, row in group_df.iterrows():
            group = store.require_group(row.h5_key)
            group.description = row.description
            group.representativeDomain = row.representative
            if row.cath_code.count(".") == 3:
                store.require_group(f"{row.h5_key}/domains")
                RealtimeLogger.info(f"Create domains for {row.h5_key}")
        store.flush()

def delete_groups(root):
    if hasattr(root, "keys"):
        for key in root.keys():
            delete_groups(root[key])
            del root[key]
    del root

def process_cath_names(job, cath_full_h5, cathcode=None, skip_cathcode=None, force=False, work_dir=None):
    """Will overwrite all files"""
    RealtimeLogger.info("process_cath_names")
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    RealtimeLogger.info(f"Creating file {cath_full_h5}")

    if is_num(force) and int(force)==3:
        RealtimeLogger.info(f"Removing all previous data from {cath_full_h5}")
        try:
            #Empty file if it has been created before
            with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
                delete_groups(store)
                store.flush()
                RealtimeLogger.info(f"Deleted {cath_full_h5}")

            with h5py.Folder(os.path.dirname(cath_full_h5)+"/", mode="a") as hparent:
                del hparent[os.path.basename(cath_full_h5)]
        except IOError:
            raise
            #not created
            pass
        except OSError:
            pass
        RealtimeLogger.info(f"Removed all previous data from {cath_full_h5}")
    else:
        RealtimeLogger.info(f"Not deleting any previous data from {cath_full_h5}")

    cath_names_file = os.path.join(work_dir, "cath-names.txt")
    if not os.path.isfile(cath_names_file):
        urllib.request.urlretrieve(
            "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-names.txt",
            cath_names_file)

    names = pd.read_csv(cath_names_file, sep="    ", header=None, comment="#",
        names=["cath_code", "representative", "description"])
    names["description"] = names["description"].str[1:]
    names = names.assign(h5_key="/"+names["cath_code"].str.replace(".","/"))
    names = names.assign(group=names["cath_code"].str.split(".", expand=True)[[0,1]].fillna("").agg('.'.join, axis=1))

    RealtimeLogger.info(f"Read cath names file")

    if cathcode is not None:
        if not isinstance(cathcode, (list, tuple)):
            cathcode = [cathcode]

        use_sfams = tuple([".".join(sfam)+"." if isinstance(sfam, (list,tuple)) else sfam.replace("/", ".")+"." \
            for sfam in cathcode])

        names_ = None
        for code in use_sfams:
            subset = names.apply(lambda r: r if code.startswith(r["cath_code"]+".") else pd.Series(dtype=str), axis=1).dropna()
            if names_ is None:
                names_ = subset
            else:
                names_ = pd.concat((names_, subset))
        names = names_.drop_duplicates().reset_index(drop=True)
        del names_

    if skip_cathcode is not None:
        if not isinstance(skip_cathcode, (list, tuple)):
            skip_cathcode = [skip_cathcode]

        skip_sfams = tuple([".".join(sfam)+"." if isinstance(sfam, (list,tuple)) else sfam.replace("/", ".")+"." \
            for sfam in skip_cathcode])

        skip_names = None
        for code in use_sfams:
            subset = names.apply(lambda r: r if code.startswith(r["cath_code"]+".") else pd.Series(), axis=1).dropna()
            if skip_names is None:
                skip_names = subset
            else:
                skip_names = pd.concat((skip_names, subset))
        skip_names = skip_names.drop_duplicates().reset_index(drop=True)
        skip_names = pd.merge(names, skip_names, how='inner', on="cath_code").index

        names = names.drop(skip_names)

    RealtimeLogger.info(f"Names Running {len(names[names['cath_code'].str.split('.').agg(len)==4])} and {len(names)} nodes")

    root_nodes = names[names["cath_code"].str.split(".").agg(len)==1]
    process_cath_names_for_group(job, ("root", root_nodes), cath_full_h5)

    groups = [g for g in names[~names.index.isin(root_nodes.index)].groupby("group")]
    map_job(job, process_cath_names_for_group, groups, cath_full_h5)

    job.addFollowOnJobFn(finish_section, cath_full_h5, "completed_names")

def setup_custom_file(cath_full_h5, pdbs, force=False):
    RealtimeLogger.info(f"Creating file {cath_full_h5}")

    if is_num(force) and int(force)==3:
        RealtimeLogger.info(f"Removing all previous data from {cath_full_h5}")
        try:
            #Empty file if it has been created before
            with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
                delete_groups(store)
                store.flush()
                RealtimeLogger.info(f"Deleted {cath_full_h5}")

            with h5py.Folder(os.path.dirname(cath_full_h5)+"/", mode="a") as hparent:
                del hparent[os.path.basename(cath_full_h5)]
        except IOError:
            raise
            #not created
            pass
        except OSError:
            pass
        RealtimeLogger.info(f"Removed all previous data from {cath_full_h5}")
    else:
        RealtimeLogger.info(f"Not deleting any previous data from {cath_full_h5}")

    with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
        pass


def create_h5_hierarchy(job, cath_full_h5, cathcode=None, skip_cathcode=None, pdbs=None, work_dir=None, force=False):
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    if pdbs is not None:
        return setup_custom_file(cath_full_h5, pdbs, force=force)

    run_names = True

    if isinstance(force, bool) or (is_num(force) and int(force)<3):
        try:
            with h5py.File(cath_full_h5, mode="r", use_cache=False) as store:
                run_names = not store.attrs.get("completed_names", False)
        except IOError:
            #Never created, ignore
            pass

    RealtimeLogger.info(f"Force is {force}, run_names={run_names}")

    if run_names:
        job.addChildJobFn(process_cath_names, cath_full_h5, cathcode=cathcode, skip_cathcode=skip_cathcode, force=force)

    job.addFollowOnJobFn(process_cath_domain_list, cath_full_h5, cathcode=cathcode, skip_cathcode=skip_cathcode, force=force)

def finish_section(job, cath_full_h5, attribute):
    with h5py.File(cath_full_h5, mode="a", use_cache=False) as store:
        store.attrs[attribute] = True

def is_num(a):
    try:
        int(a)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    force = len(sys.argv)>1 and args[1] in ["-f", "--force"]
