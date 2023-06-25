from numbers import Number

import h5pyd
import pandas as pd
from toil.realtimeLogger import RealtimeLogger

def split_dataset_at_level(job, cath_full_h5, superfamily, sfam_df, level_key, level_name,
  split_size={"train":0.8, "validation":0.1, "test":0.1}):
    """Split a dataset into train/validation/test sets, saving the splits into new h5 groups with
    links back to the the main dataset.
    
    Paramters:
    ----------
    job : toi.job.Job
        Toil job
    cath_full_h5 : str
        Path to H5 file on HSDS enpoint
    superfamily : str
        Group prefix, can be empty ('') for h5 file
    sfam_df : pd.DataFrame
        The data frame to split. Each row must a single protein and the df M\must contain 2 columns: 
            (i) "cath_domain", the column of protein domain names, must match groups of the same name in this 'superfamily' group;
            (ii) level_key, custom variable name for the name of the cluster the protein domain belongs to
    level_name : str
        Name of the column that contains cluster names
    split_size : Dict [split_name->split perecent]
        A dictionary containing the total number of splits
    """
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

        with h5pyd.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
            store.require_group(f"{superfamily}/data_splits/{level_name}")
            group = store.require_group(f"{superfamily}/data_splits/{level_name}/{split_name}")
            group.attrs["percent"] = size_pct

            for domain in subset:
                group[domain] = store[f'{superfamily}/domains/{domain}']