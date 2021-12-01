import os
import sys

import toil
import pandas as pd

from toil.realtimeLogger import RealtimeLogger

def get_file(job, prefix, path_or_fileStoreID, work_dir=None, cache=False, return_type=False):
    if isinstance(path_or_fileStoreID, str) and os.path.isfile(path_or_fileStoreID):
        if return_type:
            return path_or_fileStoreID, "path"
        else:
            return path_or_fileStoreID
    else:
        work_dir = work_dir or job.fileStore.getLocalTempDir()
        new_file = os.path.join(work_dir, prefix)

        if isinstance(path_or_fileStoreID, (toil.fileStores.FileID, str)):
            if cache:
                new_file = job.fileStore.readGlobalFile(
                    path_or_fileStoreID, userPath=new_file, cache=True)
            else:
                with job.fileStore.readGlobalFileStream(path_or_fileStoreID) as fs, open(new_file, "wb") as nf:
                    for line in fs:
                        nf.write(line)
        elif hasattr(path_or_fileStoreID, "read_input_file"):
            #Might be file store itself
            path_or_fileStoreID.read_input_file(prefix, new_file)

        else:
            raise RuntimeError("Invalid path_or_fileStoreID {} {}".format(type(path_or_fileStoreID), path_or_fileStoreID))

        if return_type:
            return new_file, "fileStoreID"
        else:
            return new_file

def filter_hdf(hdf_path, dataset, column=None, value=None, columns=None,
  drop_duplicates=False, **query):
    #assert len(query) > 0 or (column is not None and value is not None)
    if len(query) == 0 and (column is not None and value is not None):
        if not isinstance(column, (list, tuple)) or not isinstance(value, (list, tuple)):
            where = "{}={}".format(column, value)
        elif len(column) == len(value):
            where = ["{}={}".format(c, v) for c, v in zip(column, value)]
        else:
            raise RuntimeError("Cols and values must match")
    elif len(query) > 0:
        where = ["{}={}".format(c,v) for c, v in list(query.items())]
    else:
        where = None

    try:
        df = pd.read_hdf(str(hdf_path), dataset, where=where, columns=columns, mode="r")
        if df.shape[0] == 0: raise KeyError
        if drop_duplicates:
            df = df.drop_duplicates()
    except (KeyError, ValueError, SyntaxError, OSError):
        df = filter_hdf_chunks(hdf_path, dataset, column=column, value=value,
            columns=columns, drop_duplicates=drop_duplicates, **query)
    return df

def filter_hdf_chunks(hdf_path, dataset, column=None, value=None, columns=None,
  chunksize=500, drop_duplicates=False, **query):
    df = None
    store = pd.HDFStore(str(hdf_path), mode="r")
    nrows = store.get_storer(dataset).nrows

    #for _df in pd.read_hdf(str(hdf_path), dataset, chunksize=chunksize):
    for i in range(nrows//chunksize + 1):
        _df = store.select(dataset, start=i*chunksize, stop=(i+1)*chunksize)

        if len(query) > 0:
            try:
                filtered_df = _df.query("("+") & (".join(["{}=={}".format(c, v) for c, v in list(query.items())])+")")
            except SyntaxError:
                expression = None
                for c, v in list(query.items()):
                    _exp = _df[c]==v
                    if expression is None:
                        expression = _exp
                    else:
                        expression &= _exp
                filtered_df = _df[expression]
        elif column is not None and value is not None:
            if not isinstance(column, (list, tuple)) and not isinstance(value, (list, tuple)):
                filtered_df = _df[_df[column]==value].copy()
            elif len(column) == len(value):
                filtered_df = _df.query(" & ".join(["{}=={}".format(c, v) for c, v in zip(column, value)]))
            else:
                raise RuntimeError("Cols and values must match")
        else:
            filtered_df = _df.copy()
        if filtered_df.shape[0]>0:
            if columns is not None:
                filtered_df = filtered_df[columns]
            if drop_duplicates:
                filtered_df = filtered_df.drop_duplicates()
            df = pd.concat((df, filtered_df), axis=0) if df is not None else filtered_df
            if drop_duplicates:
                df = df.drop_duplicates()
            del _df
            _df = None
    if df is None:
        raise TypeError("Unable to parse HDF")

    #Make sure store is closed
    store.close()

    return df

def make_h5_tables(files, iostore):
    for f in files:
        iostore.read_input_file(f, f)
        store = pd.HDFStore(str(f))
        for key in list(store.keys()):
            df = store.get(key)
            df.to_hdf(str(f+".new"), "table", format="table",
                table=True, complevel=9, complib="bzip2", min_itemsize=1024)
        iostore.write_output_file(f+".new", f)

def split_df(df, cluster_key, split_size=0.8):
    clusters = df.groupby(cluster_key)
    sorted_cluster_indices = list(sorted(clusters.indices.keys()))

    ideal_set1_size = int(clusters.ngroups*split_size)

    last_size = [None, None]
    while True:
        set1_clusters = sorted_cluster_indices[:ideal_set1_size]
        set1 = [idx for cluster in set1_clusters for idx in clusters.get_group(cluster).index] #.indices[cluster]]
        size_pct = len(set1)/(len(df))
        print("size", len(set1), len(df), size_pct, ideal_set1_size)
        if size_pct in last_size:
            break
        if size_pct > split_size+.01:
            ideal_set1_size -= 1
        elif size_pct < split_size-.01:
            ideal_set1_size += 1
        else:
            break

        last_size[0] = last_size[1]
        last_size[1] = size_pct

    set1 = list(sorted(set1))
    set1 = df.index.isin(set1)
    set2 = ~set1

    df1, df2 = df.loc[set1], df.loc[set2]
    return df1, df2
