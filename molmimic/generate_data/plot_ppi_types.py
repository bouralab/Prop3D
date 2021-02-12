import os
import glob

import pandas as pd
from joblib import Parallel, delayed
from molmimic.generate_data.util import iter_cdd

import dask.dataframe as dd

NUM_WORKERS = 20

def get_counts(dataset_name, sfam_id):
    path = os.path.join(get_interfaces_path(dataset_name), "by_superfamily",
        str(int(sfam_id)), "{}_bsa.h5".format(int(sfam_id)))

    try:
        store = pd.HDFStore(path)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        return

    if len(list(store.keys())) == 0:
        store.close()
        return

    obs = store.get("/observed")
    obs = obs[(~obs["ppi_type"].isnull())&(obs["ppi_type"]!="unknown")]
    obs_sizes = obs.groupby("ppi_type").size().reset_index(name='count').copy()
    if obs_sizes.shape[0]>0:
        obs_sizes.loc[:, "status"] = "observed"
    del obs

    if not "/inferred" in list(store.keys()):
        store.close()
        return obs_sizes.set_index("ppi_type")

    inf = store.get("/inferred")
    inf = inf[(~inf["ppi_type"].isnull())&(inf["ppi_type"]!="unknown")]
    inf_sizes = inf.groupby("ppi_type").size().reset_index(name='count').copy()
    if inf_sizes.shape[0]>0:
        inf_sizes.loc[:, "status"] = "inferred"

    sizes = pd.concat([obs_sizes, inf_sizes], axis=0)

    changed = inf.groupby(["ppi_type", "ppi_type_obs"]).size().reset_index(name='count').copy()
    #changed = changed[changed["ppi_type"]!=changed["ppi_type_obs"]]
    #import pdb; pdb.set_trace()
    changed = changed.apply(lambda r: pd.Series({
        "ppi_type":"{}->{}".format(r.ppi_type_obs, r.ppi_type),
        "count":r["count"]}), axis=1)
    if changed.shape[0]>0:
        changed.loc[:, "status"] = "inferred"
    del inf
    store.close()
    del store

    sizes = pd.concat([sizes, changed], axis=0).reset_index(drop=True)
    return sizes

def plot_ppi_types(dataset_name, superfams=None, histograms=True, violin=False):
    if superfams is not None:
        return

    if histograms:
        _sizes = Parallel(n_jobs=NUM_WORKERS)(delayed(get_counts)(dataset_name, sfam_id) \
            for _, sfam_id in iter_cdd(use_id=True, group_superfam=True))
        return _sizes
        sizes = None
        for _size in _sizes[1:]:
            if _size is not None:
                if sizes is not None:
                    sizes = sizes.radd(_size.fillna(0.0), fill_value=0.0)
                else:
                    sizes = _size.fillna(0.0)

        return sizes.fillna(0.0)

    if violin:
        files = glob.glob(os.path.join(get_interfaces_path(dataset_name), "by_superfamily",
            "*", "*_bsa.h5"))
        obs_df = dd.read_hdf(files, "/observed")
        obs_df = obs_df.repartition(nparititon=20)
        obs_df = obs_df[["ppi_type", "bsa"]]
        obs_df = obs_df.assign({"status":"observed"})

        inf_df = dd.read_hdf(files, "/inferred")
        inf_df = inf_df.repartition(nparititon=20)
        inf_df = inf_df[["ppi_type", "bsa"]]
        inf_df = inf_df.assign({"status":"inferred"})

        df = obs_df.append(inf_df)
