import os
import pandas as pd
from joblib import Parallel, delayed
from util import get_interfaces_path, iter_cdd

NUM_WORKERS = 20

def get_counts(sfam_id):
    path = os.path.join(get_interfaces_path(dataset_name), "by_superfamily",
        str(int(sfam_id)), "{}_bsa.h5".format(int(sfam_id)))
    store = pd.HDFStore(path)

    obs = store.get("/observed")
    obs_sizes = obs.groupby("ppi_type").size().reset_index(name='observed').T
    del obs

    inf = store.get("/inferred")
    inf_sizes = inf.groupby("ppi_type").size().reset_index(name='inferred').T

    sizes = pd.concat([obs_size, inf_sizes], axis=1)

    changed = inf.groupby(["ppi_type", "ppi_type_obs"]).size().reset_index(name='count')
    #changed = changed[changed["ppi_type"]!=changed["ppi_type_obs"]]
    changed = changed.apply(lambda r: pd.Series({
        "{}->{}".format(r.ppi_type_obs, r.ppi_type):r.count}, index="inferred"), axis=1)
    del inf
    store.close()
    del store

    import pdb; pdb.set_trace()
    sizes = pd.concat([sizes, changed], axis=1)

def plot_ppi_types(superfams=None):
    if superfams is not None:
        return

    sizes = Parallel(n_jobs=NUM_WORKERS)(delayed(get_counts)(sfam_id) \
        for _, sfam_id in iter_cdd(use_id=True, group_superfam=True))
