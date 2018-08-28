from dask_jobqueue import SLURMCluster
from dask.distributed import Client

import distributed.joblib
from joblib import Parallel, parallel_backend

from skopt import BayesSearchCV

from molmimic.torch_model.torch_train import Molmimic

cluster = SLURMCluster(
    walltime="3-00:00:00",
    memory="12000M",
    cores=16,
    project="muragroup",
    queue="gpu",
    gres="gpu:p100:1",
    ntasks="1"
)
cluster.adapt()
client = Client(cluster)

space = {
        "learning_rate": Real(1e-6, 1e-1, prior='log-uniform'),
        "num_epochs": Integer(30, 500, prior='log-uniform'),
        "batch_size": Integer(1, 30, prior='log-uniform'),
        "dropout_depth": Integer(0,1), #boolean
        "dropout_width": Integer(0,1), #boolean
        "dropout_p": Real(0., 1., prior='log-uniform')
}

opt = BayesSearchCV(Molmimic(), space, n_iter=100)

with parallel_backend('dask.distributed', client=client):
    opt.fit("default")

params = {"params": opt.best_params_,
          "best_score": opt.best_score_}

import json
with open("best_molmimic.json", "w") as f:
    json.dump(params, f)
