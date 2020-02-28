def setup_dask(num_workers):
    dask.config.set(scheduler='multiprocessing')
    dask.config.set(pool=ThreadPool(num_workers))
