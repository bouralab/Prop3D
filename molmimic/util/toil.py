#From toil_rnaseq
#By John Vivian January 13 2018

import os
from math import ceil
import tarfile
import uuid
from contextlib import closing

import pandas as pd
from joblib import Parallel, delayed

from molmimic.util import safe_remove, safe_call
from molmimic.util.iostore import IOStore

from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

class RealtimeLogger_:
    @staticmethod
    def info(*args, **kwds):
        print(args, kwds)

def createToilOptions():
    parser = Job.Runner.getDefaultArgumentParser()

    if any("SLURM" in k for k in os.environ.keys()):
        tmpDir = os.path.join(os.getcwd(), "tmp-{}".format(uuid.uuid4()))
        if not os.path.isdir(tmpDir):
            os.makedirs(tmpDir)
        parser.set_defaults(workDir=tmpDir, statePollingWait=120,
            clean="always", logLevel="DEBUG", disableAutoDeployment=True)

    parser.set_defaults(clean="always", logLevel="DEBUG", targetTime=1)

    return parser

def partitions(l, partition_size):
    """
    >>> list(partitions([], 10))
    []
    >>> list(partitions([1,2,3,4,5], 1))
    [[1], [2], [3], [4], [5]]
    >>> list(partitions([1,2,3,4,5], 2))
    [[1, 2], [3, 4], [5]]
    >>> list(partitions([1,2,3,4,5], 5))
    [[1, 2, 3, 4, 5]]
    :param list l: List to be partitioned
    :param int partition_size: Size of partitions
    """
    for i in range(0, len(l), partition_size):
        yield l[i:i + partition_size]

def cleanup_ids(job, ids_to_delete):
    """
    Delete fileStoreIDs for files no longer needed

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param list ids_to_delete: list of FileStoreIDs to delete
    """
    [job.fileStore.deleteGlobalFile(x) for x in ids_to_delete if x is not None]


def map_job(job, func, inputs, *args, **kwds):
    """
    Spawns a tree of jobs to avoid overloading the number of jobs spawned by a single parent.
    This function is appropriate to use when batching samples greater than 1,000.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param function func: Function to spawn dynamically, passes one sample as first argument
    :param list inputs: Array of samples to be batched
    :param list args: any arguments to be passed to the function
    """
    # num_partitions isn't exposed as an argument in order to be transparent to the user.
    # The value for num_partitions is a tested value
    num_partitions = 100
    partition_size = int(ceil(len(inputs)/num_partitions))
    if partition_size > 1:
        RealtimeLogger.info("MAP_JOB: total: {}; paritions_size: {}".format(
            len(inputs), partition_size
        ))
        for partition in partitions(inputs, partition_size):
            job.addChildJobFn(map_job, func, partition, *args, **kwds)
    else:
        RealtimeLogger.info("MAP_JOB: Running: {}".format(len(inputs)))
        for sample in inputs:
            job.addChildJobFn(func, sample, *args, **kwds)

def map_job_follow_ons(job, func, inputs, *args, **kwds):
    # num_partitions isn't exposed as an argument in order to be transparent to the user.
    # The value for num_partitions is a tested value
    num_partitions = 100
    partition_size = int(ceil(len(inputs)/num_partitions))
    # if partition_size > 1:
    #     for partition in partitions(inputs, partition_size):
    #         job.addChildJobFn(map_job, func, partition, *args, **kwds)
    # else:
    #     for sample in inputs:
    #         job.addChildJobFn(func, sample, *args, **kwds)

    run_samples, follow_on_samples = inputs[:partition_size], inputs[partition_size:]
    RealtimeLogger.info("MAP_JOB: total: {}; paritions_size: {}; num_paritions: {}; num_run: {}".format(
        len(inputs), partition_size, len(follow_on_samples), len(run_samples)
    ))
    for sample in run_samples:
        job.addChildJobFn(func, sample, *args, **kwds)
    if partition_size > 1:
        job.addFollowOnJobFn(map_job, func, follow_on_samples, *args, **kwds)

def loop_job_rv(job, func, inputs, *args, cores=1, mem="2G", **kwds):
    if cores > 1:
        RealtimeLogger.info("Looping over domain in PARALLEL with {} cores".format(cores))
        for rv in Parallel(n_jobs=cores)(delayed(safe_call)(job, func, i, *args, **kwds) for i in inputs):
            yield rv
    else:
        RealtimeLogger.info("Looping over domain")

        for input in inputs:
            yield safe_call(job, func, input, *args, **kwds)

def loop_job(job, func, inputs, *args, cores=1, mem="2G", **kwds):
    [None for _ in loop_job_rv(job, func, inputs, *args, cores=cores, mem=mem, **kwds)]
    return

def map_job_rv(job, func, inputs, *args, **kwds):
    """
    Spawns a tree of jobs to avoid overloading the number of jobs spawned by a single parent.
    This function is appropriate to use when batching samples greater than 1,000.

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param function func: Function to spawn dynamically, passes one sample as first argument
    :param list inputs: Array of samples to be batched
    :param list args: any arguments to be passed to the function
    """
    # num_partitions isn't exposed as an argument in order to be transparent to the user.
    # The value for num_partitions is a tested value
    num_partitions = 100
    partition_size = int(ceil(len(inputs)/num_partitions))
    if partition_size > 1:
        promises = [job.addChildJobFn(map_job, func, partition, *args, **kwds).rv() \
            for partition in partitions(inputs, partition_size)]
    else:
        promises = [job.addChildJobFn(func, sample, *args, **kwds).rv() for sample in inputs]

    return promises

def map_job_to_pandas(job, func, inputs, *args, **kwds):
    # num_partitions isn't exposed as an argument in order to be transparent to the user.
    # The value for num_partitions is a tested value
    num_partitions = 100
    partition_size = int(ceil(len(inputs)/num_partitions))
    if partition_size > 1:
        promises = [job.addChildJobFn(map_job, func, partition, *args, **kwds).rv() \
            for partition in partitions(inputs, partition_size)]
    else:
        promises = [job.addChildJobFn(func, sample, *args, **kwds).rv() for sample in inputs]

    df = job.addFollowOnJobFn(pandas_concat, promises).rv()

    return df

def pandas_concat(job, df_promises):
    return pd.concat(df_promises, axis=0)

def map_job_rv_list(promises, *path):
    for p in promises:
        if p is None:
            continue
        rv = p.rv()
        if isinstance(rv, list) and len(rv) > 0:
            for p1 in map_job_rv_list(rv):
                if p1 is not None:
                    yield p1
        else:
            if p is not None:
                yield p

def finish_group(job, key, store):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get(store)

    done_file = os.path.join(work_dir, "DONE")
    done_key = "{}/DONE".format(key)

    with open(done_file, "w") as fh:
        pass

    store.write_output_file(done_file, done_key)

    safe_remove(done_file)


def consolidate_output(job, config, output):
    """
    Combines the contents of the outputs into one tarball and places in output directory or s3

    :param JobFunctionWrappingJob job: passed automatically by Toil
    :param Expando config: Dict-like object containing workflow options as attributes
    :param dict(str, str) output:
    """
    # Collect all tarballs from fileStore
    tars = {}
    for tool, filestore_id in list(output.items()):
        tars[os.path.join(config.uuid, tool)] = job.fileStore.readGlobalFile(filestore_id)

    # Consolidate tarballs into one output tar as streams (to avoid unnecessary decompression)
    out_tar = os.path.join(job.tempDir, config.uuid + '.tar.gz')
    with tarfile.open(out_tar, 'w:gz') as f_out:
        for name, tar in list(tars.items()):
            with tarfile.open(tar, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        tarinfo.name = os.path.join(name, os.path.basename(tarinfo.name))
                        f_out.addfile(tarinfo, fileobj=f_in_file)

    # Move to output location
    IOStore.get(congif.filestore).write_output_file(out_tar, config.uuid + '.tar.gz')
