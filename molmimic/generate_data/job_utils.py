#From toil_rnaseq
#By John Vivian January 13 2018

import os
from math import ceil
import tarfile
from contextlib import closing
from molmimic.generate_data.iostore import IOStore

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
        for partition in partitions(inputs, partition_size):
            job.addChildJobFn(map_job, func, partition, *args, **kwds)
    else:
        for sample in inputs:
            job.addChildJobFn(func, sample, *args, **kwds)

def map_job_rv(job, func, inputs, *args):
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
    partition_size = len(inputs) / num_partitions
    if partition_size > 1:
        promises = [job.addChildJobFn(map_job, func, partition, *args).rv() \
            for partition in partitions(inputs, partition_size)]
    else:
        promises = [job.addChildJobFn(func, sample, *args).rv() for sample in inputs]

    return promises

def map_job_rv_list(promises, *path):
    for p in promises:
        rv = p.rv()
        if isinstance(rv, list) and len(rv) > 0:
            for p1 in map_job_rv_list(rv):
                yield p1
        else:
            yield p


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
