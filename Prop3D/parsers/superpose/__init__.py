import os
import re

import pandas as pd

from Prop3D.parsers.container import Container
from Prop3D.util.toil import map_job

from toil.realtimeLogger import RealtimeLogger

def start_one_vs_all_jobs(job, superposer, file_list, table_out_file=None, work_dir=None, cores=1, mem="72G", memory="72G", **kwds):
    RealtimeLogger.info("start")

    if job is None:
        from toil.job import Job
        job = Job()

    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    if not isinstance(file_list, (list, tuple)):
        file_list_file = get_file(job, "file_list", file_list, work_dir=work_dir, cache=True)
        with open(file_list_file) as fh:
            pdbs = [pdb.rstrip() for pdb in fh]
    else:
        pdbs = file_list

    results_file = lambda pdb: os.path.join(work_dir, f"{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")
    pdbs_to_run = [(i, pdb) for i, pdb in enumerate(pdbs) if not os.path.isfile(
        os.path.join(work_dir, results_file(pdb))) and i<len(pdbs)-1]

    RealtimeLogger.info("Starting distributed start_one_vs_all_jobs 222 parallel")
    map_job(job, one_vs_all, pdbs_to_run, pdbs, superposer, work_dir=work_dir,
        cores=cores, memory=mem, **kwds)

    return job.addFollowOnJobFn(combine_all_vs_all, superposer, file_list, table_out_file=table_out_file,
        work_dir=work_dir).rv()

def one_vs_all(job, exp_pdb, file_list, superposer, *args, work_dir=None, **kwds):
    if job is None:
        from toil.job import Job
        job = Job()

    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    if not isinstance(file_list, (list, tuple)):
        file_list_file = get_file(job, "file_list", file_list, work_dir=work_dir, cache=True)
        with open(file_list_file) as fh:
            all_pdbs = [pdb.rstrip() for pdb in fh]
    else:
        all_pdbs = file_list

    in_all_vs_all = False
    if isinstance(exp_pdb, (list, tuple)):
        if len(exp_pdb) == 2:
            in_all_vs_all = True
            exp_index, curr_pdb = exp_pdb
            exp_index = int(exp_index)
            RealtimeLogger.info(f"Exp index saves {exp_index}; len pdbs {len(all_pdbs)} {all_pdbs[exp_index]} {curr_pdb}")

        else:
            raise RuntimeError("Must be path to pdb or (experent number, path)")
    else:
        try:
            exp_index = all_pdbs.index(exp_pdb)
        except IndexError:
            exp_index = None

    RealtimeLogger.info(f"Exp index is {exp_index}; len all_pdbs {len(all_pdbs)}")
    if exp_index is not None and all_pdbs[exp_index] == curr_pdb:
        #Only perform comparisons for indexes exp_index+1 or higher
        pdbs = all_pdbs[exp_index+1:]
    else:
        assert 0
        pdbs = all_pdbs

    RealtimeLogger.info(f"Testong PDBS {len(pdbs)}")

    results_file = os.path.join(work_dir, f"{os.path.basename(curr_pdb)}_one_vs_all_{superposer}.dist")

    try:
        superposer = load_class(superposer)(job=job, work_dir=work_dir)
        return superposer.one_vs_all(curr_pdb, pdbs, table_out_file=results_file, in_all_vs_all=True, **kwds)
    except ImportError:

        map_job(job, one_vs_one, pdbs_to_run, pdbs, None, superposer, work_dir=work_dir,
            cores=cores, memory=mem, **kwds)

        with open(results_file, "w") as f:
            for i, fixed_pdb in enumerate(pdbs):
                output = one_vs_one(job, curr_pdb, fixed_pdb, superposer, *args, **kwds)
                print(",".join(map(str, results)), file=f)

def one_vs_one(job, exp_pdb, target_pdb, superposer, *args, memory="72G", work_dir=None, **kwds):
    if job is None:
        from toil.job import Job
        job = Job()

    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    if isinstance(exp_pdb, (list, tuple)) and len(exp_pdb) == 2:
        target_pdb = exp_pdb[1]
        exp_pdb = exp_pdb[0]

    try:
        superposer = load_class(superposer)(job=job, work_dir=work_dir)
        return superposer.one_vs_one(exp_pdb, target_pdb, **kwds)
    except ImportError:
        raise NotImplementedError("Superposer must have a one_vs_one method")

def combine_all_vs_all(job, superposer, pdbs, table_out_file=None, work_dir=None):
    if job is None:
        from toil.job import Job
        job = Job()

    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    distances = None
    results_file = lambda pdb: os.path.join(work_dir, f"{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")
    for pdb in pdbs[:-1]:
        distance_file = results_file(pdb)
        df = pd.read_csv(distance_file, index_col=False)
        if distances is None:
            distances = df
        else:
            distances = pd.concat((distances,df), axis=0)

    if table_out_file is not None:
        distances.to_csv(table_out_file)

    return distances

SUPERPOSERS = {}
class Superpose(Container):
    RULES = {"make_file_list": "make_file_list"}

    @classmethod
    def __init_subclass__(cls, *args, **kwds):
        global SUPERPOSERS
        SUPERPOSERS[cls.__name__.rsplit(".", 1)[-1]] = cls

    def make_file_list(self, key, file_list, format_in_path=True, basename=False):
        updated_file_list = self.tempfile()

        if isinstance(file_list, str) and os.path.isfile(file_list):
            with open(file_list) in fh:
                file_list = [pdb.rstrip() for pdb in fh]

        if not isinstance(file_list, (list, tuple)):
            raise RuntimeError("file_list must be list, tuple, or path to file")

        self.n_structures = len(file_list)
        self.file_list = updated_file_list

        with open(updated_file_list, "w") as new:
            for pdb in file_list:
                if basename:
                    pdb = self.format_in_path(None, pdb.rstrip())
                    print(os.path.basename(pdb.rstrip()), file=new)
                elif format_in_path:
                    pdb = self.format_in_path(None, pdb.rstrip(),
                        move_files_to_work_dir=True)
                    print(pdb, file=new)

        if format_in_path:
            return self.format_in_path(None, updated_file_list)

        return updated_file_list

    def _distribute_all_vs_all(self, file_list, table_out_file=None, n_jobs=None, **kwds):
        from toil.common import Toil
        from toil.job import Job
        options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
        options.logLevel = "DEBUG"
        options.clean = "always"

        options.logFile = "log.txt"
        #options.targetTime = 1
        options.retryCount = 0
        if options.provisioner == None:
            options.realTimeLogging = True
            options.maxLocalJobs = str(n_jobs) if n_jobs is not None else "8" #"96" #str(cores) if cores > 1 else str(multiprocessing.cpu_count())
            options.maxNodes = str(n_jobs) if n_jobs is not None else "8" #"96" #str(cores) if cores > 1 else str(multiprocessing.cpu_count())

        with Toil(options) as workflow:
            job = Job.wrapJobFn(start_one_vs_all_jobs, fullname(self), file_list, table_out_file=table_out_file, work_dir=self.work_dir, **kwds)
            distance_file = workflow.start(job)

        return distance_file

    def all_vs_all(self, pdb_list, table_out_file=None, distributed=False, **kwds):
        raise NotImplementedError

    def one_vs_all(self, experiment, pdb_list, table_out_file=None, **kwds):
        raise NotImplementedError

    def one_vs_one(self, moving_pdb_file, fixed_pdb_file, table_out_file=None, **kwds):
        raise NotImplementedError

    def cluster(self, pdb_list, linkage="ward", centroids=True, **kwds):
        distances = self.all_vs_all(pdb_list, **kwds)
        clusters = self.cluster_from_distances(pdb_list, distances, linkage=linkage, centroids=centroids)
        if centroids:
            return distances, clusters[0], clusters[1]
        return distances

    def cluster_from_distances(self, pdb_list, distances, linkage="ward", centroids=True):
        from sklearn.cluster import AgglomerativeClustering
        from sklearn.neighbors import NearestCentroid
        cluster_assignments = AgglomerativeClustering(linkage=linkage).fit_predict(distances)
        df = pd.DataFrame({"structure":pdb_list, "cluster":cluster_assignments})
        if centroids:
            clf = NearestCentroid()
            clf.fit(distances, cluster_assignments)
            centroids = pd.DataFrame({"centroid":clf.centroids_})
            df = pd.merge(df, centroids, left_on="cluster", right_index=True)
            centroid_dist = df.apply(lambda x: distances[x.index, x.cluster], axis=1)
            df.assign(centroid_dist=centroid_dist)
            centroids = df.groupby("centroid", as_index=False)["centroid_dist"].agg(['mean', 'max', np.std])
            return df, centroids
        return df



def fullname(o):
    klass = o.__class__
    module = klass.__module__
    if module == 'builtins':
        return klass.__qualname__ # avoid outputs like 'builtins.str'
    return module + '.' + klass.__qualname__

def load_class(m):
    module, klass = m.rsplit(".",1)
    mod = __import__(module, fromlist=[klass])
    return getattr(mod, klass)
