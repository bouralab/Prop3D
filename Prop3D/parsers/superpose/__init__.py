from numbers import Real
import os
import re
import shutil
from pathlib import Path

import pandas as pd

from Prop3D.parsers.container import Container
from Prop3D.util.toil import map_job, map_job_rv, map_job_rv_list

from toil.realtimeLogger import RealtimeLogger

def start_many(job, superposer, file_lists, table_out_file=None, work_dir=None, **kwds):
    RealtimeLogger.info("RUNNING MANY")
    if not isinstance(table_out_file, (list, tuple)) or not len(file_lists) == len(table_out_file):
        raise RuntimeError("table_out_files must be a list containing the same number of file_lists")
    rvs = []
    for file_list, out_file in zip(file_lists, table_out_file):
        ova = job.addChildJobFn(start_one_vs_all_jobs, superposer, file_list, out_file,
            work_dir=work_dir, **kwds)
        rvs.append(ova.rv())

    return rvs

def start_one_vs_all_jobs(job, superposer, file_list, table_out_file=None, work_dir=None, cores=1, mem="72G", memory="72G", cluster=None, **kwds):
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

        if table_out_file is None:
            table_out_file = os.path.basename(file_list_file)+".all_dist"
    else:
        pdbs = file_list

    if table_out_file is not None:
        results_file = lambda pdb: os.path.join(work_dir, f"{table_out_file}_{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")
    else:
        results_file = lambda pdb: os.path.join(work_dir, f"{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")

    #results_file = lambda pdb: os.path.join(work_dir, f"{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")
    pdbs_to_run = [(i, pdb) for i, pdb in enumerate(pdbs) if i<len(pdbs)-1] # and os.path.isfile(results_file(pdb))]

    n_by_n = kwds.pop("n_by_n", False)

    kwds.pop("table_out_file", None)

    similarity_files = map_job(job, one_vs_all, pdbs_to_run, pdbs[:], superposer, work_dir=work_dir, table_out_file=table_out_file,
        cores=cores, memory=mem, **kwds)

    if cluster is not None or (isinstance(cluster, bool) and cluster):
        similarity_threshold = 0.3 if isinstance(cluster, bool) else cluster
        return job.addFollowOnJobFn(combine_all_vs_all_and_cluster, superposer,
            pdbs, similarity_files, similarity_threshold=similarity_threshold, table_out_file=table_out_file,
            work_dir=work_dir).rv()

    return job.addFollowOnJobFn(combine_all_vs_all, superposer, pdbs, similarity_files, table_out_file=table_out_file,
        n_by_n=n_by_n, work_dir=work_dir).rv()

def one_vs_all(job, exp_pdb, file_list, superposer, *args, table_out_file=None, work_dir=None, **kwds):
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

    if table_out_file is not None:
        results_file = f"{table_out_file}_{os.path.basename(curr_pdb)}_one_vs_all_{superposer}.dist"
    else:
        results_file = f"{os.path.basename(curr_pdb)}_one_vs_all_{superposer}.dist"

    try:
        superposer = load_class(superposer)(job=job, work_dir=work_dir, is_distributed=True)
        out_file  =  superposer.one_vs_all(curr_pdb, pdbs, table_out_file=results_file, in_all_vs_all=True, **kwds)
        return out_file
    except ImportError:
        raise
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
        superposer = load_class(superposer)(job=job, work_dir=work_dir, is_distributed=True)
        return superposer.one_vs_one(exp_pdb, target_pdb, **kwds)
    except ImportError:
        raise NotImplementedError("Superposer must have a one_vs_one method")

def combine_all_vs_all(job, superposer, pdbs, similarity_files, table_out_file=None, n_by_n=False, work_dir=None):
    if job is None:
        from toil.job import Job
        job = Job()

    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    distances = None

    similarity_files = list(map_job_rv_list(similarity_files))

    if table_out_file is not None:
        results_file = lambda pdb: os.path.join(work_dir, f"{table_out_file}_{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")
    else:
        results_file = lambda pdb: os.path.join(work_dir, f"{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")

    #results_file = lambda pdb: os.path.join(work_dir, f"{os.path.basename(pdb)}_one_vs_all_{superposer}.dist")
    assert len(pdbs)>1, pdbs
    #print(pdbs)

    for similarity_file in similarity_files: #pdb in pdbs[:-1]:
        #similarity_file = results_file(pdb)
        similarities_file_name = similarity_file[0]
        similarities_file = job.fileStore.readGlobalFile(similarity_file[1])

        assert os.path.isfile(similarities_file), similarity_file
        df = pd.read_csv(similarities_file, sep="\t")
        assert len(df[df.isna()])==0, (df[df.isna()].shape, df.shape, similarities_file)
        #TMAlign specific
        df = df.rename(columns={"#PDBchain1":"chain1", "PDBchain2":"chain2", "TM1":"fixed_tm_score", "TM2":"moving_tm_score",
                "RMSD":"rmsd"})
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        shutil.copy(similarities_file, "/media/smb-rivanna/ed4bu/UrfoldServer/urfold_runs/deepurfold_figures/lrp_cluster")
        assert len(df[df.isna()])==0, (df[df.isna()].shape, df.shape, similarities_file)
        df = df[~df.isna()]
        if distances is None:
            distances = df
        else:
            distances = pd.concat((distances,df), axis=0)

        Path(similarities_file).unlink()

    if table_out_file is not None:
        distances.to_csv(table_out_file)


    if n_by_n:
        distances = load_class(superposer).make_n_by_n(distances)
        if table_out_file is not None:
            table_out_file += ".n_by_n"
            distances.to_csv(table_out_file)
        #     return table_out_file
        # return new_similarities

    if table_out_file is not None:
        return (table_out_file, job.fileStore.writeGlobalFile(table_out_file))
    else:
        return distances

def combine_all_vs_all_and_cluster(job, superposer, pdbs, similarities_file, similarity_threshold=0.2, table_out_file=None, work_dir=None):
    combined_similarities_file = job.addChildJobFn(combine_all_vs_all, superposer, pdbs,
        similarities_file, table_out_file=table_out_file, n_by_n=True, work_dir=work_dir).rv()
    return combined_similarities_file, job.addFollowOnJobFn(cluster, superposer,
        pdbs, combined_similarities_file, similarity_threshold=similarity_threshold,
        work_dir=work_dir).rv()

def cluster(job, superposer, pdbs, similarities_file, similarity_threshold=0.2, work_dir=None):
    try:
        superposer = load_class(superposer)(job=job, work_dir=work_dir, is_distributed=True)
        return superposer.cluster_from_similarities(pdbs, similarities_file, linkage="single",
            similarity_threshold=similarity_threshold, centroids=True)
    except ImportError:
        raise NotImplementedError("Superposer must have a cluster_from_similarities method")



SUPERPOSERS = {}
class Superpose(Container):
    RULES = {"make_file_list": "make_file_list"}
    DETACH = True #Needed so it doenst get converted to list, __init_subvlass doesn't work for subsub classes?

    def __init__(self, *args, is_distributed=False, **kwds):
        super().__init__(*args, **kwds)
        self.is_distributed = is_distributed

    @classmethod
    def __init_subclass__(cls, *args, **kwds):
        global SUPERPOSERS
        SUPERPOSERS[cls.__name__.rsplit(".", 1)[-1]] = cls
        super().__init_subclass__(*args, **kwds)

    def make_file_list(self, key, file_list, format_in_path=True, basename=False):
        updated_file_list = self.tempfile()

        if isinstance(file_list, str) and os.path.isfile(file_list):
            with open(file_list) in fh:
                file_list = [pdb.rstrip() for pdb in fh]

        if not isinstance(file_list, (list, tuple)):
            raise RuntimeError("file_list must be list, tuple, or path to file")

        self.n_structures = len(file_list)
        self.file_list = updated_file_list

        used_pdbs = set()

        with open(updated_file_list, "w") as new:
            for pdb in file_list:
                if basename:
                    pdb = self.format_in_path(None, pdb.rstrip())
                    print(os.path.basename(pdb.rstrip()), file=new)
                elif format_in_path:
                    #print("old pdb:", pdb, end="")
                    pdb_new = self.format_in_path(None, pdb.rstrip(),
                        move_files_to_work_dir=True, absolute_path=False)
                    #print("->", pdb_new)
                    print(pdb_new, file=new)
                    used_pdbs.add(pdb_new)

        if format_in_path:
            return self.format_in_path(None, updated_file_list)

        return updated_file_list

    def _distribute_all_vs_all(self, file_list, table_out_file=None, n_jobs=None, many_independent_jobs=False, **kwds):
        from toil.common import Toil
        from toil.job import Job

        if many_independent_jobs:
            assert isinstance(file_list, (list, tuple)) and isinstance(file_list[0], (list, tuple))
            assert isinstance(table_out_file, (list, tuple)) and len(table_out_file)==len(file_list)
            workflow_job = start_many
            name = f"many_independent_jobs_{os.path.basename(self.tempfile())}"
        else:
            workflow_job = start_one_vs_all_jobs
            name = table_out_file if table_out_file is not None else self.tempfile()


        options = Job.Runner.getDefaultOptions(f"./toilWorkflowRun_{name}")
        options.logLevel = "DEBUG"
        #options.clean = "always"

        options.logFile = "log.txt"
        #options.targetTime = 1
        options.retryCount = 0
        if options.provisioner == None:
            options.realTimeLogging = True
            options.maxLocalJobs = str(n_jobs) if n_jobs is not None else "8" #"96" #str(cores) if cores > 1 else str(multiprocessing.cpu_count())
            options.maxNodes = str(n_jobs) if n_jobs is not None else "8" #"96" #str(cores) if cores > 1 else str(multiprocessing.cpu_count())

        distance_files = []

        with Toil(options) as workflow:
            job = Job.wrapJobFn(workflow_job, fullname(self), file_list, table_out_file=table_out_file, work_dir=self.work_dir, **kwds)
            distance_file = workflow.start(job)

            if isinstance(distance_file, (list, tuple)):
                #Should all be at least tuples (name, fileID)
                if isinstance(distance_file[0], (list, tuple)):
                    if isinstance(distance_file[0][0], (list, tuple)):
                        for f in distance_file:
                            for name, fileID in f:
                                workflow.export_file(fileID, f"file://{Path.cwd().absolute()}/{name}")
                                distance_files.append(name)
                    else:
                        for name, fileID in distance_file:
                            workflow.export_file(fileID, f"file://{Path.cwd().absolute()}/{name}")
                            distance_files.append(name)

                elif len(distance_file) == 2:
                    name, fileID = distance_file
                    workflow.export_file(fileID, f"file://{Path.cwd().absolute()}/{name}")
                    distance_files.append(name)
                else:
                    raise RuntimeError(f"Can't read distrance file {distance_file}")

        if len(distance_files) == 1:
            return distance_file[0]
        return distance_files

    def all_vs_all(self, pdb_list, table_out_file=None, distributed=False, **kwds):
        raise NotImplementedError

    def one_vs_all(self, experiment, pdb_list, table_out_file=None, **kwds):
        raise NotImplementedError

    def one_vs_one(self, moving_pdb_file, fixed_pdb_file, table_out_file=None, **kwds):
        raise NotImplementedError

    def cluster(self, pdb_list, linkage="single", similarity_threshold=0.3, table_out_file=None, **kwds):
        is_distributed = int(kwds.get("distributed", False)) != 0
        if is_distributed or kwds.get("force", False) or (table_out_file is not None or not Path(table_out_file).is_file()):
            similarities = self.all_vs_all(pdb_list, table_out_file=table_out_file, n_by_n=True,
                cluster=similarity_threshold if is_distributed else None, **kwds)
        else:
            similarities = kwds["table_out_file"]+".n_by_n"
            assert Path(similarities).is_file()

        if not is_distributed:
            clusters = self.cluster_from_similarities(pdb_list, similarities, linkage=linkage, similarity_threshold=similarity_threshold, centroids=centroids)
            # if centroids:
            #     return similarities, clusters[0], clusters[1]
            return similarities, clusters
        else:
            return similarities

    def cluster_from_similarities(self, pdb_list, similarities_file, linkage="single", similarity_threshold=0.2, centroids=True):
        from sklearn.cluster import AgglomerativeClustering
        from sklearn.neighbors import NearestCentroid
        from sklearn.metrics import pairwise_distances_argmin_min
        print("Reading similarities", similarities_file)

        if self.is_distributed:
            similarities_file_name = similarities_file[0]
            similarities_file = self.job.fileStore.readGlobalFile(similarities_file[1])

        similarities = pd.read_csv(similarities_file)
        similarities = similarities.set_index("chain1").fillna(1.0)

        distances = 1-similarities.values
        distance_threshold = 1-similarity_threshold
        cluster_assignments = AgglomerativeClustering(linkage=linkage, affinity='precomputed').fit_predict(distances)
        df = pd.DataFrame({"structure":similarities.index, "cluster":cluster_assignments})

        if centroids:
            df = df.assign(centroid=np.nan)
            for cluster in range(cluster_assignments.n_clusters_):
                domain_idx = np.where(cluster_assignments.labels_==cluster)[0]
                cluster_distances = distances[domain_idx,domain_idx]
                centroid = np.linalg.norm(cluster_distances-cluster_distances.mean(), axis=1).argmin(axis=1)
                df.loc[df[df.cluster==cluster].index] = centroid

            # closest, _ = pairwise_distances_argmin_min(km.cluster_centers_, X)
            # clf = NearestCentroid()
            # clf.fit(distances, cluster_assignments)
            # import pdb; pdb.set_trace()
            # centroids = pd.DataFrame({"centroid":clf.centroids_}, index=df.index)
            # df = pd.merge(df, centroids, left_on="cluster", right_index=True)
            # centroid_dist = df.apply(lambda x: distances[x.index, x.cluster], axis=1)
            # df.assign(centroid_dist=centroid_dist)
            # centroids = df.groupby("centroid", as_index=False)["centroid_dist"].agg(['mean', 'max', np.std])
            # return df, centroids

        cluster_file = f"{similarities_file_name}.clusters"
        df.to_csv(cluster_file)

        if self.is_distributed:
            return (cluster_file, self.job.fileStore.writeGlobalFile(cluster_file))
        return cluster_file

        return df

    @staticmethod
    def make_n_by_n(similarities):
        return similarities



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
