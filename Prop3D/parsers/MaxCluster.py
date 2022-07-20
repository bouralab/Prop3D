import os
import re
import shutil
import tempfile
from math import ceil
import multiprocessing
from io import StringIO

from joblib import Parallel, delayed
import pandas as pd

from Prop3D.parsers.container import Container
from Prop3D.parsers.superpose import TMAlign
from Prop3D.util import safe_remove
from Prop3D.util.toil import partitions, map_job, map_job_rv, map_job_rv_list, loop_job_rv
from Prop3D.util.iostore import IOStore
from Prop3D.util.hdf import get_file

from toil.fileStores import FileID
from toil.realtimeLogger import RealtimeLogger

def start_one_vs_all_jobs(job, file_list, C=1, P=10, intermediate_file_store=None, further_parallelize=True, cores=1, mem="72G", memory="72G", **kwds):
    work_dir = kwds["work_dir"] = kwds.get("work_dir", job.fileStore.getLocalTempDir())
    if not isinstance(file_list, (list, tuple)):
        file_list_file = get_file(job, "file_list", file_list, work_dir=work_dir, cache=True)
        with open(file_list_file) as fh:
            pdbs = [pdb.rstrip() for pdb in fh]
    else:
        pdbs = file_list

    pdbs_to_run = [(i, pdb) for i, pdb in enumerate(pdbs) if not os.path.isfile(os.path.join(work_dir,
        f"{os.path.splitext(os.path.basename(pdb))[0]}.dist"))]
    #pdbs = list(enumerate(pdbs))



    assert cores == 1, "Cores why"

    if not further_parallelize or cores > 1:
        #Run sequentially or parallelize via joblib if multiple cores
        RealtimeLogger.info("Starting distributed start_one_vs_all_jobs 222 loop")
        distance_promises = list(loop_job(job, one_vs_all, pdbs_to_run, pdbs,
            intermediate_file_store=intermediate_file_store,
            cores=cores, memory=memory, **kwds))
    else:
        #Parallelize with toil
        RealtimeLogger.info("Starting distributed start_one_vs_all_jobs 222 parallel")
        distance_promises = map_job(job, one_vs_all, pdbs_to_run, pdbs,
            intermediate_file_store=intermediate_file_store,
            cores=cores, memory=mem, **kwds) #list(enumerate(file_list))
    return job.addFollowOnJobFn(cluster_from_distances, file_list,
        C=C, P=P, intermediate_file_store=intermediate_file_store, cores=cores,
        mem=mem, **kwds).rv()

def one_vs_all(job, exp_pdb, file_list, *args, distance="auto", work_dir=None, intermediate_file_store=None, cores=1, memory="72G", **kwds):

    if not isinstance(file_list, (list, tuple)):
        file_list_file = get_file(job, "file_list", file_list, work_dir=work_dir, cache=True)
        with open(file_list_file) as fh:
            all_pdbs = [pdb.rstrip() for pdb in fh]
    else:
        all_pdbs = file_list

    RealtimeLogger.info("START multprocess")

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
        assert 0, "WHY THH"
        try:
            exp_index = all_pdbs.index(exp_pdb)
        except IndexError:
            exp_index = None

    RealtimeLogger.info(f"Exp index is {exp_index}; len all_pdbs {len(all_pdbs)}")
    if exp_index is not None and all_pdbs[exp_index] == curr_pdb:
        #Only perform comparisons for indexes exp_index+1 or higher
        pdbs = all_pdbs[exp_index+1:]
    else:
        assert 0, "Why here"
        pdbs = all_pdbs

    RealtimeLogger.info(f"Testong PDBS {len(pdbs)}")

    if cores == 1:
        return _one_vs_all(job, exp_pdb, pdbs, *args, distance=distance, work_dir=work_dir,
            intermediate_file_store=intermediate_file_store, memory=memory, **kwds)

    # num_partitions = 100 #cores-1
    # partition_size = int(ceil(len(pdbs)/num_partitions))
    # partition_size = 1000
    # num_partitions = int(ceil(len(pdbs)/partition_size))
    # pdb_partitions = list(partitions(pdbs, partition_size))
    #
    # RealtimeLogger.info("Npartitions {}/{}; {}".format(len(pdbs), len(all_pdbs), len(pdb_partitions)))
    #
    # return Parallel(n_jobs=cores-2)(
    #     delayed(_one_vs_all)(None, exp_pdb, part, *args, distance=distance,
    #     full_file_list=all_pdbs, work_dir=work_dir,
    #     intermediate_file_store=intermediate_file_store, num_partition=i,
    #     total_partitions=num_partitions, **kwds) \
    #         for i, part in enumerate(pdb_partitions))

def _one_vs_all(job, exp_pdb, file_list, *args, distance="auto", num_partition=None, total_partitions=None, work_dir=None, intermediate_file_store=None, memory="72G", **kwds):
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    if job is None:
        from toil.job import Job
        job = Job()

    kwds.pop("cores", None)

    RealtimeLogger.info("Start partition {} with {} pdbs".format(num_partition, len(file_list)))

    tmalign = TMAlign(job=job, work_dir=work_dir, intermediate_file_store=intermediate_file_store, cleanup_when_done=False, **kwds.get("container_kwds", {}))
    results = []
    for i, fixed_pdb in enumerate(file_list):
        RealtimeLogger.info(f"Running {i}/{len(file_list)} in {exp_pdb}")
        try:
            output = tmalign(moving_pdb_file=exp_pdb[1], fixed_pdb_file=fixed_pdb)
        except RuntimeError as e:
            if "'return_code': -8" in str(e):
                try:
                    output = tmalign(moving_pdb_file=exp_pdb[1], fixed_pdb_file=fixed_pdb)
                except RuntimeError:
                    continue
        rmsd, moving_tm_score, fixed_tm_score = tmalign.get_rmsd_tm(include_fixed=True, stdout=output)
        results.append([
            os.path.basename(exp_pdb[1])[:7],
            os.path.basename(fixed_pdb)[:7],
            rmsd,
            moving_tm_score,
            fixed_tm_score])

    results = pd.DataFrame(results, columns=["PDB1", "PDB2", "rmsd", "TM-Score1", "TM-Score2"])

    results_file = os.path.join(work_dir, "{}.dist".format(
        os.path.splitext(os.path.basename(exp_pdb[1]))[0]))
    RealtimeLogger.info(f"Saving to {results_file}")
    results.to_csv(results_file)

    return

    #mc = MaxCluster(job=job, work_dir=work_dir, intermediate_file_store=intermediate_file_store, **kwds.get("container_kwds", {}))


    # try:
    #     return mc.one_vs_all(exp_pdb, file_list, *args, distance=distance, num_partition=num_partition, total_partitions=total_partitions, cores=1, **kwds)
    # except (SystemExit, KeyboardInterrupt):
    #     raise
    # except:
    #     raise
    #     store = IOStore.get(intermediate_file_store)
    #     import traceback as tb
    #     message = tb.format_exc()
    #     errfile = mc.tempfile()
    #     with open(errfile, "w") as f:
    #         print(message, file=f)
    #         print("", file=f)
    #         print(exp_pdb, file=f)
    #         print("", file=f)
    #         print("vs", file=f)
    #         print(len(file_list), file=f)
    #         for pdb in file_list:
    #             print(pdb, file=f)
    #
    #     exp_pdb =  exp_pdb[1] if isinstance(exp_pdb, (list, tuple)) else exp_pdb
    #     key = "maxcluster_intermediates/dist/FAILED_{}_part{}of{}.dist".format(
    #         os.path.splitext(os.path.basename(exp_pdb))[0],
    #         num_partition, total_partitions)
    #     store.write_output_file(errfile, key)
    #
    #     safe_remove(errfile)
    #
    #     RealtimeLogger.info("MaxCluster FAILED {}".format(key))
    #
    #     del store, key, errfile

def cluster_from_distances(job, pdbs, C=1, P=10, work_dir=None, intermediate_file_store=None, cores=20, memory="72G", **kwds):
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    distances = None
    for pdb in pdbs[:-1]:
        distance_file = os.path.join(work_dir, f"{os.path.splitext(os.path.basename(pdb))[0]}.dist")
        df = pd.read_csv(distance_file)
        if distances is None:
            distances = df
        else:
            distances = pd.concat((distances,df), axis=0)

    all_dist = os.path.join(work_dir, "all_distances.csv")
    distances.to_csv(all_dist)

    if C==0:
        return all_dist

    mc = MaxCluster(job=job, work_dir=work_dir, intermediate_file_store=intermediate_file_store)
    return mc.cluster_from_distances(file_list, distances, C=C, P=P, **kwds)

class MaxCluster(Container):
    IMAGE = 'docker://edraizen/maxcluster:latest'
    RULES = {"make_file_list": "make_file_list"}
    ARG_START = "-"
    ENTRYPOINT = "/opt/maxcluster/maxcluster"
    PARAMETERS = [
        #Structure Comparison
        (":l", "make_file_list"),
        (":e", "path:in"),
        (":p", "path:in"),

        (":L", "str"),
        (":d", "str"),
        (":N", "str"),
        (":rmsd", "store_true"),
        (":gdt", "str"),
        (":av", "store_true"),
        (":TSF", "str"),
        (":rankMS", "store_true"),
        (":noscore", "store_true"),
        (":bb", "store_true"),
        (":urmsd", "store_true"),
        (":maxRMSD", "store_true"),
        (":TM", "store_true"),
        (":seed", "str"),
        (":all", "store_true"),
        (":super", "str"),
        (":nourmsd", "store_true"),
        (":noalign", "store_true"),
        (":lws", "str"),
        (":sequence_independent", "store_true", ["-in"]),
        (":matrix", "store_true"),
        (":opt", "str"),
        (":open", "str"),
        (":aff", "str"),

        #Clustering
        (":C", "str"),
        (":T", "str"),
        (":Tm", "str"),
        (":a", "str"),
        (":is", "str"),
        (":ms", "str"),
        (":s", "str"),
        (":P", "str"),
        (":R", "path:out"),
        (":Rl", "path:out"),
        (":Ru", "path:out"),
        (":F", "path:in"),
        (":M", "path:in"),

        (":log", "path:out"),

        #Contact Maps
        (":contact", "store_true"),
        (":dist", "str"),
        (":gap", "str"),

        #MaxSubDom
        (":O", "str"),
        (":mS", "str"),
        (":mP", "str"),
        (":i", "str"),
    ]

    @classmethod
    def create_output_name(cls, kwds, prefix="MaxCluster", suffix=".log", work_dir=None):
        if work_dir is None and isinstance(cls, MaxCluster):
            work_dir = cls.work_dir
        elif work_dir is None:
            work_dir = ""

        fname = prefix
        for k, v in kwds.items():
            if k in ["e", "p", "F", "M"] or (k=="l" and isinstance(v, str) and os.path.isfile(v)):
                fname += "_{}={}".format(k, os.path.splitext(os.path.basename(v))[0])
            elif k == "l":
                #A list
                fname += "_l={}".format(len(v))
            elif k in ["log", "R", "Rl", "Ru"]:
                #Skip log
                pass
            else:
                fname += "_{}={}".format(k, v)
        fname += suffix
        return os.path.join(work_dir, fname)

    def __call__(self, *args, **kwds):
        log_file = self.create_output_name(kwds, suffix="")

        original_return_files = self.return_files

        for k, ext in [("log", ".log"), ("R", ".dist"), ("Rl", ".distlite"), ("Ru", ".udist")]:
            if k in kwds and isinstance(kwds[k], bool):
                if kwds[k]:
                    kwds[k] = log_file+ext
                else:
                    del kwds[k]

        if "R" in kwds or "Rl" in kwds or "Ru" in kwds:
            self.return_files = True
            if "log" not in kwds:
                kwds["log"] = log_file+".log"

        if "log" in kwds:
            self.return_files = True
            self.log_file = kwds["log"]
            try:
                output = super().__call__(*args, **kwds)
            except (RuntimeError, AssertionError) as e:
                if "chown: changing ownership" not in str(e):
                    with open(self.log_file) as f:
                        msg = f.read()
                    raise RuntimeError("MaxCluster Failed {}: {}".format(str(e), msg))


        else:
            #Container handles exception
            output = super().__call__(*args, **kwds)
            self.log_file = io.StringIO(output)

        self.return_files = original_return_files

        return output

    def make_file_list(self, key, file_list, format_in_path=True):
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
                RealtimeLogger.info(f"Adding {pdb}")
                #print(self.format_in_path(None, pdb.rstrip(),
                #    move_files_to_work_dir=False), file=new)

        if format_in_path:
            return self.format_in_path(None, updated_file_list)

        return updated_file_list

    def all_vs_all(self, file_list, C=1, P=10, distributed="auto", distance="auto", cores=1, memory="1G", **kwds):
        if isinstance(file_list, str):
            with open(file_list) as fh:
                n_structures = sum(1 for pdb in fh)
        else:
             n_structures = len(file_list)

        if distributed == "auto":
            distributed = n_structures > 10000

        RealtimeLogger.info("Start all vs all: distributed={}".format(distributed))

        if distributed:
            # RealtimeLogger.info("MAKING FILE LIST")
            # file_list_key = "maxcluster_intermediates/dist/file_list"
            #
            # if self.intermediate_file_store is not None:
            #     save_file_list = True
            #     if self.intermediate_file_store.exists(file_list_key):
            #         file_list = self.tempfile()
            #         self.intermediate_file_store.read_input_file(file_list_key, file_list)
            #         save_file_list = False
            #     else:
            #         file_list = self.make_file_list(None, file_list, format_in_path=False)
            #
            #     if save_file_list:
            #         self.intermediate_file_store.write_output_file(file_list, file_list_key)
            #
            # RealtimeLogger.info("DONE MAKING FILE LIST {}".format(file_list))

            if False: #cores > 1:
                #Parallelize with joblib
                kwds["C"] = C
                kwds["P"] = P
                kwds["l"] = file_list
                kwds["log"] = kwds.get("log", True)
                kwds["work_dir"] = self.work_dir
                kwds["distance"] = distance
                kwds["intermediate_file_store"] = self.intermediate_file_store
                kwds["container_kwds"] = dict(
                    return_files=self.return_files,
                    force_local=self.force_local,
                    fallback_local=self.fallback_local,
                    detach=self.detach)
                log_files = Parallel(n_jobs=cores, backend="multiprocessing")(
                    delayed(one_vs_all)(self.job, (i, f), file_list, **kwds) \
                        for i, f in enumerate(file_list))
                return log_files

            else:
                #Parallelize through toil
                kwds["C"] = C
                kwds["P"] = P
                kwds["distance"] = distance
                kwds["log"] = kwds.get("log", True)
                kwds["work_dir"] = self.work_dir
                kwds["intermediate_file_store"] = self.intermediate_file_store
                kwds["cores"] = cores

                if kwds.get("further_parallize", True):
                    kwds["mem"] = memory
                    kwds["memory"] = "1G"
                else:
                    kwds["memory"] = memory

                RealtimeLogger.info("Starting distributed toil start_one_vs_all_jobs")

                if not hasattr(self.job, "fileStore"):
                    from toil.common import Toil
                    from toil.job import Job
                    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
                    options.logLevel = "INFO"
                    options.clean = "always"
                    #options.targetTime = 1
                    options.retryCount = 0
                    options.logFile = "out.txt"
                    if options.provisioner == None:
                        options.realTimeLogging = True
                        options.maxLocalJobs = "96" #str(cores) if cores > 1 else str(multiprocessing.cpu_count())
                        options.maxNodes = "96" #str(cores) if cores > 1 else str(multiprocessing.cpu_count())
                    kwds["cores"] = 1

                    with Toil(options) as workflow:
                        job = Job.wrapJobFn(start_one_vs_all_jobs, file_list, **kwds)
                        distance_file = workflow.start(job)

                    return distance_file
                else:
                    assert 0, "Start here??"
                    file_list = self.job.fileStore.writeGlobalFile(file_list)
                    return self.job.addChildJobFn(start_one_vs_all_jobs, file_list, **kwds).rv()
        else:
            kwds["C"] = C
            kwds["P"] = P
            kwds["l"] = file_list
            kwds["log"] = kwds.get("log", True)
            if distance=="auto" or distance.lower() not in ["maxsub", "tm", "rmsd"]:
                #Default to maxsub, no options, unless:
                if n_structures > 2000:
                    kwds["rmsd"] = True
            elif distance.lower() == "tm":
                kwds["TM"] = True
            elif distance.lower() == "rmsd":
                kwds["rmsd"] = True
            else:
                #Maxsub was chosen no options
                pass

            output = self(**kwds)

            if isinstance(output, dict):
                if self.intermediate_file_store is not None or hasattr(self.job, "fileStore"):
                    for arg, path in output.items():
                        RealtimeLogger.info("OUTFILE: {}={}".format(arg, path))
                        try:
                            self.intermediate_file_store.write_output_file(path,
                                "maxcluster_intermediates/{}/{}".format(arg, os.path.basename(path)))
                        except:
                            output[arg] = self.job.fileStore.writeGlobalFile(path)
            elif isinstance(output, (list, tuple)):
                if self.intermediate_file_store is not None or hasattr(self.job, "fileStore"):
                    for i, path in enumerate(output[:]):
                        try:
                            self.intermediate_file_store.write_output_file(path,
                                "maxcluster_intermediates/{}/{}".format(arg, ps.path.basename(path)))
                        except:
                            output[i] = self.job.fileStore.writeGlobalFile(path)


            return output

    def one_vs_all(self, exp_pdb, file_list, distance="auto", save_distance_file=False, num_partition=None, total_partitions=None, full_file_list=None, cores=1, **kwds):
        # if cores > 1:
        #     RealtimeLogger.info("Will start multprocess")
        #     #Run with multiple processes
        #     return one_vs_all(self.job, exp_pdb, file_list, distance=distance, work_dir=self.work_dir,
        #         intermediate_file_store=self.intermediate_file_store.store_string if self.intermediate_file_store is not None else None,
        #         cores=cores, **kwds)

        in_all_vs_all = False
        if isinstance(exp_pdb, (list, tuple)):
            if len(exp_pdb) == 2:
                in_all_vs_all = True
                exp_index, exp_pdb = exp_pdb
            else:
                raise RuntimeError("Must be path to pdb or (experent number, path)")
        else:
            exp_index = None

        if not os.path.isfile(exp_pdb):
            return

        # if not isinstance(file_list, (list, tuple)):
        #     fname = "file_list"+"_{}".format(num_partition) if num_partition is not None else ""
        #     file_list = get_file(self.job, fname, file_list, work_dir=self.work_dir, cache=True)
        #     with open(file_list) as fh:
        #         file_list = [pdb.rstrip() for pdb in fh]
        #
        # if full_file_list is not None and not isinstance(full_file_list, (list, tuple)):
        #     full_fname = "full_file_list"+"_{}".format(num_partition) if num_partition is not None else ""
        #     file_list = get_file(self.job, full_fname, file_list, work_dir=self.work_dir, cache=True)
        #     with open(full_file_list) as fh:
        #         full_file_list = [pdb.rstrip() for pdb in fh]

        if exp_index is None:
            try:
                exp_index = file_list.index(exp_pdb)
            except (IndexError, ValueError):
                pass

        kwds["e"] = exp_pdb
        if num_partition is None and in_all_vs_all and file_list[exp_index] == exp_pdb:
            #Only perform comparisons for indexes exp_index+1 or higher
            RealtimeLogger.info("MaxCLuster {}: {} {}".format(exp_index, len(file_list[exp_index+1:]), len(file_list)))
            if len(file_list[exp_index+1:]) == 0:
                #It is the last one, no need to compare against itself
                return
            kwds["l"] = file_list[exp_index+1:]
            file_list[exp_index+1:]
        else:
            kwds["l"] = file_list

        if in_all_vs_all or save_distance_file:
            if num_partition is not None:
                #Make sure log files have different names
                kwds["log"] = self.create_output_name(kwds, suffix=".log{}".format(num_partition))
                RealtimeLogger.info(f"LOGFILE {kwds['log']}")
                if os.path.isfile(kwds["log"]):
                    return
            else:
                kwds["log"] = self.create_output_name(kwds, suffix=".log") #True
                RealtimeLogger.info(f"LOGFILE {kwds['log']}")
                if os.path.isfile(kwds["log"]):
                    return

            if num_partition is not None:
                dist_file_name = os.path.join(self.work_dir, "{}_part{}of{}.dist".format(
                    os.path.splitext(os.path.basename(exp_pdb))[0],
                    num_partition, total_partitions))
            else:
                dist_file_name = os.path.join(self.work_dir, "{}.dist".format(
                    os.path.splitext(os.path.basename(exp_pdb))[0]))

            if os.path.isfile(dist_file_name):
                return dist_file_name

        if "R" in kwds:
            del kwds["R"]

        # if distance=="auto" or distance.lower() not in ["maxsub", "tm", "rmsd"]:
        #     #Default to maxsub, no options, unless:
        #     if len(file_list) > 2000:
        #         kwds["rmsd"] = True
        # elif distance.lower() == "tm":
        #     kwds["TM"] = True
        # elif distance.lower() == "rmsd":
        #     kwds["rmsd"] = True
        # else:
        #     #Maxsub was chosen no options
        #     pass

        if distance.lower() == "rmsd":
            kwds["rmsd"] = True
        use_rmsd = kwds.get("rmsd", False)


        output = self(**kwds)

        if in_all_vs_all:
            try:
                if hasattr(self.log_file, 'read'):
                    log_file_name = self.create_output_name(kwds, suffix="_{}.log".format(num_parition))
                    with open(log_file_name, "w") as fh:
                        shutil.copyfileobj(self.log_file, fh)
                    self.log_file.seek(0)
                else:
                    log_file_name = self.log_file

                if full_file_list is not None:
                    distances = self.get_distances(log_file_name, full_file_list,
                        use_rmsd=use_rmsd)
                else:
                    distances = self.get_distances(log_file_name, file_list,
                        use_rmsd=use_rmsd)



                distances.to_csv(dist_file_name, index=False)
                RealtimeLogger.info(f"distances are {distances}")
                RealtimeLogger.info(f"Saved to {dist_file_name}")
                if self.intermediate_file_store is not None:
                    dist_key = f"maxcluster_intermediates/dist/{os.path.basename(dist_file_name)}"
                    self.intermediate_file_store.write_output_file(dist_file_name, dist_key)

                #safe_remove(dist_key)
                #safe_remove(log_file_name)

                return dist_file_name

            except (SystemExit, KeyboardInterrupt):
                if not hasattr(self.log_file, 'read'):
                    save_remove(self.log_file)
                raise
            except:
                raise
                import traceback as tb
                msg = tb.format_exc()
                RealtimeLogger.info("FAILED MAXCLUSTER All VS ALL: {}".format(msg))
                return output

        return output

    def one_vs_one(self, exp_pdb, ref_pdb, **kwds):
        pass

    def cluster_from_distances(self, file_list, distance_file, C=1, P=10, distributed=False, **kwds):
        # if isinstance(distance_file, (list, tuple)):
        #     #Flatten promised return values
        #     from toil.fileStores import FileID
        #     distance_fileIDs = list(map_job_rv_list(distance_file))
        #     if len(distance_fileIDs) == 0:
        #         #Using intermediate_file_store
        #         use_intermediate_file_store = True
        #         distance_fileIDs = self.intermediate_file_store.list_input_directory("maxcluster_intermediates/dist")
        #
        #     distance_file = os.path.join(self.work_dir, "MaxCluster-all_vs_all.dist")
        #     with open(distance_file, "w") as fh:
        #         if isinstance(distance_file[0], FileID):
        #             for fileID in distance_fileIDs:
        #                 shutil.copyfileobj(
        #                     self.job.fileStore.readGlobalFileStream(fileID), fh)
        #                 self.job.fileStore.deleteGlobalFile(fileID)
        #         elif use_intermediate_file_store:
        #             download_file_store = not hasattr(self.intermediate_file_store, "path_prefix")
        #
        #             for distance_fname in distance_fileIDs:
        #                 if download_file_store:
        #                     distance_fname_local = os.path.join(self.work_dir, os.path.basename(distance_fname))
        #                     self.intermediate_file_store.read_input_file(distance_fname, distance_fname_local)
        #                 else:
        #                     distance_fname_local = distance_fname_local
        #                 with open(distance_fname_local) as dist_file:
        #                     for distance_fname in distance_fileIDs:
        #                         with open(distance_fname) as dist_file:
        #                             shutil.copyfileobj(dist_file, fh)
        #                 if download_file_store:
        #                     safe_remove(distance_fname_local)
        #                     self.intermediate_file_store.remove_file(distance_fname)
        #         else:
        #             with open(distance_file) as fh:
        #                 for distance_fname in distance_fileIDs:
        #                     with open(distance_fname) as dist_file:
        #                         shutil.copyfileobj(dist_file, fh)

        kwds["C"] = C
        kwds["P"] = P
        kwds["l"] = file_list
        if kwds.get("rmsd", False):
            kwds["F"] = distance_file
        else:
            kwds["M"] = distance_file

        return self(**kwds)

    @staticmethod
    def get_distances(log_file, file_list, use_rmsd=False):
        if use_rmsd:
            pdbs_file_and_dist_re = re.compile("^INFO  : \d+\. (?P<PDB1>.+) vs. (?P<PDB2>.+)  RMSD=(?P<RMSD>[ \d\.]+) (Pairs=(?P<Pairs>.+), rRMSD=(?P<rRMSD>[\d\.]+) ((?P<rRMSD_Zscore>[ \d\.]+))), URMSD=(?P<URMSD>[ \d\.]+) (rURMSD=(?P<rURMSD>[ \d\.]+))")
            dist_cols = ["PDB1", "PDB2", "RMSD", "Pairs", "rRMSD", "rRMSD_Zscore", "URMSD", "rURMSD"]
        else:
            pdbs_file_and_dist_re = re.compile("^INFO  : \d+\. (?P<PDB1>.+) vs. (?P<PDB2>.+)  Pairs=(?P<Pairs>.+), RMSD=(?P<RMSD>[ \d\.]+), MaxSub=(?P<MaxSub>[ \d\.]+), TM=(?P<TM>[ \d\.]+), MSI=(?P<MSI>[ \d\.]+)")
            dist_cols = ["PDB1", "PDB2", "Pairs", "RMSD", "MaxSub", "TM", "MSI"]

        results = []
        with open(log_file) as log:
            for i, line in enumerate(log):
                m = pdbs_file_and_dist_re.match(line)
                if m:
                    groups = m.groupdict()
                    groups["PDB1"] = os.path.splitext(os.path.basename(groups["PDB1"].rstrip()))[0]
                    groups["PDB2"] = os.path.splitext(os.path.basename(groups["PDB2"].rstrip()))[0]
                    results.append([groups[k].strip() if isinstance(groups[k], str) \
                        else groups[k] for k in dist_cols])

        distances = pd.DataFrame(results, columns=dist_cols)

        return distances


    @classmethod
    def get_centroid(cls, log_file=None):
        if hasattr(cls.log_file):
            log_file = cls.log_file
        if log_file is None:
            raise RuntimeError("Invalid log file")

        parse_centroids = False
        best_centroid_size = None
        best_centroid_file = None

        with open(log_file) as log:
            for line in log:
                job.log("LINE: {}".format(line.rstrip()))
                if not parse_centroids and line.startswith("INFO  : Centroids"):
                    parse_centroids = True
                    next(log)
                    next(log)
                elif parse_centroids and line.startswith("INFO  : ="):
                    break
                elif parse_centroids:
                    job.log("PARSING LINE: {}".format(line.rstrip()))
                    fields = line.rstrip().split()
                    size, pdb_file = fields[-3], fields[-1]
                    if best_centroid_size is None or int(size)>best_centroid_size:
                        best_centroid_size = int(size)
                        best_centroid_file = pdb_file

        return best_centroid_file

    @classmethod
    def get_clusters(cls, log_file=None):
        if hasattr(cls, "log_file"):
            log_file = cls.log_file
        if log_file is None:
            raise RuntimeError("Invalid log file")

        clusters = []
        centroids = []
        parse_clusters = False
        parse_centroids = False
        with open(log_file) as log:
            for line in log:
                if not parse_clusters and "Clusters @ Threshold" in line:
                    parse_clusters = True
                    next(log)
                    next(log)
                elif parse_clusters and line.startswith("INFO  : ="):
                    parse_clusters = False
                elif not parse_centroids and "Centroids" in line:
                    parse_centroids = True
                    next(log)
                    next(log)
                elif parse_centroids and line.startswith("INFO  : ="):
                    parse_centroids = False
                    break
                elif parse_clusters:
                    fields = line.rstrip().split()
                    pdb_num, cluster, pdb_file = fields[-4], fields[-2], fields[-1]
                    clusters.append((int(pdb_num), pdb_file, int(cluster)))
                elif parse_centroids:
                    fields = line.rstrip().split()
                    cluster, pdb_num, spread, pdb_file = fields[-6], fields[-4], fields[-2], fields[-1]
                    centroids.append((int(cluster), spread, int(pdb_num), pdb_file))

        clusters = pd.DataFrame(clusters, columns=["item", "pdb_file", "structure_clusters"])
        centroids = pd.DataFrame(centroids, columns=["structure_clusters", "spread", "centroid_item", "centroid_pdb_file"])
        clusters = pd.merge(clusters, centroids, on="structure_clusters")

        return clusters

    @classmethod
    def get_distances_from_dist_file(cls, distance_file=None):
        if hasattr(cls, "return_files") and isinstance(cls.return_files, dict) and cls.return_files.get("R"):
            distance_file = cls.return_files["R"]
        if distance_file is None:
            raise RuntimeError("Invalid log file")

        distances = []
        parse_distances = False
        with open(distance_file) as dist:
            for line in dist:
                if not parse_distances and "# Distance records" in line:
                    parse_distances = True
                    next(dist)
                    next(dist)
                elif parse_distances:
                    #DIST :     79     82  1000.0000
                    fields = line.rstrip().split()
                    item1, item2, distance = fields[-3:]
                    distances.append((int(item1), int(item2), float(distance)))
                    distances.append((int(item2), int(item1), float(distance)))

            distances = pd.DataFrame(distances, columns=["item1", "item2", "distances"])
            distances = distances.set_index(["item1", "item2"])

            return distances


    @classmethod
    def get_hierarchical_tree(cls, log_file=None):
        #https://stackoverflow.com/questions/31033835/newick-tree-representation-to-scipy-cluster-hierarchy-linkage-matrix-format
        if hasattr(cls, "log_file"):
            log_file = cls.log_file
        if log_file is None:
            raise RuntimeError("Invalid log file")

        import networkx as nx
        import pandas as pd
        nx_tree = nx.DiGraph()
        tree = nx.Graph()
        Z = []
        parse_tree = False

        pdbs = {}
        nodes = {}

        class Node(object):
            def __init__(self, pdb1, pdb2, distance):
                self.pdb1 = pdb1
                self.pdb2 = pdb2
                self.distance = distance

            def __str__(self):
                return "({}, {})".format(self.pdb1, self.pdb2)


        with open(log_file) as log:
            for line in log:
                if not parse_tree and line.startswith("INFO  : Hierarchical Tree"):
                    parse_tree = True
                    next(log)
                    next(log)
                elif parse_tree and line.startswith("INFO  : ="):
                    break
                elif parse_tree:
                    #print(line)
                    _, node, info = line.rstrip().split(":")
                    node = int(node.strip())
                    nx_tree.add_node(node)

                    fields = [f.strip() for f in info.split()]
                    item1, item2 = list(map(int, fields[:2]))
                    distance = float(fields[2])

                    if item1>0 and item2>0:
                        pdb1, pdb2 = fields[3:5]
                        pdbs[item1] = pdb1[:7]
                        pdbs[item2] = pdb2[:7]
                        nodes[node] = Node(pdb1[:7], pdb2[:7], distance)
                        #nx_tree.add_node(item1, pdb=pdb1)
                        #nx_tree.add_node(item2, pdb=pdb2)
                    elif item1>0 and item2<=0:
                        #Node 2 already in graph
                        pdb1 = fields[3]
                        pdbs[item1] = pdb1[:7]
                        nodes[node] = Node(pdb1[:7], nodes[item2], distance)
                        #nx_tree.add_node(item1, pdb=pdb1)
                    elif item1<=0 and item2>0:
                        #Node 1 already in graph
                        pdb2 = fields[3]
                        pdbs[item2] = pdb2[:7]
                        nodes[node] = Node(nodes[item1], pdb2[:7], distance)
                        #nx_tree.add_node(item2, pdb=pdb2)
                    else:
                        nodes[node] = Node(nodes[item1], nodes[item2], distance)

                    #Z.append([])

                    Z.append([node, item1, item2, distance])

                    nx_tree.add_edge(node, item1)
                    nx_tree.add_edge(node, item2)
                    tree.add_edge(item1, item2)

            root_node = min(nodes.items(), key=lambda x: x[0])

            Z = pd.DataFrame(Z, columns=["Node", "Item1", "Item2", "Distance"])

            Z = Z.assign(count=0)
            for i in range(len(Z)):
                row = Z.iloc[i]
                count = Z[Z["Node"]==row["Item1"]].iloc[0]["count"] if row["Item1"]<=0 else 1
                count += Z[Z["Node"]==row["Item2"]].iloc[0]["count"] if row["Item2"]<=0 else 1
                Z.loc[Z["Node"]==row["Node"], "count"] = count

            #print(Z)

            Z2 = Z.copy()
            top = Z[["Item1", "Item2"]].max().max()+1
            #print("TOP", top)
            update_nodes = -1*Z[Z[["Item1", "Item2"]]<=0]+top
            new_item1 = update_nodes.dropna(subset=["Item1"])["Item1"]
            new_item2 = update_nodes.dropna(subset=["Item2"])["Item2"]
            Z.loc[new_item1.index, "Item1"] = new_item1
            Z.loc[new_item2.index, "Item2"] = new_item2
            Z = Z.drop(columns=["Node"])

            import numpy as np
            dist = np.zeros((top-1, top-1))
            _dist = Z2[(Z2["Item1"]>0)&(Z2["Item2"]>0)]
            dist[_dist["Item1"], _dist["Item2"]] = _dist["Distance"]

        return nx_tree, tree, Z, dist, root_node

    def get_distances_and_transformations(self, distance_file):
        pdbs = {}
        dist_mat1 = dist_mat2 = None
        with open(distance_file) as f:
            for line in f:
                if line.startswith("PDB"):
                    num, pdb = line.strip().split()[2:]
                    pdbs[num] = pdb
                elif line.startswith("# Transformtion"):
                    dist_mat1 = dist_mat2 = pd.empty((len(pdbs), len(pdbs)))
                elif line.startswith("TRANS") and dist_mat1 is not None:
                    fields = line.rstrip().split()
                    dist_mat1[int(fields[1]), int(fields[2])] = float(fields[3])
                    dist_mat1[int(fields[2]), int(fields[1])] = float(fields[3])
                    dist_mat2[int(fields[1]), int(fields[2])] = float(fields[4])
                    dist_mat2[int(fields[1]), int(fields[2])] = float(fields[4])

        from sklearn.cluster import KMeans
        kmeans1_2 = KMeans(init='k-means++', n_clusters=2, n_init=10)


def get_aligned_residues(file_list, transformation_file):
    with open(file_list) as fl:
        pdb_fnames = [l.strip() for l in fl]

    positions = {}
    parse_pos = False
    with open(transformation_file) as log:
        for line in log:
            job.log("LINE: {}".format(line.rstrip()))
            if not parse_pos and line.startswith("POS  : <ID> <ID>"):
                parse_pos = True
                next(log)
            elif parse_pos and line.startswith("#"):
                break
            elif parse_pos:
                _, info = line.rstrip().split(":")
                fields = [f.strip() for f in info.split()]
                pdbs = tuple([pdb_fnames[int(i)] for i in fields[:2]])
                positions[pdbs] = {[int(p.strip()) for p in pos.split()] for pos in fields[2:]}
    return positions

def get_msa(file_list, transformation_file, log_file):
    positions = get_aligned_residues(file_list, transformation_file)



def make_distance_file(distances, file_list, fh, header=True, entries=True):
    if header:
        from scipy.special import comb

        if isinstance(file_list, str) and os.path.isfile(file_list):
            with open(file_list) in fh:
                file_list = [pdb.rstrip() for pdb in fh]

        print("""###################################
# Maxcluster Distance Matrix File
###################################""", file=fh)
        print("SIZE :", int(comb(len(file_list), 2)), file=fh)
        print("""##########################
# Chain records
# PDB  : <ID> <Filename>
##########################""", file=fh)

    for i, pdb in enumerate(file_list):
        print("PDB : {0: >6}".format(i), pdb, file=fh)

    print("""###################################
# Distance records
# DIST : <ID1> <ID2> <Distance>
###################################""", file=fh)

    if entries:
        for row in distances.itertuples():
            print("DIST : {0: >6} {1: >6} {2}".format(row.PDB1,
                row.PDB2, row.RMSD), file=fh)

def parallel_cluster():
    kwds["M"] = make_distance_file()
    run_maxcluster(*args, **kwds)

def parallel(file_list, args, kwds, joblib=True):
    with open(file_list) as fh:
        pdbs = [pdb.rstrip() for pdb in fh]

    from joblib import Parallel, delayed

    if joblib:
        all_vs_all = [one_vs_one for _vs_all in Parallel(n_jobs=-1)(delayed(one_vs_all) \
            (exp_index, exp_pdb, file_list, args, kwds) for exp_index, exp_pdb in \
            enumerate(pdbs) for one_vs_one in _vs_all)]

# def run_maxcluster(*args, **kwds):
#     work_dir = kwds.pop("work_dir", None)
#     docker = kwds.pop("docker", True)
#     job = kwds.pop("job", None)
#
#     if work_dir is None:
#         work_dir = os.getcwd()
#
#     if "file_list" in kwds and not "l" in kwds:
#         kwds["l"] = kwds.pop("file_list")
#     else:
#         kwds.pop("file_list", None)
#
#     log = kwds.get("log", False)
#     if log and not isinstance(log, str):
#         f = tempfile.NamedTemporaryFile(dir=work_dir, suffix=".log", delete=False)
#         f.close()
#         kwds["log"] = f.name
#
#     file_kwds = ["log", "e", "p", "l", "R", "Rl", "Ru", "F", "M"]
#     in_file_kwds = ["e", "p", "l", "F", "M"]
#     parameters = ["-"+a for a in args]
#     for k, v in list(kwds.items()):
#         if k not in file_kwds:
#             parameters += ["-{}".format(k), str(v)]
#     job.log("ORIG PARAMS: {}".format(parameters))
#     file_parameters = {k:v for k, v in list(kwds.items()) if k in file_kwds}
#
#     if docker and apiDockerCall is not None and job is not None:
#         for k,f in list(file_parameters.items()):
#             if k in in_file_kwds and not os.path.abspath(os.path.dirname(f)) == os.path.abspath(work_dir):
#                 shutil.copy(f, work_dir)
#             job.log("BASENAMING: {}".format(os.path.basename(f)))
#             parameters += ["-{}".format(k), os.path.basename(f)]
#
#         oldcwd = os.getcwd()
#         os.chdir(work_dir)
#         try:
#             out = apiDockerCall(job,
#                           'edraizen/maxcluster:latest',
#                           working_dir="/data",
#                           volumes={work_dir:{"bind":"/data", "mode":"rw"}},
#                           parameters=parameters
#                           )
#         except (SystemExit, KeyboardInterrupt):
#             raise
#         except Exception as e:
#             job.log("Error: {} {}".format(type(e), e))
#             job.log("FILE LIST IS [{}]".format(open(file_parameters["l"]).read()))
#             return None
#             #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
#             #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
#             #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)
#         os.chdir(oldcwd)
#     else:
#         file_args = []
#         for k,f in list(file_parameters.items()):
#             parameters += ["-{}".format(k), f]
#         args = [maxcluster_path]+file_args+parameters
#         try:
#             out = subprocess.check_output(args)
#         except (SystemExit, KeyboardInterrupt):
#             raise
#         except Exception as e:
#             raise
#             #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))
#
#     if "log" in kwds and os.path.isfile(kwds["log"]):
#         return kwds["log"]
#     return out
