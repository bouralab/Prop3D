import os
import shutil

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess
    maxcluster_path = os.path.dirname(os.path.dirname(subprocess.check_output(["which", "maxcluster"])))

def run_maxcluster(*args, **kwds):
    work_dir = kwds.get("work_dir", None)
    docker = kwds.get("docker", True)
    job = kwds.get("job", None)

    if work_dir is None:
        work_dir = os.getcwd()

    file_kwds = ["log", "e", "p", "l", "R", "Rl", "Ru", "F", "M"]
    in_file_kwds = ["e", "p", "l", "F", "M"]
    parameters = ["-"+a for a in args]+["-{}={}".format(k,v) for k,v in \
        kwds.iteritems() if k not in file_kwds]
    file_parameters = {k:v for k, v in kwds.iteritems() if k in file_kwds}

    assert any(os.path.isfile(os.path.join(work_dir, f)) for f in ("new.html", "run.cns"))

    if docker and apiDockerCall is not None and job is not None:
        for k,f in file_parameters.iteritems():
            if k in in_file_kwds and not os.path.abspath(os.path.dirname(f)) == os.path.abspath(work_dir):
                shutil.copy(full, work_dir)
            parameters.append("-{}={}".format(k, os.path.basename(f)))

        oldcwd = os.getcwd()
        os.chdir(work_dir)
        try:
            out = apiDockerCall(job,
                          'edraizen/maxcluster:latest',
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters
                          )
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
            #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
            #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)
        os.chdir(oldcwd)
    else:
        args = [maxcluster_path]+parameters
        try:
            out = subprocess.check_output(args)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))
    return out

def get_centroid(file_list, work_dir=None, docker=True, job=None):
    log_file = os.path.splitext(file_list)[0]+".log"
    run_maxcluster(file_list=file_list, log=file_list+".log", work_dir=work_dir,
        docker=docker, job=job)

    parse_centroids = False
    best_centroid_size = None
    best_centroid_file = None

    with open(log_file) as log:
        for line in log:
            if not parse_centroids and line.startswith("INFO  : Centroids"):
                next(log)
                next(log)
            elif parse_centroids and line.startswith("INFO  : ="):
                break
            elif parse_centroids:
                line = line[8:].rstrip()
                cluster, centroid, size, spread, pdb_file = line.split()
                if best_centroid_size is None or size>best_centroid_size:
                    best_centroid_size = size
                    best_centroid_file = pdb_file

    return best_centroid_file
