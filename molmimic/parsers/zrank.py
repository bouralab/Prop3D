import os
import tempfile

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess
    zrank_path = subprocess.check_output(["which", "zrank"])

def run_zrank(complex_path, refinement=False, work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    _parameters = ["-R"] if refinement else []

    if not isinstance(complex_path, (list, tuple)):
        complex_path = [complex_path]

    listfile = tempfile.NamedTemporaryFile(dir=work_dir, prefix="listfile", suffix=".txt", delete=False)
    for pdb in complex_path:
        print >> listfile, os.path.basename(pdb)
    listfile.close()

    with open(listfile.name) as f:
        print f.read()

    if docker and apiDockerCall is not None and job is not None:
        parameters = _parameters + [os.path.basename(listfile.name)]
        try:
            out = apiDockerCall(job,
                          image='edraizen/zrank:latest',
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters)
            print out
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            print "Error running docker for {} becuase {}".format(complex_path, e)
            raise
            return run_zrank(complex_path, refinement=refinement, work_dir=work_dir, docker=False)
    else:
        cmd = [zrank_path] + _parameters + [os.path.join("/data", listfile.name)]
        try:
            subprocess.call(cmd)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise RuntimeError("Cannot run zrank for {}. Error: {}".format(complex_path, e))

    assert os.path.isfile(listfile.name+".zr.out"), "No output for zrank"

    with open(listfile.name+".zr.out") as f:
        scores = dict(line.rstrip().split() for line in f)

    if len(complex_path) == 1:
        scores = scores.values()[0]

    for f in (listfile.name, listfile.name+".zr.out"):
        try:
            os.remove(f)
        except OSError:
            pass

    return scores
