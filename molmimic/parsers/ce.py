try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess
    haddock_path = os.path.dirname(os.path.dirname(subprocess.check_output(["which", "RunHaddock.py"])))

def align(fixed_file, moving_file, work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    outFile = os.path.join(work_dir, "aligned.pdb")

    if docker and apiDockerCall is not None and job is not None:
        try:
            apiDockerCall(job,
                          image='edraizen/ce:latest',
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=[
                            "--file1", os.path.basename(fixed_file),
                            "--file2", os.path.basename(moving_file),
                            "-outputPDB",
                            "-outFile", os.path.basename(outFile)])
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise

    else:
        raise RuntimeError("Only docker works at the moment")

    assert os.path.isfile(outFile)
    return outFile
