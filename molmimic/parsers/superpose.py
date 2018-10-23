try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess
    haddock_path = os.path.dirname(os.path.dirname(subprocess.check_output(["which", "RunHaddock.py"])))

from molmimic.generate_data.util import PDB_TOOLS

def align(fixed_file, fixed_chain, moving_file, moving_chain, method="tmalign", work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    outFile = os.path.join(work_dir, "aligned.pdb")

    if method == "tmalign":
        image = "edraizen/tmalign:latest"
        parameters = lambda f, m: [m, f, "-o", outfile+".sup"]
    elif method == "ce":
        image = "edraizen/ce:latest:
        parameters = lambda f, m, o: ["--file1", f, "--file2", m, "-outputPDB", "-outFile", outfile+".sup"])

    if docker and apiDockerCall is not None and job is not None:
        try:
            apiDockerCall(job,
                          image=image,
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters(os.path.basename(fixed_file), os.path.basename(moving_file))
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise

    else:
        raise RuntimeError("Only docker works at the moment")

    with open(outFile, "w") as out:
        subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"),
            "-{}".format(moving_chain), outfile+".sup"], stdout=out)

    os.remove(outfile+".sup")

    assert os.path.isfile(outFile)
    return outFile
