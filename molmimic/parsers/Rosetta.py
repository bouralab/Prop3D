import os, sys
import time
import shutil
from joblib import Memory
from molmimic.generate_data.util import silence_stdout, silence_stderr, number_of_lines

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
	import subprocess

memory = Memory(verbose=0)

def run_minimize(pdb_file, min_type="lbfgs_armijo_nonmonotone", min_tolerance=0.001,
  ignore_zero_occupancy=False, parameters=None, work_dir=None):
    if work_dir is None:
        work_dir = os.getcwd()

    prefix = os.path.splitext(os.path.basename(pdb_file))[0]
    parameters = list(parameters) if isinstance(parameters, (list, tuple)) else []
    parameters += [
        "-run:min_type", min_type,
        "-run:min_tolerance", str(min_tolerance),
        "-ignore_zero_occupancy", str(ignore_zero_occupancy).lower()]

    if apiDockerCall is not None:
        if not os.path.abspath(os.path.dirname(pdb_file)) == os.path.abspath(work_dir):
            shutil.copy(pdb_file, work_dir)

        minimized_file = os.path.join("/data", "{}.pdb_0001.pdb".format(prefix))
        score_file = os.path.join("/data", "{}.sc".format(prefix))
        pdb_file_in = os.path.join("/data", "{}.pdb".format(prefix))

		parameters += [
            "-s", pdb_file_in,
            "-out:file:scorefile", score_file,
            "-out:path:pdb", "/data",
            "-out:path:score", "/data"]

        try:
            apiDockerCall(job,
                          image='edraizen/minimize:latest',
                          working_dir=work_dir,
                          parameters=parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            return run(pdb_file, min_type=min_type, min_tolerance=min_tolerance,
              ignore_zero_occupancy=ignore_zero_occupancy, parameters=parameters,
              work_dir=work_dir)

	else:
        work_dir = os.path.dirname(os.path.abspath(pdb_file))
        minimized_file = os.path.join(workdir, "{}.pdb_0001.pdb".format(prefix))
        score_file = os.path.join(work_dir, "{}.sc".format(prefix))

		parameters += [
            "-out:file:scorefile", score_file,
            "-out:path:pdb", work_dir,
            "-out:path:score", work_dir]

        command = ["minimize.static.linuxgccrelease"]+parameters
        command += [
            "-out:file:scorefile", score_file,
            "-out:path:pdb", work_dir,
            "-out:path:score", work_dir
        ]

        try:
            subprocess.check_output(["minimize.static.linuxgccrelease",
                "-s", pqr_file,
                "-run:min_type", "lbfgs_armijo_nonmonotone",
                "-run:min_tolerance", "0.001",
                "-overwrite", "false", #Pandas apply calls first row twice so this is needed
                "-ignore_zero_occupancy", "false",
                "-out:file:scorefile", score_file,
                "-out:path:pdb", pdb_path,
                "-out:path:score", pdb_path],
                stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            raise RuntimeError("Unable to minimize file {}".format(pqr_file))

    attempts = 0
    while not os.path.isfile(minimized_file) or number_of_lines(minimized_file) == 0:
        if attempts >= 10:
            raise RuntimeError("Unable to minimize file {}".format(pdb_file))
        time.sleep(1)
        attempts += 1

    shutil.move(minimized_file, minimized_file+".rosetta")
    return minimized_file+".rosetta", score_file
