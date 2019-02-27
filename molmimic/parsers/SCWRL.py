import os, sys
import subprocess
import shutil

from joblib import Memory

from molmimic.util import silence_stdout, silence_stderr
from molmimic.parsers.psize import Psize

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

def run_scwrl(pdb_file, output_prefix=None, framefilename=None, sequencefilename=None,
  paramfilename=None, in_cystal=False, remove_hydrogens=False, remove_h_n_term=False,
  work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    full_pdb_path = pdb_file
    pdb_path = os.path.basename(full_pdb_path)
    output_file = "{}.scwrl".format(pdb_path)

    if output_prefix is None:
        output_prefix = os.path.splitext(full_pdb_path)[0]

    output_file = "{}.scwrl.pdb".format(output_prefix)

    _parameters = [p for p, use in [("-#", in_cystal), ("-h", remove_hydrogens), \
        ("-t", remove_h_n_term)] if use]

    if docker and apiDockerCall is not None and job is not None:
        #Docker can only read from work_dir
        if not os.path.abspath(os.path.dirname(pdb_file)) == os.path.abspath(work_dir):
            shutil.copy(pdb_file, work_dir)

        parameters = ["-i", "/data/{}".format(os.path.basename(pdb_file))]
        parameters += ["-o", "/data/{}".format(os.path.basename(output_file))]
        parameters += _parameters

        if framefilename is not None and os.path.isfile(framefilename):
            if not os.path.abspath(os.path.dirname(framefilename)) == os.path.abspath(work_dir):
                shutil.copy(framefilename, work_dir)
            parameters += ["-f", "/data/{}".format(os.path.basename(framefilename))]

        if sequencefilename is not None and os.path.isfile(sequencefilename):
            if not os.path.abspath(os.path.dirname(sequencefilename)) == os.path.abspath(work_dir):
                shutil.copy(sequencefilename, work_dir)
            parameters += ["-s", "/data/{}".format(os.path.basename(sequencefilename))]

        if paramfilename is not None and os.path.isfile(paramfilename):
            if not os.path.abspath(os.path.dirname(paramfilename)) == os.path.abspath(work_dir):
                shutil.copy(paramfilename, work_dir)
            parameters += ["-p", "/data/{}".format(os.path.basename(paramfilename))]

        try:
            apiDockerCall(job,
                          image='edraizen/scwrl4:latest',
                          working_dir=work_dir,
                          parameters=parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
            #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
            #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)

        output_file = os.path.join(work_dir, os.path.basename(output_file))
        #if not os.path.abspath(os.path.dirname(output_file)) == os.path.abspath(work_dir):
        #    shutil.move(os.path.join(work_dir, os.path.basename(output_file)),
        #        os.path.abspath(os.path.dirname(output_file)))

    else:
        parameters = ["scwrl4", "-i", pdb_file, "-o", output_file]+_parameters
        if framefilename is not None and os.path.isfile(framefilename):
            parameters += ["-f", framefilename]
        if sequencefilename is not None and os.path.isfile(sequencefilename):
            parameters += ["-s", sequencefilename]
        if paramfilename is not None and os.path.isfile(paramfilename):
            parameters += ["-p", paramfilename]

        try:
            subprocess.call(parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    assert os.path.isfile(output_file)
    return output_file
