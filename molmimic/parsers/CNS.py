##################################################
# File: CNS.py
# Author: Eli Draizen
# Date: 9-25-2018
##################################################

#Standard libraries
import os
import re

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
	import subprocess

def load_cns_environment():
    if "CNS_TOPPAR" not in os.environ:
        try:
            cns = subprocess.check_output(["which", "cns_solve"])
        except:
            raise RuntimeError("Unable to load CNS")
        cns_dir = os.path.dirname(os.path.dirname(os.path.dirname(cns)))
        env_file = os.path.join(cns_dir, "cns_solve_env.zsh")
        output = subprocess.check_output("source {}; env -0".format(env_file),
            shell=True, executable="/bin/zsh")
        new_env = dict([line.partition('=')[::2] for line in output.split('\0')])
        os.environ.clear()
        os.environ.update(new_env)

class CNS(object):
    def __init__(self, input_file, work_dir=None, docker=True, job=None):
    	if docker and apiDockerCall is not None:
            if not os.path.abspath(os.path.dirname(input_file)) == os.path.abspath(work_dir):
                shutil.copy(input_file, work_dir)

            try:
                parameters = [input_file]
                apiDockerCall(job,
                              image='edraizen/cns:latest',
                              working_dir=work_dir,
                              parameters=parameters)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                return self.__class__(full_pqr_file, input_file=input_file,
                    work_dir=work_dir, docker=False)
        else:
            load_cns_environment()
            with open(input_file) as inp:
                return subprocess.check_output(["cns"], stdin=inp)

class Minimize(CNS):
    """"Minimize a single protein.

    Parameters
    -----------
    pdb_file : str
        Path to pdb to minimize

    Returns
    -------
    structure_file : str
        Path to minimed structure
    """
    def __init__(self, pdb_file, output_file_prefix=None, work_dir=None, docker=True, job=None):
        """Create the CNS input file"""
        self.pdb_file = pdb_file
        self.output_file_prefix = output_file_prefix or os.path.splitext(os.path.basename(pdb_file)))
        self.work_dir = work_dir or os.getcwd()
        self.docker = docker
        input_file, minimized_file = self.__generate_minimization_input()
        super(Minimize, self).__init__(input_file, work_dir=work_dir, docker=docker, job=job)
        return minimized_file

    def __generate_minimization_input(self):
        """"Minimize a single protein.

        Returns
        -------
        structure_file : str
            Path to minimed structure
        """

        temp_minimization_file = os.path.join(script_dir, "CNS-Templates", "model_minimize.inp")
        model_minimization_base = "{}_model_minimize.inp".format(self.output_file_prefix)
        model_minimization_base_file = os.path.join(self.work_dir, model_minimization_base)

        pdb_file_base = os.path.basename(self.pdb_file)
        output_file_base = "{}.min".format(self.output_file_prefix)
        minimized_file = os.path.join(self.work_dir, output_file_base)

        if docker:
            if os.path.abspath(os.path.dirname(self.pdb_file)) == os.path.abspath(self.work_dir):
                shutil.copy(self.pdb_file, self.work_dir)

            pdb_file = os.path.join("/data", pdb_file_base)
            minimized_file = os.path.join("/data", pdb_file_base)
            output_file = os.path.join("/data", output_file_base)
        else:
            pdb_file = self.pdb_file
            output_file = minimized_file

        with open(model_minimization_file, "w") as inp, open(model_minimization_file, "w") as temp:
            for line in temp:
                line = re.sub(r'INSERT_PDB_HERE', pdb_file, line)
                line = re.sub(r'INSERT_OUTPUT_FILE', output_file, line)
                inp.write(line)

        return minimized_file
