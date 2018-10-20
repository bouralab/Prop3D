##################################################
# File: CNS.py
# Author: Eli Draizen
# Date: 9-25-2018
##################################################

#Standard libraries
import os
import re
import shutil

from molmimic.generate_data.util import izip_missing

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

script_dir = os.path.dirname(__file__)

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

def CNS(input_file, prefix, work_dir=None, docker=True, job=None, template=True, **template_kwds):
    work_dir = work_dir or os.getcwd()
    if docker and apiDockerCall is not None and job is not None:
        if not template and not os.path.abspath(os.path.dirname(input_file)) == os.path.abspath(work_dir):
            shutil.copy(input_file, work_dir)
            inp = input_file
        else:
            updated_templates = {}
            for k, v in template_kwds.iteritems():
                if os.path.isfile(v) and not os.path.abspath(os.path.dirname(v)) == os.path.abspath(work_dir):
                    shutil.copy(v, work_dir)
                updated_templates[k] = os.path.join("/data", os.path.basename(v))
            inp = generate_input(input_file, prefix, work_dir, **updated_templates)

        try:
            parameters = [os.path.join("/data", os.path.basename(inp))]
	    output = apiDockerCall(job,
	        image='edraizen/cns:latest',
	        working_dir=work_dir,
		parameters=parameters)
	except (SystemExit, KeyboardInterrupt):
	    raise
        except:
            raise
            return CNS.run(input_file, work_dir=work_dir, docker=False, job=job)
    else:
        load_cns_environment()
        inp = generate_input(input_file, prefix, work_dir, **template_kwds)
        with open(input_file) as inp:
            output = subprocess.check_output(["cns"], stdin=inp)
    return output

def Minimize(pdb_file, output_file_prefix=None, work_dir=None, docker=True, job=None):
    """Minimize a single protein.

    Parameters
    -----------
    pdb_file : str
        Path to pdb to minimize

    Returns
    -------
    structure_file : str
        Path to minimed structure
    """
    work_dir = work_dir or os.getcwd()
    output_file_prefix = output_file_prefix or os.path.splitext(os.path.basename(pdb_file))[0]
    temp_minimization_file = os.path.join(script_dir, "CNS-Templates", "model_minimize.inp")
    minimized_file = os.path.join(work_dir, "{}.min".format(output_file_prefix))
    cns = CNS(temp_minimization_file, output_file_prefix, work_dir=work_dir, docker=docker, job=job, INSERT_PDB_HERE=pdb_file, INSERT_OUTPUT_FILE=minimized_file)

    if not os.path.isfile(minimized_file):
        raise RuntimeError("CNS minimization has failed!! Check output:\n{}".format(cns))

    m = re.search("%CREAD-ERR: residue ID and insertion character ([0-9a-zA-Z]+) exceed 4 characters.")
    if m:
        if len(m.group(1)) < 6:
            #Strange bug in CNS where resn is 4 and icode is 1, valid PDB but CNS chokes
            #CNS corrects them by removing with first character
            #Fix by copying the original PDB cols before the XYZ coordiantes
            new_minimized_file = os.path.join(work_dir, "{}.fixed.min".format(output_file_prefix))
            with open(pdb_file) as old, open(minimized_file) as minf, open(new_minimized_file, "w") as new:
                 for oline, mline in izip_missing(old, minf, fillvalue=None)
                    if oline is None or mline is None:
                        raise RuntimeError("CNS minimization has failed!! {}. Check output: {}".format(m.group(0), cns))
                    new.write(mline[0:22]+oline[22:27]+mline[27:])
            minimized_file = new_minimized_file
        else:
            raise RuntimeError("CNS minimization has failed!! {}. Check output: {}".format(m.group(0), cns))

    return minimized_file

def generate_input(template_file, prefix, work_dir, **kwds):
    """Minimize a single protein.

    Returns
    -------
    structure_file : str
        Path to minimed structure
    """
    input_file = os.path.join(work_dir, "{}_{}".format(prefix, os.path.basename(template_file)))

    with open(template_file) as temp, open(input_file, "w") as inp:
        for line in temp:
            for key, val in kwds.iteritems():
                line = re.sub(key, val, line)
            inp.write(line)

    return input_file
