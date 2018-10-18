import os, sys
import json
import shutil
import subprocess
from molmimic.generate_data.util import silence_stdout, silence_stderr

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

modeller_file = """import os, sys, json

print os.getcwd()
os.chdir("{work_dir}")
print os.getcwd()

#new_target = open(os.devnull, "w")
#old_target, sys.stdout = sys.stdout, new_target
#old_err, sys.stderr = sys.stderr, new_target

from modeller import *
from modeller.automodel import *

env = environ()
env.io.hydrogen = env.io.hetatm = env.io.water = True
env.io.atom_files_directory = ["{work_dir}"]
a = automodel(env,
          alnfile  = "{pir}",
          knowns   = ("{template}"),  # Use just one template.
          sequence = "{model}",
          assess_methods=(assess.DOPE))
a.starting_model= 1
a.ending_model = {num_models}
a.make()

#sys.stdout = old_target
#sys.stderr = old_err

output = {{x["name"]:x["DOPE score"] for x in a.outputs if x['failure'] is None}}
print json.dumps(output)
"""
pir_text = ">P1;{template}\nstructure:{template}:.:{chain}:.:{chain}::::\n{seq}*\n" + \
           ">P1;{target}\nsequence:{target}:.:{chain}:.:{chain}::::\n{seq}*\n"

def run_modeller(pir, template, model, num_models=5, work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    python_file = os.path.join(work_dir, "run_modeller.py")

    if docker and apiDockerCall is not None and job is not None:
        #Docker can only read from work_dir
        if not os.path.abspath(os.path.dirname(pir)) == os.path.abspath(work_dir):
            shutil.copy(pir, work_dir)

        if not os.path.abspath(os.path.dirname(template)) == os.path.abspath(work_dir):
            shutil.copy(template, work_dir)

        with open(python_file, "w") as f:
            f.write(modeller_file.format(
                pir = os.path.join("/data", os.path.basename(pir)),
                template = os.path.basename(template).rsplit(".", 1)[0],
                model=model, num_models=num_models, work_dir="/data"))
        print open(python_file).read()
        parameters = [os.path.join("/data", os.path.basename(python_file))]

        try:
            outputs = apiDockerCall(job,
                          image='edraizen/modeller:latest',
                          working_dir=work_dir,
                          parameters=parameters)
            job.log(outputs)
            outputs = outputs.splitlines()[-1]
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
            #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
            #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)

    else:
        with open(python_file, "w") as f:
            f.write(modeller_file.format(
                pir = pir, template = template, model=model, num_models=num_models, work_dir=work_dir))

        try:
            outputs = subprocess.check_output([sys.executable, python_file])
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))
    outputs = json.loads(outputs)
    best_pdb, best_dope = min(outputs.iteritems(), key=lambda x: x[1])
    best_pdb = os.path.join(work_dir, best_pdb)
    assert os.path.isfile(best_pdb)
    for f in outputs.keys():
	path = os.path.join(work_dir, f)
        if path != best_pdb:
            try:
                os.remove(path)
            except OSError:
                pass
    return best_pdb

def run_ca2model(ca_only_model, chain, num_models=5, work_dir=None, docker=True, job=None):
    """Convert a CA only model to a predicted full atom structure using modeller
    """
    if work_dir is None:
        work_dir = os.getcwd()

    prefix = os.path.basename(ca_only_model).split(".", 1)[0]
    template_file_prefix = "{}_ca_model".format(prefix)
    model_file_prefix = "{}_full_model".format(prefix)
    template_file = "{}.pdb".format(template_file_prefix)
    full_template_file = os.path.join(work_dir, template_file)

    shutil.copyfile(ca_only_model, full_template_file)

    #Auto-scaling on AWS with toil has trouble finding modules? Heres the workaround
    PDB_TOOLS = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdb_tools")

    with open(full_template_file) as f:
        seq = "\n".join(subprocess.check_output([sys.executable, os.path.join(PDB_TOOLS,
            "pdb_toseq.py")], stdin=f).split("\n")[1:])

    pir_file = os.path.join(work_dir, template_file_prefix+".pir")
    with open(pir_file, "w") as pir:
        pir.write(pir_text.format(template = template_file_prefix, target = model_file_prefix, chain=chain, seq=seq))
    print open(pir_file).read()
    return run_modeller(pir_file, full_template_file, model_file_prefix, num_models=num_models,
        work_dir=work_dir, docker=docker, job=job)
