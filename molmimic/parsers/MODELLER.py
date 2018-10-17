import os, sys
import json
from molmimic.generate_data.util import silence_stdout, silence_stderr

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

modeller_file = """import os, sys, json

new_target = open(os.devnull, "w")
old_target, sys.stdout = sys.stdout, , new_target
old_err, sys.stderr = sys.stderr, new_target

from modeller import *
from modeller.automodel import *

env = environ()
env.io.hydrogen = env.io.hetatm = env.io.water = True
env.io.atom_files_directory = ["/data"]
a = automodel(env,
          alnfile  = "{pir}",
          knowns   = ("{template}"),  # Use just one template.
          sequence = "{model}",
          assess_methods=(assess.DOPE))
a.starting_model= 1
a.ending_model = {num_models}
a.make()

sys.stdout = old_target
sys.stderr = old_err

output = {x["name"]:x["DOPE score"] for x in a.outputs if x['failure'] is None}
print json.dumps(output)
"""

def run_modeller(pir, template, model, num_models=5, work_dir=None, job=None):
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
            python_file.write(modeller_file.format(
                pir = os.path.join("/data", os.path.basename(pir)),
                template = os.path.join("/data", os.path.basename(template)),
                model=model, num_models=num_models))

        parameters = [os.path.join("/data", os.path.basename(python_file))]

        try:
            with silence_stdout(), silence_stderr():
            outputs = apiDockerCall(job,
                          image='edraizen/modeller:latest',
                          working_dir=work_dir,
                          parameters=parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_scwrl(pdb_file, output_prefix=output_prefix, framefilename=framefilename,
            #    sequencefilename=sequencefilename, paramfilename=paramfilename, in_cystal=in_cystal,
            #    remove_hydrogens=remove_hydrogens, remove_h_n_term=remove_h_n_term, work_dir=work_dir, docker=False)

    else:
        with open(python_file, "w") as f:
            python_file.write(modeller_file.format(
                pir = pir, template = template, model=model, num_models=num_models))

        try:
            outputs = subprocess.check_output([sys.executable, python_file])
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    ouputs = json.loads(outputs)
    best_pdb, best_dope = min(outputs.iteritems(), key=lambda x: x[1])
    assert os.path.isfile(best_pdb)
    return best_pdb

def run_ca2model(ca_only_model, pdb, chain, num_models=5, work_dir=None, docker=True, job=None):
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
        seq = subprocess.check_output([sys.executable, os.path.join(PDB_TOOLS,
            "pdb_toseq.py")], stdin=f)

    pir_file = template_file_prefix+".pir"
    with open(pir_file) as pir:
        pir.write(">P1;{template}\nstructure:{template}:.:{chain}:.:{chain}::::\n{seq}*\n" + \
                  ">P1;{target}\nsequence:{target}:.:{chain}:.:{chain}::::\n{seq}*\n".format(
                  template = template_file_prefix,
                  target = model_file_prefix,
                  chain=chain, seq=seq))

    return run_modeller(pir_file, full_template_file, model_file_prefix, num_models=num_models,
        work_dir=work_dir, job=job)
