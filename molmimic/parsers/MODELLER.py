import os, sys
import json
import shutil
import subprocess
from molmimic.util import silence_stdout, silence_stderr

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    import subprocess

modeller_file = """import os, sys, json

os.chdir("{work_dir}")

from modeller import *
from modeller.automodel import *

#Add New Class Here

env = environ()
env.io.hydrogen = env.io.hetatm = env.io.water = True
env.io.atom_files_directory = ["{work_dir}"]
a = automodel(env,
          alnfile  = "{pir}",
          knowns   = ("{template}"),  # Use just one template.
          sequence = "{model}",
          assess_methods=(assess.DOPE))
a.starting_model = 1
a.ending_model = {num_models}
a.make()

output = {{x["name"]:x["DOPE score"] for x in a.outputs if x['failure'] is None}}
print json.dumps(output)
"""

pir_text = ">P1;{template}\nstructure:{template}:.:{chain}:.:{chain}::::\n{seq}*\n" + \
           ">P1;{target}\nsequence:{target}:.:{chain}:.:{chain}::::\n{seq}*\n"


ca_model = """# This will transfer residue numbers and chain ids from model2 to model.
class CaModel(automodel):
    def special_patches(self, aln):
        if not hasattr(self, "ca_template"):
            self.ca_template = model(env, file="{template}")
        self.res_num_from(self.ca_template, aln)
"""

class MODELLER(Container):
    IMAGE = "docker://edraizen/modeller:latest"
    LOCAL = [sys.executable]
    PARAMETERS = [("modeller_file", "path:in")]

    def __call__(self, *args, **kwds):
        pass

    def automodel(self, template_pdb, model_sequence, pir_file, num_models=5,
      make_model_callback=None):
        python_file = os.path.join(self.work_dir, "run_modeller.py")

        with open(python_file, "w") as f:
            f.write(make_modeller_file(
                pir = os.path.join("/data", os.path.basename(pir)),
                template = os.path.basename(template).rsplit(".", 1)[0],
                model=model, num_models=num_models, work_dir="/data",
                callback=make_model_callback))

    def make_ca_model(self, template_pdb, chain, num_models=5):
        assert is_ca_model(template_pdb)

        prefix = os.path.basename(template_pdb).split(".", 1)[0]
        template_file_prefix = "{}_ca_model".format(prefix)
        model_file_prefix = "{}_full_model".format(prefix)
        template_file = "{}.pdb".format(template_file_prefix)

        with open(template_pdb) as f:
            seq = "\n".join(subprocess.check_output(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_toseq.py")],
                stdin=f).decode("utf-8").split("\n")[1:])

        pir_file = self._make_pir(template_file_prefix, model_file_prefix,
            chain, seq)

        def ca_callback(modeller_file, pir, template, model, num_models, work_dir):
            _ca_model = ca_model.format(template=template)
            f = modeller_file.replace("#Add New Class Here", _ca_model)
            return f.replace("a = automodel", "a = CaModel")

        return self.automodel(pir_file, template_file, model_file_prefix,
            num_models=num_models, make_model_callback=ca_callback)

    def _make_pir(self, template_file_prefix, model_file_prefix, chain, seq):
        pir_file = os.path.join(self.work_dir, template_file_prefix+".pir")
        with open(pir_file, "w") as pir:
            pir.write(pir_text.format(
            template = template_file_prefix,
            target = model_file_prefix,
            chain=chain,
            seq=seq))
        return pir_file






def make_modeller_file(pir, template, model, num_models=5, work_dir="", callback=None):
    global _modeller_file
    _modeller_file = modeller_file.format(pir=pir, template=template, model=model,
        num_models=num_models, work_dir=work_dir)
    if callable(callback):
        _modeller_file = callback(_modeller_file, pir, template, model, num_models, work_dir)
    return _modeller_file



def run_modeller(pir, template, model, num_models=5, work_dir=None, docker=True,
  job=None, make_model_callback=None):
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
            f.write(make_modeller_file(
                pir = os.path.join("/data", os.path.basename(pir)),
                template = os.path.basename(template).rsplit(".", 1)[0],
                model=model, num_models=num_models, work_dir="/data",
                callback=make_model_callback))

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
            f.write(make_modeller_file(
                pir = pir, template = template, model=model, num_models=num_models,
                work_dir=work_dir, callback=make_model_callback))

        try:
            outputs = subprocess.check_output([sys.executable, python_file])
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))
    outputs = json.loads(outputs)
    best_pdb, best_dope = min(iter(list(outputs.items())), key=lambda x: x[1])
    best_pdb = os.path.join(work_dir, best_pdb)
    assert os.path.isfile(best_pdb)
    for f in list(outputs.keys()):
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
        pir.write(pir_text.format(
        template = template_file_prefix,
        target = model_file_prefix,
        chain=chain,
        seq=seq))

    return run_modeller(pir_file, full_template_file, model_file_prefix, num_models=num_models,
        work_dir=work_dir, docker=docker, job=job, make_model_callback=ca_callback)
