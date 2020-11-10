import os, sys
import glob
import json
import shutil
import subprocess

from molmimic.parsers.container import Container
from molmimic.util import SubprocessChain, safe_remove, silence_stdout, silence_stderr
from molmimic.util.pdb import is_ca_model, get_first_chain, get_all_chains, \
    PDB_TOOLS
from molmimic.common.ProteinTables import three_to_one

modeller_text = """import os, sys, json

os.chdir("{work_dir}")

from modeller import *
from modeller.automodel import *

#Add New Class Here

env = environ()
env.io.hydrogen = env.io.hetatm = env.io.water = True
env.io.atom_files_directory = ["{work_dir}", os.getcwd()]
a = {automodel}(env,
          alnfile  = "{pir_file}",
          knowns   = {template_id},
          sequence = "{target_id}",
          assess_methods=(assess.DOPE))
a.starting_model = 1
a.ending_model = {num_models}
a.make()

results = {{"scores":[], "errors":[]}}
results["scores"] = [{{"name": x["name"], "score":x["DOPE score"]}} for x in \
    a.outputs if x['failure'] is None]
results["errors"] = [{{"name": x["name"], "score":str(x["failure"])}} for x in \
    a.outputs if x['failure'] is not None]
if len(results["scores"]) == 0:
    del results["scores"]
if len(results["errors"]) == 0:
    del results["errors"]

with open("{target_id}.dope_scores", "w") as fh:
    json.dump(results, fh)
"""

pir_text = """>P1;{template_id}
structure:{template_id}:.:{template_chain}:.:{template_chain}::::
{template_seq}*
>P1;{target_id}
sequence:{target_id}:.:{target_chain}:.:{target_chain}::::
{target_seq}*
"""


ca_model_text = """# This will transfer residue numbers and chain ids from model2 to model.
class CaModel(automodel):
    def special_patches(self, aln):
        if not hasattr(self, "ca_template"):
            self.ca_template = model(env, file="{template_pdb}")
        self.res_num_from(self.ca_template, aln)
"""

class MODELLER(Container):
    IMAGE = "docker://edraizen/modeller:latest"
    LOCAL = [sys.executable]
    PARAMETERS = [("modeller_file", "path:in")]

    def automodel(self, template_id, target_id, pir_file, num_models=5,
      extra_modeller_code=None, automodel_command="automodel", return_best=True,
      clean=1):
        """
        Parameters
        ----------
        template_id : str, list, or tuple
            Name of pdb file name(s) to be used as the template. Do not include
            full path or extenison. This name(s) must in one of the PIR structures.
        target_id : str
            Name of the new modelled pdb file, only the basename of the file
            without extension.
        pir_file : str
            Full Path to the PIR file.
        num_models : int
            Number of models to create. Best one can be inferred using the DOPE
            score.
        extra_modeller_code: str or None
            Addition code to add to the python modeller script.

        Returns
        -------
        The path to the best scoring model if return_best is True, else a
        dictionary containing the scores or errors of each model.
        """
        self.template_id = template_id if isinstance(template_id, (list, tuple)) \
            else [template_id]
        self.target_id = target_id
        self.pir_file = pir_file
        self.num_models = num_models
        self.extra_modeller_code = extra_modeller_code
        self.automodel_command = automodel_command

        modeller_file = self.make_modeller_file()

        self.running_automodel = True
        with silence_stdout(), silence_stderr():
            outputs = self(modeller_file=modeller_file)
        self.running_automodel = False

        try:
            with open(os.path.join(self.work_dir, "{}.dope_scores".format(target_id))) as fh:
                dope_scores = json.load(fh)
        except (IOError, FileNotFoundError):
            raise RuntimeError("Error running modeller, no dope score file found")

        files_to_remove = []
        if clean>0:
            files_to_remove.append(modeller_file)
            files_to_remove += glob.glob(os.path.join(self.work_dir,
                "{}.*".format(target_id)))
        if clean>1:
            files_to_remove.append(self.pir_file)
            for t in self.template_id:
                files_to_remove += glob.glob(os.path.join(self.work_dir,
                    "{}.*".format(t)))
        files_to_remove = set(files_to_remove)

        if return_best:
            output = self.select_best_model(dope_scores)
            files_to_remove -= set([output])
        else:
            output = dope_scores

        safe_remove(files_to_remove)

        return output

    def select_best_model(self, results):
        if len(results["scores"]) == 0:
            return None

        best_model = min(results["scores"], key=lambda x: x["score"])
        best_pdb = os.path.join(self.work_dir, best_model["name"])
        assert os.path.isfile(best_pdb)

        #Remove other lower scoring models
        for model in results["scores"]:
            if model["name"] == best_model["name"]:
                continue

            try:
                os.remove(os.path.join(self.work_dir, model["name"]))
            except OSError:
                pass

        return best_pdb

    def local(self, *args, **kwds):
        self.is_local = True
        if hasattr(self, "running_automodel") and self.running_automodel:
            modeller_file = self.make_modeller_file()
            if len(args) >= 1:
                args = [modeller_file]+args[1:]
            kwds["modeller_file"] = modeller_file
        super().local(*args, **kwds)

    def make_ca_model(self, template_pdb, chain=None, num_models=5, return_best=True):
        return self.remodel_structure(template_pdb, chain=chain,
            num_models=num_models, return_best=return_best)

    def remodel_structure(self, template_pdb, chain=None, num_models=5, return_best=True, clean=True):
        if chain is None:
            chain = get_first_chain(template_pdb)

        if chain is None:
            raise ValueError("Invalid PDB")

        prefix = os.path.splitext(os.path.basename(template_pdb))[0]
        target_id = "{}_full_model".format(prefix)

        #Rename template to specify CA model
        template_id = "{}_ca_model".format(prefix)
        full_template_file = os.path.join(self.work_dir, "{}.pdb".format(template_id))
        shutil.copyfile(template_pdb, full_template_file)

        if len(get_all_chains(full_template_file)) > 1:
            cmd = [
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), \
                    "-{}".format(chain), full_template_file],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_toseq.py")]]
        else:
            cmd = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_toseq.py"),
                full_template_file]]

        seq = "\n".join(SubprocessChain(cmd)[0].decode("utf-8").split("\n")[1:]).strip()

        pir_file = self.make_pir(template_id, chain, seq, target_id, chain, seq)

        return self.automodel(template_id, target_id, pir_file,
            num_models=num_models, extra_modeller_code=ca_model_text,
            automodel_command="CaModel", return_best=return_best, clean=2 if clean else 0)

    def mutate_structure(template_pdb, mutants, chain=None, num_models=5, return_best=True, clean=True):
        if chain is None:
            chain = get_first_chain(template_pdb)

        if chain is None:
            raise ValueError("Invalid PDB")

        prefix = os.path.splitext(os.path.basename(template_pdb))[0]
        target_id = "{}_full_model".format(prefix)

        #Rename template to specify CA model
        template_id = "{}_ca_model".format(prefix)
        full_template_file = os.path.join(self.work_dir, "{}.pdb".format(template_id))

        if len(get_all_chains(full_template_file)) > 1:
            with open(full_template_file, "w") as fh:
                SubprocessChain([[sys.executable, os.path.join(PDB_TOOLS,
                    "pdb_selchain.py"), "-{}".format(chain), full_template_file]],
                    fh)
        else:
            shutil.copyfile(template_pdb, full_template_file)

        original_seq = ""
        mutant_seq = ""
        for resn, resi in get_pdb_residues(full_template_file, include_resn=True, use_raw_resi=True):
            original_seq += three_to_one(resi)
            if resi in mutants:
                mutant_seq += mutants[resi] if len(mutants[resi]) == 1 else three_to_one(mutants[resi])
            else:
                mutant_seq += three_to_one(resi)

        pir_file = self.make_pir(template_id, chain, original_seq, target_id,
            chain, mutant_seq)

        return self.automodel(template_id, target_id, pir_file,
            num_models=num_models, extra_modeller_code=ca_model_text,
            automodel_command="CaModel", return_best=return_best, clean=2 if clean else 0)

    def make_pir(self, template_id, template_chain, template_seq, target_id,
      target_chain, target_seq):
        pir_file = os.path.join(self.work_dir, template_id+".pir")
        with open(pir_file, "w") as pir:
            pir.write(pir_text.format(
                template_id    = template_id,
                template_chain = template_chain,
                template_seq   = template_seq,
                target_id      = target_id,
                target_chain   = target_chain,
                target_seq     = target_seq))
        return pir_file

    def make_modeller_file(self, template_id=None, target_id=None, pir_file=None,
      num_models=5, extra_modeller_code="", automodel_command="automodel"):
        """Make container specific paths written in modeller file"""
        if pir_file is None:
            template_id = self.template_id
            target_id = self.target_id
            pir_file = self.pir_file
            num_models = self.num_models
            extra_modeller_code = self.extra_modeller_code
            automodel_command = self.automodel_command
            work_dir = self.work_dir
        else:
            work_dir = os.path.dirname(pir_path)

        pir_file = self.format_in_path(None, pir_file)
        template_pdb = self.format_in_path(None, "{}.pdb".format(template_id[0]))

        modeller_file = os.path.join(self.work_dir, "run_modeller_{}.py".format(
            target_id))

        _modeller_text = modeller_text.replace("#Add New Class Here", extra_modeller_code)

        with open(modeller_file, "w") as f:
            _modeller_text = _modeller_text.format(
                template_id  = str(template_id),
                template_pdb = template_pdb,
                target_id    = target_id,
                pir_file     = pir_file,
                num_models   = num_models,
                work_dir     = work_dir,
                automodel    = automodel_command)

            f.write(_modeller_text)

        return modeller_file
