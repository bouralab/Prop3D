import os, sys
import glob
import json
import shutil
import subprocess

from Prop3D.parsers.container import Container
from Prop3D.util import SubprocessChain, safe_remove, silence_stdout, silence_stderr
from Prop3D.util.pdb import is_ca_model, get_first_chain, get_all_chains, \
    PDB_TOOLS
from Prop3D.common.ProteinTables import three_to_one

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
    PARAMETERS = [("modeller_file", "path:in", "")]
    ENTRYPOINT = "python"

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
        best_pdb = os.path.join(self.work_dir, os.path.basename(best_model["name"]))
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
        elif hasattr(self, "running_modloop") and self.running_modloop:
            modeller_file = self.make_loop_modeller_file()
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

    def mutate_structure(self, template_pdb, mutants, chain=None, num_models=5, return_best=True, clean=True):
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

    def model_loops(self, template_pdb, loops, ss_restraints=None, sequences=None, chain=None, num_models=5, return_best=True, clean=True):
        if chain is None:
            chain = get_first_chain(template_pdb)

        if chain is None:
            raise ValueError("Invalid PDB")

        self.template_pdb = template_pdb
        self.loops = loops
        self.ss_restraints = ss_restraints
        self.num_models = num_models
        self.prefix = os.path.splitext(os.path.basename(template_pdb))[0]
        self.target_id = f"{self.prefix}_target"


        self.pir_file = None
        if sequences is not None:
            template_pdb = f"{self.prefix}_template"
            shutil.copy(self.template_pdb, f"{template_pdb}.pdb")
            self.pir_file = self.make_pir(template_pdb, chain, sequences[0],
                self.target_id, chain, sequences[1])
            self.template_pdb = template_pdb


        modeller_file = self.make_loop_modeller_file()

        self.running_modloop = True
        with silence_stdout(), silence_stderr():
            outputs = self(modeller_file=modeller_file)
        self.running_modloop = False

        #target_id = os.path.splitext(os.path.basename(template_pdb))[0]

        try:
            with open(os.path.join(self.work_dir, "{}.dope_scores".format(self.target_id))) as fh:
                dope_scores = json.load(fh)
        except (IOError, FileNotFoundError):
            raise RuntimeError("Error running modeller, no dope score file found")

        if return_best:
            output = self.select_best_model(dope_scores)

            files_to_remove = []
            clean = int(clean)
            if clean>0:
                files_to_remove.append(modeller_file)
                files_to_remove += glob.glob(os.path.join(self.work_dir,
                    "{}.*".format(self.target_id)))
            if clean>1:
                files_to_remove += glob.glob(os.path.join(self.work_dir,
                        "{}.*".format(self.template_pdb)))
            files_to_remove = set(files_to_remove)
            files_to_remove -= set([output])
        else:
            output = dope_scores

        #safe_remove(files_to_remove)

        return output

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
        work_dir = os.path.dirname(template_pdb)

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

    def make_loop_modeller_file(self):
        """Modified from modloop"""
        residue_range = []
        for start, stop in self.loops:
            residue_range.append("           self.residue_range"
                                 "('%s:%s', '%s:%s')," % (*start, *stop))
        residue_range = "\n".join(residue_range)

        alpha_range = []
        beta_range = []
        if self.ss_restraints is not None:
            for ss_type, start, stop in self.ss_restraints:
                if ss_type == "E":
                    beta_range.append("            self.residue_range"
                                "('%s:%s', '%s:%s')," % (*start, *stop))
                else:
                    alpha_range.append("           self.residue_range"
                                             "('%s:%s', '%s:%s')," % (*start, *stop))
        alpha_range = "\n".join(alpha_range)
        beta_range = "\n".join(beta_range)

        template_pdb = self.format_in_path(None, self.template_pdb)
        target_id = self.target_id
        modeller_file = os.path.join(self.work_dir, "run_modloop_{}.py".format(
            self.prefix))
        work_dir = os.path.dirname(template_pdb) #self.work_dir
        template_pdb = os.path.splitext(os.path.basename(template_pdb))[0]

        if self.pir_file is not None:
            pir_file = self.format_in_path(None, self.pir_file)

        with open(modeller_file, "w") as f:
            f.write("""
# Run this script with something like
#    python loop.py N > N.log
# where N is an integer from 1 to the number of models.
#
# ModLoop does this for N from 1 to 300 (it runs the tasks in parallel on a
# compute cluster), then returns the single model with the best (lowest)
# value of the Modeller objective function.

import sys, os, json

# to get different starting models for each task
taskid = int(sys.argv[1]) if len(sys.argv) > 1 else 5

os.chdir('{work_dir}')

from modeller import *
from modeller.automodel import *

env = environ()
env.io.hydrogen = env.io.hetatm = env.io.water = True
env.io.atom_files_directory = ['{work_dir}', os.getcwd()]

class MyLoop(loopmodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=[c.name for c in aln[self.sequence].chains])

    def select_loop_atoms(self):
        rngs = (
            {residue_range}
        )
        for rng in rngs:
            if len(rng) > 30:
                raise ModellerError("loop too long")
            s = selection(rngs)
            if len(s.only_no_topology()) > 0:
                raise ModellerError("some selected residues have no topology")
            return s

    def special_restraints(self, aln):
        rsr = self.restraints

        for a in (
            {alpha_range}
        ):
            rsr.add(secondary_structure.alpha(a))
        for b in (
            {beta_range}
        ):
            rsr.add(secondary_structure.strand(b))
""".format(**locals()))
            if self.pir_file is None:
                f.write("""
m = MyLoop(env, inimodel='{template_pdb}',
    sequence='{target_id}',
    assess_methods=(assess.DOPE),
    loop_assess_methods=(assess.DOPE)
)
""".format(**locals()))
            else:
                f.write("""
m = MyLoop(env, alnfile = '{pir_file}',
    knowns = '{template_pdb}',
    sequence = '{target_id}',
    assess_methods=(assess.DOPE),
    loop_assess_methods=(assess.DOPE)
)
m.starting_model= 1
m.ending_model  = 1
""".format(**locals()))

            f.write("""
m.loop.md_level = refine.slow
m.loop.starting_model = 1
m.loop.ending_model = taskid
m.make()

results = {{"scores":[], "errors":[]}}
results["scores"] = [{{
    "name": x["name"],
    "score":x["DOPE score"]}} for x in m.loop.outputs if x['failure'] is None]
results["errors"] = [{{
    "name": x["name"],
    "score":str(x["failure"])}} for x in m.loop.outputs if x['failure'] is not None]
if len(results["scores"]) == 0:
    del results["scores"]
if len(results["errors"]) == 0:
    del results["errors"]

with open("{target_id}.dope_scores", "w") as fh:
    json.dump(results, fh)
""".format(**locals()))
        return modeller_file

# # Run this script with something like
# #    python loop.py N > N.log
# # where N is an integer from 1 to the number of models.
# #
# # ModLoop does this for N from 1 to 300 (it runs the tasks in parallel on a
# # compute cluster), then returns the single model with the best (lowest)
# # value of the Modeller objective function.
#
# import sys, os, json
#
# # to get different starting models for each task
# taskid = int(sys.argv[1]) if len(sys.argv) > 1 else 5
#
# os.chdir('{work_dir}')
#
# from modeller import *
# from modeller.automodel import *
#
# env = environ()
# env.io.hydrogen = env.io.hetatm = env.io.water = True
# env.io.atom_files_directory = ['{work_dir}', os.getcwd()]
#
# class MyLoop(loopmodel):
#     def special_patches(self, aln):
#         # Rename both chains and renumber the residues in each
#         self.rename_segments(segment_ids=[c.name for c in aln[self.sequence].chains])
#
#     def fake_select_atoms(self):
#         rngs = [
# {alpha_range}
#         ] + [
# {beta_range}
#         ]
#         s = selection(rngs)
#         if len(s.only_no_topology()) > 0:
#             raise ModellerError("some selected residues have no topology")
#         return s
#
#     def special_restraints(self, aln):
#         rsr = self.restraints
#         for a in (
# {alpha_range}
#         ):
#             rsr.add(secondary_structure.alpha(a))
#         for b in (
# {beta_range}
#         ):
#             rsr.add(secondary_structure.strand(b))
# """.format(**locals()))
#             if self.pir_file is None:
#                 f.write("""
# m = MyLoop(env, inimodel='{template_pdb}',
#     sequence='{target_id}',
#     assess_methods=(assess.DOPE),
# #    loop_assess_methods=(assess.DOPE)
# )
# """.format(**locals()))
#             else:
#                 f.write("""
# m = MyLoop(env, alnfile = '{pir_file}',
#     knowns = '{template_pdb}',
#     sequence = '{target_id}',
#     assess_methods=(assess.DOPE),
#     #loop_assess_methods=(assess.DOPE)
# )
# m.starting_model= 1
# m.ending_model  = 1
# """.format(**locals()))
#
#             f.write("""
# m.loop.md_level = refine.slow
# m.loop.starting_model = 1
# m.loop.ending_model = taskid
# m.very_fast = True
# m.md_level = refine.very_fast
# m.make(exit_stage=2)
#
# results = {{"scores":[], "errors":[]}}
# results["scores"] = [{{
#     "name": x["name"],
#     "score":x["DOPE score"]}} for x in m.loop.outputs if x['failure'] is None]
# results["errors"] = [{{
#     "name": x["name"],
#     "score":str(x["failure"])}} for x in m.loop.outputs if x['failure'] is not None]
# if len(results["scores"]) == 0:
#     del results["scores"]
# if len(results["errors"]) == 0:
#     del results["errors"]
#
# with open("{target_id}.dope_scores", "w") as fh:
#     json.dump(results, fh)
# """.format(**locals()))
#         return modeller_file
