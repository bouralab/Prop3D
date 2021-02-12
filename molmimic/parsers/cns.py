##################################################
# File: CNS.py
# Author: Eli Draizen
# Date: 9-25-2018
##################################################

#Standard libraries
import os
import sys
import re
import shutil

import pandas as pd

from molmimic.parsers.container import Container
from molmimic.util import izip_missing, SubprocessChain
from molmimic.util.pdb import PDB_TOOLS

script_dir = os.path.dirname(__file__)

class CNS(Container):
    IMAGE = 'edraizen/cns:latest'
    LOCAL = ["cns_solve"]
    PARAMETERS = [("input_file", "path:in:stdin", "")]

    def _generate_input(self, template_file, prefix, **kwds):
        """Minimize a single protein.

        Returns
        -------
        structure_file : str
            Path to minimed structure
        """
        input_file = os.path.join(self.work_dir, "{}_{}".format(prefix, os.path.basename(template_file)))

        with open(template_file) as temp, open(input_file, "w") as inp:
            for line in temp:
                for key, val in list(kwds.items()):
                    line = re.sub(key, val, line)
                inp.write(line)

        return input_file

    def set_local_env(self):
        if "CNS_TOPPAR" not in os.environ:
            pass

        try:
            cns = subprocess.check_output(["which", "cns_solve"])
        except:
            return None
            #raise RuntimeError("Unable to load CNS")
        cns_dir = os.path.dirname(os.path.dirname(os.path.dirname(cns)))
        env_file = os.path.join(cns_dir, "cns_solve_env.zsh")
        output = subprocess.check_output("source {}; env -0".format(env_file),
            shell=True, executable="/bin/zsh")
        new_env = dict([line.partition('=')[::2] for line in output.split('\0')])
        os.environ.clear()
        os.environ.update(new_env)


class CNSMinimize(CNS):
    LOCAL = ["cns"]
    PARAMETERS = [("pdb_file", "path:in", ""), ("output_file", "path:out")]
    RETURN_FILES = True
    cns4char_re = re.compile("%CREAD-ERR: residue ID and insertion character ([0-9a-zA-Z]+) exceed 4 characters.")
    cns_output_re = re.compile("\| Etotal =([0-9\.]+)\s+grad\(E\)=([0-9\.]+)\s+E\(BOND\)=([0-9\.]+)\s+E\(ANGL\)=([0-9\.]+)\s+\|\n \| E\(DIHE\)=([0-9\.]+)\s+E\(IMPR\)=([0-9\.]+)\s+E\(VDW \)=([0-9\.]+)\s+E\(ELEC\)=([0-9\.]+)\s+\|\n -{79}\n LBFGS: normal termination - NSTEP limit reached")


    def __call__(self, pdb_file, output_file_prefix=None):
        return self.run(pdb_file=pdb_file, output_file_prefix=output_file_prefix)

    def run(self, pdb_file, output_file_prefix=None):
        output_file_prefix = output_file_prefix or os.path.splitext(pdb_file)[0]
        minimized_file = "{}.min".format(output_file_prefix)

        _pdb_file = self.format_in_path(None, pdb_file)
        output_file_prefix = self.format_in_path(None, output_file_prefix)

        temp_minimization_file = os.path.join(script_dir, "CNS-Templates", "model_minimize.inp")
        input_file = self._generate_input(temp_minimization_file, output_file_prefix,
            INSERT_PDB_HERE=_pdb_file, INSERT_OUTPUT_FILE=minimized_file)

        #Call CNS
        self.stdout = super()(input_file=input_file)

        minimized_file = self.check_output()
        minimized_file = self.fix_output(pdb_file, minimized_file)

        energy_results = Energy.parse_energy(self.stdout)

        return out_files

    def fix_output(self, pdb_file, minimized_file):
        m = cns4char_re.search(self.stdout)
        if m:
            if len(m.group(1)) < 6:
                #Strange bug in CNS where resn is 4 and icode is 1, valid PDB but CNS chokes
                #CNS corrects them by removing with first character
                #Fix by copying the original PDB cols before the XYZ coordiantes
                new_minimized_file = os.path.splitext(minimized_file)[0]+".fixed.min" #os.path.join(work_dir, "{}.fixed.min".format(output_file_prefix))
                with open(pdb_file) as old, open(minimized_file) as minf, open(new_minimized_file, "w") as new:
                     for oline, mline in izip_missing(old, minf, fillvalue=None):
                        if oline is None or mline is None:
                            raise RuntimeError("CNS minimization has failed!! {}. Check output: {}".format(m.group(0), cns))
                        new.write(mline[0:22]+oline[22:27]+mline[27:])
                minimized_file = new_minimized_file
            else:
                raise RuntimeError("CNS minimization has failed!! Check output: {}".format(cns))

        cmds = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), minimized_file]]
        tidy_file = minimized_file+".tidy"
        with open(tidy_file, "w") as f:
            SubprocessChain(cmds, f)
        minimized_file = tidy_file

        return minimized_file

        m = cns_output_re.search(cns)
        if m:
            results = pd.Series({
                 "cns_Etotal": m.group(1),
                 "cns_grad(E)": m.group(2),
                 "cns_E(BOND)": m.group(3),
                 "cns_E(ANGL)": m.group(4),
                 "cns_E(DIHE)": m.group(5),
                 "cns_E(IMPR)": m.group(6),
                 "cns_E(VDW )": m.group(7),
                 "cns_E(ELEC)": m.group(8)})
        else:
            raise RuntimeError("CNS minimization has failed!! Check output: {}".format(cns))

        return minimized_file, results

class Energy(CNS):
    LOCAL = ["cns"]
    PARAMETERS = [("pdb_file", "path:in")]
    RETURN_FILES = True
    cns_energy_re = re.compile("\| Etotal =([0-9\.]+)\s+grad\(E\)=([0-9\.]+)\s+E\(BOND\)=([0-9\.]+)\s+E\(ANGL\)=([0-9\.]+)\s+\|\n \| E\(DIHE\)=([0-9\.]+)\s+E\(IMPR\)=([0-9\.]+)\s+E\(VDW \)=([0-9\.]+)\s+E\(ELEC\)=([0-9\.]+)\s+\|")

    def __call__(self, pdb_file):
        return self.run(pdb_file=pdb_file)

    def run(self, pdb_file):
        _pdb_file, output_file_prefix = self.format_parameters(pdb_file=pdb_file)

        temp_energy_file = os.path.join(script_dir, "CNS-Templates", "model_energy.inp")

        input_file = self._generate_input(temp_energy_file, pdb_file+".out",
            INSERT_PDB_HERE=_pdb_file)

        #Call CNS
        self.stdout = super(self, Energy)(input_file=input_file)

        energy_results = self.parse_energy(self.stdout)

        return out_files

    @staticmethod
    def parse_energy(self, stdout):
        for m in self.cns_energy_re.finditer(cns):
            #Get last match
            pass

        if m:
            job.log("GROUPS: {}".format(m.groups()))
            results = pd.Series({
                 "cns_Etotal": m.group(1),
                 "cns_grad(E)": m.group(2),
                 "cns_E(BOND)": m.group(3),
                 "cns_E(ANGL)": m.group(4),
                 "cns_E(DIHE)": m.group(5),
                 "cns_E(IMPR)": m.group(6),
                 "cns_E(VDW )": m.group(7),
                 "cns_E(ELEC)": m.group(8)})
        else:
            raise RuntimeError("CNS minimization has failed!! Check output: {}".format(cns))

        return results
