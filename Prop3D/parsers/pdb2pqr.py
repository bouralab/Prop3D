import os
import sys
import shutil

from Prop3D.util import safe_remove, SubprocessChain
from Prop3D.parsers.container import Container
from Prop3D.util.pdb import remove_ter_lines as _remove_ter_lines
from Prop3D.util.pdb import get_all_chains, extract_chains, \
    replace_occ_b, PDB_TOOLS

class MissingAtomsError(ValueError):
    pass

class Pdb2pqr(Container):
    IMAGE = 'docker://edraizen/pdb2pqr:latest'
    LOCAL = ["/usr/share/pdb2pqr/pdb2pqr.py"]
    PARAMETERS = [
        (":force_field:amber", None, "ff"),
        (":whitespace", "store_true"),
        (":chain", "store_true"),
        (":noopt", "store_true"),
        (":apbs_input", "store_true", ["--apbs-input"]),
        ("in_file", "path:in", ["{}"]),
        ("out_file", "path:out", ["{}"])]
    RETURN_FILES = True
    ARG_START = "--"
    ARG_SEP = "="

    def __call__(self, *args, **kwds):
        try:
            return Container.__call__(self, *args, **kwds)
        except Exception as e:
            self.clean()
            if "ValueError: This PDB file is missing too many" in str(e):
                raise MissingAtomsError(str(e)[13:])
            raise

    def create_pqr(self, pdb_file, remove_ter_lines=True, whitespace=False,
      chain=False, **kwds):
        """Run pdb2pqr for a given pdb_file.

        Parameters
        ----------
        pdb_file : str
            Path to pdb file
        remove_ter_lines : bool
            Remove TER lines before running since pdb2pqr may choke on them.
            Default True.
        **kwds:
            Paramters to pass to pdb2pqr.

        Returns
        -------
        A path the the new pqr file
        """
        pqr_file = "{}.pqr".format(pdb_file)

        if remove_ter_lines:
            pdb_file = _remove_ter_lines(pdb_file)

        #Whitespace must be False to ensure correct PDB file
        if whitespace:
            kwds["whitespace"] = True

        if chain:
            kwds["chain"] = True

        try:
            output = self(pdb_file, pqr_file, **kwds)
        except MissingAtomsError:
            if remove_ter_lines:
                safe_remove(pdb_file)
            raise

        if remove_ter_lines:
            safe_remove(pdb_file)

        return pqr_file

    def debump_add_h(self, pdb_file, remove_ter_lines=True, keep_occ_and_b=True, **kwds):
        """Run pdb2pqr for a given pdb_file.

        Parameters
        ----------
        pdb_file : str
            Path to pdb file
        remove_ter_lines : bool
            Remove TER lines before running since pdb2pqr may choke on them.
            Default True.
        keep_occ_and_b : None, bool, or 2-tuple
            Replace the new partial charge and radius fields with original
            occupancy and bfactor values and create new occupancy and bfactor
            values for added hydrogens in order to create a valid PDB file with
            standard column sizes.
                None, False: keep the partial charge and radius fields
                True: Replace fields with originals and new hydogens get an
                    occupancy of 1.0 and bfactor of 20.0
                2-tuple (occupancy, bfactor): Replace fields with originals and
                    new hydogens get the first tuple value, bfactors get the
                    second tuple value.
        **kwds:
            Paramters to pass to pdb2pqr.

        Returns
        -------
        A path the the new pqr file
        """
        save_chains = "".join(list(sorted(get_all_chains(pdb_file))))

        pqr_file = self.create_pqr(pdb_file, remove_ter_lines=remove_ter_lines,
            keep_occ_and_b=keep_occ_and_b, whitespace=False, chain=True, **kwds)

        new_pdb = extract_chains(pqr_file, save_chains)

        cmds = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), new_pdb]]
        tidy_file = new_pdb+".tidy"
        with open(tidy_file, "w") as f:
            SubprocessChain(cmds, f)

        safe_remove([new_pdb, pqr_file])
        shutil.move(tidy_file, new_pdb)

        if keep_occ_and_b:
            occ_b_file = replace_occ_b(pdb_file, new_pdb)
            safe_remove(new_pdb)
            shutil.move(occ_b_file, new_pdb)

        return new_pdb

    def get_charge_from_pdb_file(self, pdb_file, remove_ter_lines=True, **kwds):
        pqr_file = self.create_pqr(pdb_file, remove_ter_lines=remove_ter_lines, **kwds)
        return dict(self.parse_pqr(pqr_file))

    @staticmethod
    def parse_pqr(pqr_file):
        with open(pqr_file) as pqr:
            for line in pqr:
                if not line.startswith("ATOM  ") or line.startswith("HETATM"): continue

                fields = line.rstrip().split()
                if len(fields) == 11:
                    recordName, serial, atomName, residueName, chainID, residueNumber, X, Y, Z, charge, radius = fields
                elif len(fields) == 10:
                    recordName, serial, atomName, residueName, residueNumber, X, Y, Z, charge, radius = fields
                else:
                    raise RuntimeError("Invalid PQR File")

                try:
                    resseq = int("".join([i for i in residueNumber if i.isdigit()]))
                except ValueError:
                    continue

                icode = "".join([i for i in residueNumber if not i.isdigit()])
                if icode == "":
                    icode = " "

                if recordName == "HETATM":  # hetero atom flag
                    if residueName in ["HOH", "WAT"]:
                        hetero_flag = "W"
                    else:
                        hetero_flag = "H_{}".format(residueName)
                else:
                    hetero_flag = " "

                residue_id = (hetero_flag, resseq, icode)

                key = (residue_id, (atomName.strip(), ' '))

                yield key, float(charge)
                #result[key] = (float(charge), electrostatic_potential)
