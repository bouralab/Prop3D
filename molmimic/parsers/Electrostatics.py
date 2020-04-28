import os, sys
import subprocess
import shutil
from itertools import cycle

from joblib import Memory
from docker.errors import ImageNotFound, ContainerError

from molmimic.util import silence_stdout, silence_stderr, SubprocessChain, safe_remove
from molmimic.util.pdb import get_all_chains, extract_chains, get_atom_lines, \
    replace_occ_b, PDB_TOOLS, split_xyz
from molmimic.util.pdb import remove_ter_lines as _remove_ter_lines
from molmimic.util.pdb import get_first_chain, replace_chains
from molmimic.parsers.psize import Psize

from molmimic.parsers.container import Container

from toil.realtimeLogger import RealtimeLogger

memory = Memory(verbose=0)

class APBS(Container):
    IMAGE = 'docker://edraizen/apbs:latest'
    LOCAL = ["apbs"]
    PARAMETERS = [("in_file", "path:in")]#, (":out_file", "path:out")]
    RETURN_FILES = True

    def atom_potentials_from_pdb(self, pdb_file, force_field="amber", with_charge=True):
        remove_pdb = False
        if get_first_chain(pdb_file).isdigit():
            #APBS requires chains A-Z
            pdb_file = replace_chains(pdb_file, new_chain="A")
            remove_pdb = True

        pdb2pqr = Pdb2pqr(work_dir=self.work_dir, job=self.job)
        pqr_file = pdb2pqr.create_pqr(pdb_file, force_field=force_field,
            whitespace=True, chain=True, reres=True)
        atom_pot_file = self.atom_potentials_from_pqr(pqr_file)

        if with_charge:
            output = {atom:(charge, potential) for (atom, charge), potential in \
                zip(pdb2pqr.parse_pqr(pqr_file),
                    self.parse_atom_pot_file(atom_pot_file))}
        else:
            output = {atom:potential for (atom, charge), potential in \
                zip(pdb2pqr.parse_pqr(pqr_file),
                    self.parse_atom_pot_file(atom_pot_file))}

        self.files_to_remove += [pqr_file, atom_pot_file]
        if remove_pdb:
            self.files_to_remove.append(pdb_file)
        self.clean()

        return output

    def atom_potentials_from_pqr(self, pqr_file):
        self.full_pqr_path = pqr_file

        pqr_path = self.format_in_path(None, pqr_file)

        input_file, out_file = self.write_apbs_electrostatics_input(pqr_file)
        self.files_to_remove.append(input_file)

        try:
            self.running_atom_potentials = True
            atom_pot_file = self(in_file=input_file)
            self.running_atom_potentials = False
        except RuntimeError as e:
            if "Problem opening virtual socket" in str(e):
                try:
                    with open(pqr_file) as f:
                        raise RuntimeError("{} \n PQR file: {}".format(str(e), f.read()))
                except:
                    pass
            raise e

        return atom_pot_file

    def local(self, *args, **kwds):
        self.is_local = True
        if hasattr(self, "running_atom_potentials") and self.running_atom_potentials:
            pqr_path = self.format_in_path(None, self.full_pqr_path)
            input_file, out_file = self.write_apbs_electrostatics_input(self.full_pqr_path)
            kwds["in_file"] = input_file
        super(self, APBS).local(*args, **kwds)

    @staticmethod
    def parse_atom_pot_file(atom_pot_file):
        with open(atom_pot_file) as atom_pots:
            #Skip first 4 rows of atompot file
            for _ in range(4):
                next(atom_pots)

            for atom_pot in atom_pots:
                yield float(atom_pot)

    def write_apbs_electrostatics_input(self, pqr_file, ouput_input_prefix=None,
      output_prefix=None):
        base = os.path.basename(pqr_file)
        ouput_input_prefix = os.path.join(self.work_dir, ouput_input_prefix or base+".apbs_input")
        output_prefix = os.path.join(self.work_dir, output_prefix or base+".apbs_output.txt")
        pot_path = self.format_out_path(None, output_prefix)
        with open(ouput_input_prefix, "w") as f:
            contents = self.make_apbs_electrostatics_input(pqr_file, pot_path[:-4])
            f.write(contents)
        return ouput_input_prefix, pot_path

    def make_apbs_electrostatics_input(self, pqr_file, output_prefix):
        ps = Psize()
        ps.runPsize(pqr_file)
        pqr_path = self.format_in_path(None, pqr_file)

        cglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getCoarseGridDims())
        fglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getFineGridDims())
        dime = "{:d} {:d} {:d}".format(*ps.getFineGridPoints())

        contents = """read
        mol pqr {pqr}
    end
    elec
        mg-auto # Specify the mode for APBS to run
        dime {dime} # The grid dimensions
        cglen {cglen}
        fglen {fglen}
        cgcent mol 1
        fgcent mol 1
        mol 1 # Perform the calculation on molecule 1
        lpbe # Solve the linearized Poisson-Boltzmann equation
        bcfl sdh # Use all single moments when calculating the potential
        pdie 2.00 # Solute dielectric
        sdie 2.00 #78.54 # Solvent dielectric
        chgm spl2 # Spline-based discretization of the delta functions
        srfm smol # Molecular surface definition
        srad 1.4 # Solvent probe radius (for molecular surface)
        swin 0.3 # Solvent surface spline window (not used here)
        sdens 10.0 # Sphere density for accessibility object
        temp 298.15 # Temperature
        calcenergy total # Calculate energies
        calcforce no # Do not calculate forces
        write atompot flat {out}
    end
    quit
    """.format(
        pqr=pqr_path,
        dime=dime,
        cglen=cglen,
        fglen=fglen,
        out=output_prefix
        )

        return contents

class MissingAtomsError(ValueError):
    pass

class Pdb2pqr(Container):
    IMAGE = 'docker://edraizen/pdb2pqr:latest'
    LOCAL = ["/usr/share/pdb2pqr/pdb2pqr.py"]
    PARAMETERS = [
        (":force_field:amber", None, "--ff={}"),
        (":whitespace", "store_true"),
        (":chain", "store_true"),
        ("in_file", "path:in"),
        ("out_file", "path:out")]
    RETURN_FILES = True

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

# def run_pdb2pqr(pdb_file, whitespace=True, ff="amber", parameters=None,
#   remove_ter_lines=True, tidy=False, chain=True, work_dir=None, docker=True,
#   job=None, attempts=3):
#     if work_dir is None:
#         work_dir = os.getcwd()
#
#     save_chains = get_all_chains(pdb_file)
#
#     full_pdb_path = _remove_ter_lines(pdb_file) if remove_ter_lines else pdb_file
#     pdb_path = os.path.basename(full_pdb_path)
#     pqr_file = "{}.pqr".format(pdb_path)

    # _parameters = list(parameters) if isinstance(parameters, (list, tuple)) else []
    # _parameters.append("--ff={}".format(ff))
    # if whitespace:
    #     _parameters.append("--whitespace")
    # if chain and "--chain" not in _parameters:
    #     _parameters.append("--chain")
    #
    # if docker and apiDockerCall is not None and job is not None:
    #     #Docker can only read from work_dir
    #     if not os.path.abspath(os.path.dirname(full_pdb_path)) == os.path.abspath(work_dir):
    #         shutil.copy(full_pdb_path, work_dir)
    #
    #     _parameters += ["/data/{}".format(pdb_path), "/data/{}.pqr".format(pdb_path)]
    #     try:
    #         output = apiDockerCall(job,
    #                       image='edraizen/pdb2pqr:latest',
    #                       working_dir=work_dir,
    #                       parameters=_parameters)
    #         pqr_file = os.path.join(work_dir, pqr_file)
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except ValueError as e:
    #         if str(e).startswith("This PDB file is missing too many"):
    #             return None
    #         raise
        # except ContainerError as e:
        #     if "ValueError: This PDB file is missing too many" in str(e):
        #         return None
        #     raise
        # except ImageNotFound:
        #     if attempts > 0:
        #         return run_pdb2pqr(pdb_file, whitespace=whitespace, ff=ff,
        #             parameters=parameters, remove_ter_lines=remove_ter_lines,
        #             tidy=tidy, chain=chain, work_dir=work_dir, docker=docker,
        #             job=job, attempts=attempts-1)
        #     raise
        #except:
            #raise
            #return run_pdb2pqr(pdb_file, whitespace=whitespace, ff=ff,
            #    parameters=parameters, work_dir=work_dir, docker=False)

    # else:
    #     pqr_file = os.path.join(work_dir, full_pdb_path)
    #     command = ["/usr/share/pdb2pqr/pdb2pqr.py"]+parameters
    #     command += [full_pdb_path, pqr_file]
    #
    #     try:
    #         with silence_stdout(), silence_stderr():
    #             subprocess.call(command)
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except ValueError as e:
    #         if str(e).startswith("This PDB file is missing too many"):
    #             return None
    #         raise
        #except Exception as e:
            #raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    # if tidy:
    #     pqr_chains = extract_chains(pqr_file, save_chains)
    #     cmds = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), pqr_chains]]
    #     tidy_file = pqr_file+".tidy"
    #     with open(tidy_file, "w") as f:
    #         SubprocessChain(cmds, f)
    #     try:
    #         os.remove(pqr_file)
    #     except OSError:
    #         pass
    #     shutil.move(tidy_file, pqr_file)
    #
    # assert os.path.isfile(pqr_file)
    # return pqr_file

#@memory.cache
# def run_pdb2pqr_APBS(pdb_path, pdb2pqr_whitespace=False, pdb2pqr_ff="amber",
#   pdb2pqr_parameters=None, apbs_input_file=None, apbs_keep_input=False,
#   work_dir=None, docker=True, job=None, clean=True, only_charge=False):
#     """Run PDB2PQR and APBS to get charge and electrostatics for each atom
#
#     Parameters
#     ----------
#     pdb_path : str
#     Path to PDB file
#     job : Toil.job.Job
#     The job that is calling this method if using toil
#     clean : bool
#     Remove the resulting PQR and APBS files. Default is True.
#
#     Return
#     ------
#     A dictionary mapping atoms bu Bio.PDB identifiers to a tuple of floats:
#     charge, and
#     electrostatic potential
#     """
#     pqr_path = run_pdb2pqr(pdb_path, whitespace=pdb2pqr_whitespace,
#         ff=pdb2pqr_ff, parameters=pdb2pqr_parameters, work_dir=work_dir,
#         docker=docker, job=job)
#
#     if pqr_path is None:
#         return {}
#
#     assert os.path.isfile(pqr_path)
#
#     if not only_charge:
#         atom_pot_file = run_apbs(pqr_path, input_file=apbs_input_file,
#             keep_input=apbs_keep_input, work_dir=work_dir, docker=docker, job=job)
#         atom_pot = open(atom_pot_file)
#     else:
#         #Set all atom potentials to NaN
#         atom_pot = cycle([np.nan])
#
#     result = {}
#     with open(pqr_path) as pqr:
#         #Skip first 4 rows of atompot file
#         for _ in range(4):
#             next(atom_pot)
#
#         for line in pqr:
#             if not line.startswith("ATOM  ") or line.startswith("HETATM"): continue
#
#             electrostatic_potential = float(next(atom_pot))
#
#             fields = line.rstrip().split()
#             if len(fields) == 11:
#                 recordName, serial, atomName, residueName, chainID, residueNumber, X, Y, Z, charge, radius = fields
#             elif len(fields) == 10:
#                 recordName, serial, atomName, residueName, residueNumber, X, Y, Z, charge, radius = fields
#             else:
#                 raise RuntimeError("Invalid PQR File")
#
#             try:
#                 resseq = int("".join([i for i in residueNumber if i.isdigit()]))
#             except ValueError:
#                 continue
#
#             icode = "".join([i for i in residueNumber if not i.isdigit()])
#             if icode == "":
#                 icode = " "
#
#             if recordName == "HETATM":  # hetero atom flag
#                 if residueName in ["HOH", "WAT"]:
#                     hetero_flag = "W"
#                 else:
#                     hetero_flag = "H_{}".format(residueName)
#             else:
#                 hetero_flag = " "
#
#             residue_id = (hetero_flag, resseq, icode)
#
#             key = (residue_id, (atomName.strip(), ' '))
#             result[key] = (float(charge), electrostatic_potential)
#
#     if not only_charge:
#         atom_pot.close()
#         if clean:
#             try:
#                 os.remove(atom_pot_file)
#             except OSError:
#                 pass
#     del atom_pot
#
#     if clean:
#         try:
#             os.remove(pqr_path)
#         except OSError:
#             pass
#
#     return result

# def run_apbs(pqr_file, input_file=None, keep_input=False, work_dir=None,
#   docker=True, job=None, attempts=3):
#     """Run APBS. Calculates correct size using Psize and defualt from Chimera
#     """
#     if work_dir is None:
#         work_dir = os.getcwd()
#
#     full_pqr_path = pqr_file
#     pqr_path = os.path.basename(full_pqr_path)
#     file_prefix = os.path.splitext(pqr_path)[0]
#     output_prefix = os.path.join(work_dir, "{}.apbs_output".format(file_prefix))
#     if input_file is not None and os.path.isfile(input_file):
#         keep_input = True
#     else:
#         input_file_contents = make_apbs_input(full_pqr_path, "{}.apbs_output".format(file_prefix))
#
#
#     if docker and apiDockerCall is not None and job is not None:
#         input_file_name = os.path.join(work_dir, "{}.apbs_input".format(file_prefix))
#         input_file_short = "{}.apbs_input".format(file_prefix)
#         output_prefix = "{}.apbs_output".format(file_prefix)
#
#         if input_file is not None:
#             if not os.path.abspath(os.path.dirname(input_file)) == os.path.abspath(work_dir):
#                 shutil.copy(input_file, os.path.join(work_dir, input_file_short))
#             else:
#                 input_file_short = os.path.basename(input_file)
#         else:
#             input_file = input_file_name
#             with open(input_file, "w") as f:
#                 f.write(input_file_contents)

    #     if not os.path.abspath(os.path.dirname(pqr_file)) == os.path.abspath(work_dir):
    #         shutil.copy(pqr_file, work_dir)
    #
    #     try:
    #         parameters = ["/data/{}".format(input_file_short)]
    #         apiDockerCall(job,
    #                       image='edraizen/apbs:latest',
    #                       working_dir="/data",
    #                       volumes={work_dir:{"bind":"/data", "mode":"rw"}},
    #                       parameters=parameters)
    #         output_prefix = os.path.join(work_dir, output_prefix)
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except ImageNotFound:
    #         if attempts > 0:
    #             return run_apbs(pqr_file, input_file=input_file, keep_input=keep_input,
    #                 work_dir=work_dir, docker=docker, job=job, attempts=attempts-1)
    #         raise
    #     # except:
    #     #     raise
    #     #     return run_apbs(full_pqr_file, input_file=input_file_name,
    #     #         keep_input=keep_input, work_dir=work_dir, docker=False)
    #
    # else:
    #     input_file = os.path.join(work_dir, "{}.apbs_input".format(file_prefix))
    #     output_prefix = os.path.join(work_dir, "{}.apbs_output".format(file_prefix))
    #     with open(input_file, "w") as f:
    #         f.write(input_file_contents)
    #
    #     try:
    #         subprocess.call(["apbs", input_file])
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except Exception as e:
    #         raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))
    #
    # if not keep_input and input_file is not None and os.path.isfile(input_file):
    #     os.remove(input_file)
    #
    # out_file = output_prefix+".txt"
    # assert os.path.isfile(out_file), "Outfile not found: {}".format(os.listdir(work_dir))
    # return out_file
