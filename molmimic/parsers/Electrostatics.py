import os, sys
import subprocess
import shutil
from itertools import cycle

from joblib import Memory
from docker.errors import ImageNotFound, ContainerError

from molmimic.generate_data.util import silence_stdout, silence_stderr, \
    get_all_chains, extract_chains, SubprocessChain, PDB_TOOLS
from molmimic.generate_data.util import remove_ter_lines as _remove_ter_lines
from molmimic.parsers.psize import Psize

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    pdb2pqr_src = os.path.join(os.path.dirname(subprocess.check_output(["which", "pdb2pqr"])), "src")
    if pdb2pqr_src:
        sys.path.append(pdb2pqr_src)

memory = Memory(verbose=0)

def make_apbs_input(pqr_file, output_prefix=None):
    ps = Psize()
    ps.runPsize(pqr_file)
    cglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getCoarseGridDims())
    fglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getFineGridDims())
    dime = "{:d} {:d} {:d}".format(*ps.getFineGridPoints())

    output_prefix = output_prefix or pqr_file+"_apbs_output"

    return """read
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
    pqr=os.path.basename(pqr_file),
    dime=dime,
    cglen=cglen,
    fglen=fglen,
    out=output_prefix
    )

def run_apbs(pqr_file, input_file=None, keep_input=False, work_dir=None,
  docker=True, job=None, attempts=3):
    """Run APBS. Calculates correct size using Psize and defualt from Chimera
    """
    if work_dir is None:
        work_dir = os.getcwd()

    full_pqr_path = pqr_file
    pqr_path = os.path.basename(full_pqr_path)
    file_prefix = os.path.splitext(pqr_path)[0]
    output_prefix = os.path.join(work_dir, "{}.apbs_output".format(file_prefix))
    if input_file is not None and os.path.isfile(input_file):
        keep_input = True
    else:
        input_file_contents = make_apbs_input(full_pqr_path, "{}.apbs_output".format(file_prefix))


    if docker and apiDockerCall is not None and job is not None:
        input_file_name = os.path.join(work_dir, "{}.apbs_input".format(file_prefix))
        input_file_short = "{}.apbs_input".format(file_prefix)
        output_prefix = "{}.apbs_output".format(file_prefix)

        if input_file is not None:
            if not os.path.abspath(os.path.dirname(input_file)) == os.path.abspath(work_dir):
                shutil.copy(input_file, os.path.join(work_dir, input_file_short))
            else:
                input_file_short = os.path.basename(input_file)
        else:
            input_file = input_file_name
            with open(input_file, "w") as f:
                f.write(input_file_contents)

        if not os.path.abspath(os.path.dirname(pqr_file)) == os.path.abspath(work_dir):
            shutil.copy(pqr_file, work_dir)

        try:
            parameters = ["/data/{}".format(input_file_short)]
            apiDockerCall(job,
                          image='edraizen/apbs:latest',
                          working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
                          parameters=parameters)
            output_prefix = os.path.join(work_dir, output_prefix)
        except (SystemExit, KeyboardInterrupt):
            raise
        except ImageNotFound:
            if attempts > 0:
                return run_apbs(pqr_file, input_file=input_file, keep_input=keep_input,
                    work_dir=work_dir, docker=docker, job=job, attempts=attempts-1)
            raise
        # except:
        #     raise
        #     return run_apbs(full_pqr_file, input_file=input_file_name,
        #         keep_input=keep_input, work_dir=work_dir, docker=False)

    else:
        input_file = os.path.join(work_dir, "{}.apbs_input".format(file_prefix))
        output_prefix = os.path.join(work_dir, "{}.apbs_output".format(file_prefix))
        with open(input_file, "w") as f:
            f.write(input_file_contents)

        try:
            subprocess.call(["apbs", input_file])
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    if not keep_input and input_file is not None and os.path.isfile(input_file):
        os.remove(input_file)

    out_file = output_prefix+".txt"
    assert os.path.isfile(out_file), "Outfile not found: {}".format(os.listdir(work_dir))
    return out_file

def run_pdb2pqr(pdb_file, whitespace=True, ff="amber", parameters=None,
  remove_ter_lines=True, tidy=False, chain=True, work_dir=None, docker=True,
  job=None, attempts=3):
    if work_dir is None:
        work_dir = os.getcwd()

    save_chains = get_all_chains(pdb_file)

    full_pdb_path = _remove_ter_lines(pdb_file) if remove_ter_lines else pdb_file
    pdb_path = os.path.basename(full_pdb_path)
    pqr_file = "{}.pqr".format(pdb_path)

    _parameters = list(parameters) if isinstance(parameters, (list, tuple)) else []
    _parameters.append("--ff={}".format(ff))
    if whitespace:
        _parameters.append("--whitespace")
    if chain and "--chain" not in _parameters:
        _parameters.append("--chain")

    if docker and apiDockerCall is not None and job is not None:
        #Docker can only read from work_dir
        if not os.path.abspath(os.path.dirname(full_pdb_path)) == os.path.abspath(work_dir):
            shutil.copy(full_pdb_path, work_dir)

        _parameters += ["/data/{}".format(pdb_path), "/data/{}.pqr".format(pdb_path)]
        try:
            output = apiDockerCall(job,
                          image='edraizen/pdb2pqr:latest',
                          working_dir=work_dir,
                          parameters=_parameters)
            pqr_file = os.path.join(work_dir, pqr_file)
        except (SystemExit, KeyboardInterrupt):
            raise
        except ValueError as e:
            if str(e).startswith("This PDB file is missing too many"):
                return None
            raise
        except ContainerError as e:
            if "ValueError: This PDB file is missing too many" in str(e):
                return None
            raise
        except ImageNotFound:
            if attempts > 0:
                return run_pdb2pqr(pdb_file, whitespace=whitespace, ff=ff,
                    parameters=parameters, remove_ter_lines=remove_ter_lines,
                    tidy=tidy, chain=chain, work_dir=work_dir, docker=docker,
                    job=job, attempts=attempts-1)
            raise
        #except:
            #raise
            #return run_pdb2pqr(pdb_file, whitespace=whitespace, ff=ff,
            #    parameters=parameters, work_dir=work_dir, docker=False)

    else:
        pqr_file = os.path.join(work_dir, full_pdb_path)
        command = ["/usr/share/pdb2pqr/pdb2pqr.py"]+parameters
        command += [full_pdb_path, pqr_file]

        try:
            with silence_stdout(), silence_stderr():
                subprocess.call(command)
        except (SystemExit, KeyboardInterrupt):
            raise
        except ValueError as e:
            if str(e).startswith("This PDB file is missing too many"):
                return None
            raise
        #except Exception as e:
            #raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    if tidy:
        pqr_chains = extract_chains(pqr_file, save_chains)
        cmds = [[sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), pqr_chains]]
        tidy_file = pqr_file+".tidy"
        with open(tidy_file, "w") as f:
            SubprocessChain(cmds, f)
        try:
            os.remove(pqr_file)
        except OSError:
            pass
        shutil.move(tidy_file, pqr_file)

    assert os.path.isfile(pqr_file)
    return pqr_file

@memory.cache
def run_pdb2pqr_APBS(pdb_path, pdb2pqr_whitespace=False, pdb2pqr_ff="amber",
  pdb2pqr_parameters=None, apbs_input_file=None, apbs_keep_input=False,
  work_dir=None, docker=True, job=None, clean=True, only_charge=False):
    """Run PDB2PQR and APBS to get charge and electrostatics for each atom

    Parameters
    ----------
    pdb_path : str
    Path to PDB file
    job : Toil.job.Job
    The job that is calling this method if using toil
    clean : bool
    Remove the resulting PQR and APBS files. Default is True.

    Return
    ------
    A dictionary mapping atoms bu Bio.PDB identifiers to a tuple of floats:
    charge, and
    electrostatic potential
    """
    pqr_path = run_pdb2pqr(pdb_path, whitespace=pdb2pqr_whitespace,
        ff=pdb2pqr_ff, parameters=pdb2pqr_parameters, work_dir=work_dir,
        docker=docker, job=job)

    if pqr_path is None:
        return {}

    assert os.path.isfile(pqr_path)

    if not only_charge:
        atom_pot_file = run_apbs(pqr_path, input_file=apbs_input_file,
            keep_input=apbs_keep_input, work_dir=work_dir, docker=docker, job=job)
        atom_pot = open(atom_pot_file)
    else:
        #Set all atom potentials to NaN
        atom_pot = cycle([np.nan])

    result = {}
    with open(pqr_path) as pqr:
        #Skip first 4 rows of atompot file
        for _ in range(4):
            next(atom_pot)

        for line in pqr:
            if not line.startswith("ATOM  ") or line.startswith("HETATM"): continue

            electrostatic_potential = float(next(atom_pot))

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
            result[key] = (float(charge), electrostatic_potential)

    if not only_charge:
        atom_pot.close()
        if clean:
            try:
                os.remove(atom_pot_file)
            except OSError:
                pass
    del atom_pot

    if clean:
        try:
            os.remove(pqr_path)
        except OSError:
            pass

    return result
