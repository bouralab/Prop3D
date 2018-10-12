import os, sys
import subprocess

from joblib import Memory

from molmimic.util import silence_stdout, silence_stderr
from molmimic.parsers.psize import Psize

try:
    from toil.lib.docker import apiDockerCall
except ImportError:
    apiDockerCall = None
    pdb2pqr_src = os.path.join(os.path.dirname(subprocess.check_output(["which", "pdb2pqr"])), "src")
    if pdb2pqr_src:
        sys.path.append(pdb2pqr_src)

memory = Memory(verbose=0)

def make_apbs_input(pqr_file):
    ps = Psize()
    ps.runPsize(pqr_file)
    cglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getCoarseGridDims())
    fglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getFineGridDims())
    dime = "{:d} {:d} {:d}".format(*ps.getFineGridPoints())

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
    pqr=pqr_file,
    dime=dime,
    cglen=cglen,
    fglen=fglen,
    out=output_prefix
    )

@memory.cache
def run_apbs(pqr_file, input_file=None, keep_input=False, work_dir=None, docker=True, job=None):
    """Run APBS. Calculates correct size using Psize and defualt from Chimera
    """
    if work_dir is None:
        work_dir = os.getcwd()

    full_pqr_path = pqr_file
    pqr_path = os.path.basename(full_pqr_path)
    file_prefix = os.path.splitext(pqr_path)[0]
    if input_file is not None and os.path.isfile(input_file):
        keep_input = True
    else:
        input_file_contents = make_apbs_input(full_pqb_path)


    if docker and apiDockerCall is not None and job is not None:
        input_file_name = os.path.join(work_dir, "{}.apbs_input".format(file_prefix))
        input_file_short = "{}.apbs_input".format(file_prefix)
        output_prefix = "{}.apbs_output".format(file_prefix)

        if input_file is not None and not os.path.abspath(os.path.dirname(input_file)) == os.path.abspath(work_dir):
            shutil.copy(input_file, work_dir)
        else:
    		with open(input_file_name, "w") as f:
    			f.write(input_file_contents)

        try:
            parameters = ["/data/{}".format(input_file_short)]
            apiDockerCall(job,
                          image='edraizen/apbs:latest',
                          working_dir=work_dir,
                          parameters=parameters)
            output_prefix = os.path.join(work_dir, output_prefix)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            return run_apbs(full_pqr_file, input_file=input_file_name,
                keep_input=keep_input, work_dir=work_dir, docker=False)

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

    if not keep_input:
        os.remove(input_file)

    out_file = output_prefix+".txt"
    assert os.path.isfile(out_file)
    return out_file

def run_pdb2pqr(pdb_file, whitespace=True, ff="amber", parameters=None, work_dir=None, docker=True, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    full_pdb_path = pdb_file
    pdb_path = os.path.basename(full_pdb_path)
    pqr_file = "{}.pqr".format(pdb_path)

    _parameters = list(parameters) if isinstance(parameters, (list, tuple)) else []
    _parameters.append("--ff={}".format(ff))
    if whitespace:
        _parameters.add("--whitespace")

    if docker and apiDockerCall is not None and job is not None:
        #Docker can only read from work_dir
        if not os.path.abspath(os.path.dirname(pdb_file)) == os.path.abspath(work_dir):
            shutil.copy(pdb_file, work_dir)

        _parameters += ["/data/{}".format(pdb_path), "/data/{}.pqr".format(pdb_path)]
        try:
            apiDockerCall(job,
                          image='edraizen/pdb2pqr:latest',
                          working_dir=work_dir,
                          parameters=_parameters)
            pqr_file = os.path.join(work_dir, pqr_file)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_pdb2pqr(pdb_file, whitespace=whitespace, ff=ff,
            #    parameters=parameters, work_dir=work_dir, docker=False)

    else:
        pqr_file = os.path.join(work_dir, pqr_file)
        command = ["/usr/share/pdb2pqr/pdb2pqr.py"]+parameters
        command += [full_pdb_path, pqr_file]

        try:
            with silence_stdout(), silence_stderr():
                subprocess.call(command)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    assert os.path.isfile(pqr_file)
    return pqr_file

@memory.cache
def run_pdb2pqr_APBS(pdb_path, pdb2pqr_whitespace=False, pdb2pqr_ff="amber",
  pdb2pqr_parameters=None, apbs_input_file=None, apbs_keep_input=False, work_dir=None, docker=True, job=None, clean=True):
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
    pqr_path = run_pdb2pqr(pdb_path, whitespace=pdb2pdb2pqr_whitespace,
        ff=pdb2pqr_ff, parameters=pdb2pqr_parameters, work_dir=work_dir, docker=docker, job=job)
    atom_pot_file = run_apbs(pqr_path, input_file=apbs_input_file,
        keep_input=apbs_keep_input, work_dir=work_dir, docker=docker, job=job)

    result = {}
    with open(pqr_path) as pqr, open(atom_pot_file) as atom_pot:
        #Skip first 4 rows of atompot file
        for _ in xrange(4):
            next(atom_pot)

            for line in pqr:
                if not line.startswith("ATOM  ") or line.startswith("HETATM"): continue

                electrostatic_potential = float(next(atompot))

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

    if clean:
        os.remove(pqr_path)
        os.remove(atom_pot_file)

    return result
