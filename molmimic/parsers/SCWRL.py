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

def run_scwrl(pdb_file, output_prefix=None, framefilename=None, sequencefilename=None,
  paramfilename=None, in_cystal=False, remove_hydrogens=False, remove_h_n_term=False,
  work_dir=None, job=None):
    if work_dir is None:
        work_dir = os.getcwd()

    full_pdb_path = pdb_file
    pdb_path = os.path.basename(full_pdb_path)
    output_file = "{}.pqr".format(pdb_path)

    if output_prefix is None:
        output_prefix = os.path.splitext(full_pdb_path)[0]

    output_file = "{}.scwrl.pdb".format(output_prefix)

    _parameters = [p for p, use in [("-#", in_cystal), ("-h", remove_hydrogens), \
        ("-t", remove_h_n_term)] if use]

    if docker and apiDockerCall is not None and job is not None:
        #Docker can only read from work_dir
        if not os.path.abspath(os.path.dirname(pdb_file)) == os.path.abspath(work_dir):
            shutil.copy(pdb_file, work_dir)

        parameters = ["-i", "/data/{}".format(os.path.basename(pdb_file))]
        parameters += ["-o", "/data/{}".format(os.path.basename(output_file))]
        parameters += _parameters

        if framefilename is not None and os.path.isfile(framefilename):
            if not os.path.abspath(os.path.dirname(framefilename)) == os.path.abspath(work_dir):
                shutil.copy(framefilename, work_dir)
            parameters += ["-f", "/data/{}".format(os.path.basename(framefilename))]

        if sequencefilename is not None and os.path.isfile(sequencefilename):
            if not os.path.abspath(os.path.dirname(sequencefilename)) == os.path.abspath(work_dir):
                shutil.copy(sequencefilename, work_dir)
            parameters += ["-s", "/data/{}".format(os.path.basename(sequencefilename))]

        if paramfilename is not None and os.path.isfile(paramfilename):
            if not os.path.abspath(os.path.dirname(paramfilename)) == os.path.abspath(work_dir):
                shutil.copy(paramfilename, work_dir)
            parameters += ["-p", "/data/{}".format(os.path.basename(paramfilename))]

        try:
            apiDockerCall(job,
                          image='edraizen/pdb2pqr:latest',
                          working_dir=work_dir,
                          parameters=parameters)
            pqr_file = os.path.join(work_dir, pqr_file)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            raise
            #return run_pdb2pqr(pdb_file, whitespace=whitespace, ff=ff,
            #    parameters=parameters, work_dir=work_dir, docker=False)

        if not os.path.abspath(os.path.dirname(output_file)) == os.path.abspath(work_dir):
            shutil.move(os.path.join(work_dir, os.path.basename(output_file)),
                os.path.abspath(os.path.dirname(output_file))

    else:
        parameters = ["scwrl4", "-i", pdb_file, "-o", output_file]+_parameters
        if framefilename is not None and os.path.isfile(framefilename):
            parameters += ["-f", framefilename]
        if sequencefilename is not None and os.path.isfile(sequencefilename):
            parameters += ["-s", sequencefilename]
        if paramfilename is not None and os.path.isfile(paramfilename):
            parameters += ["-p", paramfilename]

        try:
            subprocess.call(parameters)
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            raise
            #raise RuntimeError("APBS failed becuase it was not found in path: {}".format(e))

    assert os.path.isfile(output_file)
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
