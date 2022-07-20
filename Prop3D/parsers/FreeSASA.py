import json
import subprocess
from joblib import Memory
from Bio.PDB import Structure
from Prop3D.util import silence_stdout, silence_stderr

from Prop3D.parsers.container import Container

try:
    from toil.lib.docker import apiDockerCall
    freesasa = None
except ImportError:
    apiDockerCall = None

try:
    import freesasa
except ImportError:
    freesasa = None

memory = Memory(verbose=0)

from Prop3D.util import silence_stdout, silence_stderr

class FreeSASAException(Exception):
    pass

class FreeSASA(Container):
    PARAMETERS = [("pdb_file", "path:in", "format"), (":out_format:json", str, "--format={}")]
    LOCAL = ["freesasa"]

    def __init__(self, bioPDB, chain):
        self.bioPDB = bioPDB
        self.chain = chain
        self.sasa, self.sasa_struct = run_freesasa_biopython(bioPDB)
        self.atom_asa = {}
        self.residue_asa = {}

    def __call__(self, *args, **kwds):
        if freesasa is not None:
            #Try python first
            pdb_path = self.format_parameters(args, kwds)[0]
            with silence_stdout(), silence_stderr():
                sasa_struct = freesasa.Structure(pdb_path)
                sasa = freesasa.calc(sasa_struct)
        else:
            super(self, FreeSASA).__call__(*args, **kwds)

        return

    def parse(self, freesasa_out):
        freesasa = "{"+freesasa_out.split("\n",1)[1].strip().replace("-nan", "NaN").replace("nan", "NaN")

        try:
            return json.loads(freesasa)
        except:
            raise FreeSASAException("Unable to freesasa convert to JSON: {}".format(freesasa))

    def get_atom_rasa(self, atom, acc_threshold=.05):
        """Get an atom's relative accescibale surface area.
        """
        atom_asa = self.get_atom_asa(atom)
        residue_asa = self.get_residue_asa(atom)
        atom_rasa = atom_asa/residue_asa
        return atom_rasa, atom_rasa <= acc_threshold

    def get_atom_asa(self, atom):
        atom_id = [
            self.chain,
            "".join(atom.get_parent().get_id()[1:]), #residue number and insertion code
            atom.get_id()[0] #Atom name
        ]

        if atom_id in self.atom_asa:
            return self.atom_asa[atom_id]

        try:
            selection = "sele, chain {} and resi {} and name {}".format(*atom_id)
            with silence_stdout(), silence_stderr():
                selections = freesasa.selectArea([selection], sasa_struct, sasa)
                atom_area = selections["sele"]
                self.atom_asa[atom_id] = atom_area
                return atom_area
        except (KeyError, AssertionError, AttributeError, TypeError):
            raise
            self.atom_asa[atom_id] = np.NaN
            return np.NaN

    def get_residue_asa(self, atom_or_residue):
        if isinstance(atom_or_residue, PDB.Atom.Atom):
            is_atom = True
            atom = atom_or_residue
            residue = atom_or_residue.get_parent()
        elif isinstance(atom_or_residue, PDB.Residue.Residue):
            is_atom = False
            residue = atom_or_residue
        else:
            raise RuntimeError("Input must be Atom or Residue")

        residue_id = [
            self.chain,
            "".join(atom.get_parent().get_id()[1:]), #residue number and insertion code
        ]

        if residue_id in self.residue_asa:
            return self.residue_asa[residue_id]

        try:
            selection = "sele, chain {} and resi {}".format(*residue_id)
            with silence_stdout(), silence_stderr():
                selections = freesasa.selectArea([selection], sasa_struct, sasa)
                residue_area = selections["sele"]
            self.residue_asa[residue_id] = residue_area
            return residue_area
        except (KeyError, AssertionError, AttributeError, TypeError):
            raise
            self.residue_asa[residue_id] = np.Nan
            return np.NaN


@memory.cache
def run_freesasa_biopython(pdb_path):
    global freesasa
    if freesasa is None:
        try:
            import freesasa
        except ImportError:
            raise RuntimeError("Cannot use this method. Please save the pdb file and rerun with docker")

    with silence_stdout(), silence_stderr():
        #Automatically removes hydrogens
        sasa_struct = freesasa.Structure(pdb_path)
        sasa = freesasa.calc(sasa_struct)

    return sasa, sasa_struct

def run_freesasa_subprocess(pdb_file, parameters=None, format="json"):
    assert format in ("json", None)

    parameters = parameters if isinstance(parameters, (list, tuple)) else []
    if format=="json":
        parameters.append("--format=json")
        parameters.append(pdb_file)

    FNULL = open(os.devnull, 'w')
    try:
        freesasa = subprocess.check_output(command, stderr=FNULL)
    except subprocess.CalledProcessError:
        raise FreeSASAException("Failed running free sasa: {}".format(command))

    return parse_freesasa(freesasa)

def run_freesasa_docker(pdb_file, job, parameters=None, format="json"):
    assert apiDockerCall is not None
    assert format in ("json", None)

    work_dir = job.fileStore.getLocalTempDir()
    full_pdb_path = pdb_file
    pdb_path = os.path.basename(full_pdb_path)

    _parameters = ["freesasa"]
    if isinstance(parameters, (list, tuple)):
        _parameters += parameters

    if format=="json":
        _parameters.append("--format=json")
        _parameters.append("/data/{}".format(pdb_path))

    try:
        freesasa = apiDockerCall(job,
          image='edraizen/freesasa:latest',
          working_dir=work_dir,
          parameters=parameters)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as e:
        raise FreeSASAException("Unable to run docker {}".format(e))

    return parse_freesasa(freesasa)

def parse_freesasa(freesasa_out):
    freesasa = "{"+freesasa_out.split("\n",1)[1].strip().replace("-nan", "NaN").replace("nan", "NaN")

    try:
        return json.loads(freesasa)
    except:
        raise FreeSASAException("Unable to freesasa convert to JSON: {}".format(freesasa))

# @memory.cache
# def run_freesasa(pdb, job=None, parameters=None):
#     if isinstance(pdb, str) and os.path.isfile(pdb) and job is not None:
#     if apiDockerCall is not None:
#     try:
#     run_freesasa_docker(pdb, job, parameters)
#     except (SystemExit, KeyboardInterrupt):
#     raise
#     except:
#     run_freesasa_subprocess(pdb, parameters)
#     else:
#     run_freesasa_subprocess(pdb, parameters)
#     if isinstance(pdb, str) and os.path.isfile(pdb) and job is None:
#     run_freesasa_subprocess(pdb, parameters)
#     elif isinstance(pdb, Structure):
#     run_freesasa_biopython(pdb)
#     else:
#     raise FreeSASAException("Unknown pdb {}".format(pdb))
