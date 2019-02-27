import json
import subprocess
from joblib import Memory
from Bio.PDB import Structure
from molmimic.util import silence_stdout, silence_stderr

try:
	from toil.lib.docker import apiDockerCall
	freesasa = None
except ImportError:
	apiDockerCall = None

	try:
		import freesasa
	except ImportError:
		freesasa = None
		if not subprocess.check_output(["which", "freesasa"]):
			raise RuntimeError("Cannot use FreeSASA because it is not in the users path, nor is Docker enabled")

memory = Memory(verbose=0)

from molmimic.util import silence_stdout, silence_stderr

class FreeSASAException(Exception):
	pass

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

@memory.cache
def run_freesasa(pdb, job=None, parameters=None):
	if isinstance(pdb, str) and os.path.isfile(pdb) and job is not None:
		if apiDockerCall is not None:
			try:
				run_freesasa_docker(pdb, job, parameters)
			except (SystemExit, KeyboardInterrupt):
				raise
			except:
				run_freesasa_subprocess(pdb, parameters)
		else:
			run_freesasa_subprocess(pdb, parameters)
	if isinstance(pdb, str) and os.path.isfile(pdb) and job is None:
		run_freesasa_subprocess(pdb, parameters)
	elif isinstance(pdb, Structure):
		run_freesasa_biopython(pdb)
	else:
		raise FreeSASAException("Unknown pdb {}".format(pdb))
