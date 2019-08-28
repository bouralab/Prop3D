import os, sys

try:
	from toil.lib.docker import apiDockerCall
except ImportError:
	apiDockerCall = None
	import subprocess

def run_usearch(parameters, work_dir=None, docker=True, job=None):
	if work_dir is None:
		work_dir = os.getcwd()

	if docker and apiDockerCall is not None and job is not None:
		_parameters = []
		outfiles = []
		for i, p in enumerate(parameters):
			if p.startswith("{i}"):
				full = p[3:]
				short = os.path.basename(full)
				if not os.path.abspath(os.path.dirname(full)) == os.path.abspath(work_dir):
					shutil.copy(full, work_dir)
				_parameters.append(os.path.join("/data", short))
			elif p.startswith("{o}"):
				full = p[3:]
				short = os.path.basename(full)
				_parameters.append(os.path.join("/data", short))
				outfiles.append(os.path.join(work_dir, short))
			else:
				_parameters.append(p)

		try:
			apiDockerCall(job,
						  image='edraizen/usearch:latest',
						  working_dir="/data",
                          volumes={work_dir:{"bind":"/data", "mode":"rw"}},
						  parameters=_parameters,
						  entrypoint="/opt/usearch/usearch")
		except (SystemExit, KeyboardInterrupt):
			raise
		except:
			raise
	else:
		new_parameters = []
		outfiles = []
		for i, p in enumerate(parameters):
			if p.startswith("{i}"):
				new_parameters.append(p[3:])
			elif p.startswith("{o}"):
				new_parameters.append(p[3:])
				outfiles.append(p[3:])
			else:
				new_parameters.append(p)
		command = ["usearch"]+new_parameters

		try:
			subprocess.check_output(command, stderr=subprocess.PIPE)
		except subprocess.CalledProcessError:
			raise RuntimeError("Unable to minimize file {}".format(pqr_file))

	return outfiles
