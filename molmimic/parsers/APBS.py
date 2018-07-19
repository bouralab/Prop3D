import os, sys
from tempfile import mkdtemp
import subprocess

pdb2pqr_src = os.path.join(os.path.dirname(subprocess.check_output(["which", "pdb2pqr"])), "src")
sys.path.append(pdb2pqr_src)

from psize import Psize
from joblib import Memory

from molmimic.util import silence_stdout, silence_stderr

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

@memory.cache
def run_apbs(pqr_file):
	"""Run APBS. Calculates correct size using Psize and defualt from Chimera
	"""
	file_prefix = os.path.splitext(pqr_file)[0]
	input_file = "{}.apbs_input".format(file_prefix)
	output_prefix = "{}.apbs_output".format(file_prefix)

	ps = Psize()
	ps.runPsize(pqr_file)
	cglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getCoarseGridDims())
	fglen = "{:.2f} {:.2f} {:.2f}".format(*ps.getFineGridDims())
	dime = "{:d} {:d} {:d}".format(*ps.getFineGridPoints())

	with open(input_file, "w") as f:
		print >> f, """read
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

	subprocess.call(["apbs", input_file])

	os.remove(input_file)

	return output_prefix+".txt"

@memory.cache
def run_pdb2pqr_APBS(struct, modified=False):
	"""Run PDB2PQR and APBS to get charge and electrostatics for each atom for each atom
	"""
	remove_pqr = True

	if modified:
		pdbfd, tmp_pdb_path = tempfile.mkstemp()
		with os.fdopen(pdbfd, 'w') as tmp:
			struct.save_pdb(tmp, True)
	else:
		tmp_pdb_path = struct.path

	_, tmp_pqr_path = tempfile.mkstemp()

	with silence_stdout(), silence_stderr():
		subprocess.call(["/usr/share/pdb2pqr/pdb2pqr.py", "--ff=amber", "--whitespace", tmp_pdb_path, tmp_pqr_path])

	if modified:
		os.remove(tmp_pdb_path)

	atompot_file = run_apbs(tmp_pqr_path)

	result = {}
	with open(tmp_pqr_path) as pqr, open(atompot_file) as atompot:
		#Skip first 4 rows of atompot file
		for _ in xrange(4):
			next(atompot)

		for line in pqr:
			for line in pqr:
				if not line.startswith("ATOM  ") or line.startswith("HETATM"): continue

				electrostatic_potential = float(next(atompot))

				fields = line.rstrip().split()
				if len(fields) == 11:
					recordName, serial, atomName, residueName, chainID, residueNumber, X, Y, Z, charge, radius = fields
				elif len(fields) == 10:
					recordName, serial, atomName, residueName, residueNumber, X, Y, Z, charge, radius = fields
				else:
					print fields
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

	if remove_pqr:
		os.remove(tmp_pqr_path)

	return result
