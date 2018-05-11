import os, sys
sys.path.append("/usr/share/pdb2pqr/src")

from psize import Psize

import subprocess

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
