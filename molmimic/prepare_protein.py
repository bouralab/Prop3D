import os, sys
import subprocess
import shutil
import gzip
import glob
import re

from util import get_first_chain

PDB_PATH = CUSTOM_DATABASE = "/data/draizene/molmimic/pdb2pqr_structures/pdbs/"
#os.environ.get("MOLMIMIC_PDB_PATH", os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdbs"))

def run_protein(pdb_file, chain=None):
    name_format = "^([A-Za-z0-9]{4}).pdb"
    match = re.match(name_format, pdb_file)
    if match:
        pdb = match.group(1)
        pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
    else:
        print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
        pdb_path = os.getcwd()

    #Cleanup PDB file and add TER lines in between gaps
    if pdb_file.endswith(".gz"):
        with gzip.open(pdb_file, 'rt') as f:
            cleaned_pdb = subprocess.Popen(["pdb_tidy.py"], stdin=f, stdout=subprocess.PIPE)
        file_prefix = os.path.basename(os.path.splitext(pdb_file[:-2])[0]) 
    else:
        cleaned_pdb = subprocess.Popen(["pdb_tidy.py", pdb_file], stdout=subprocess.PIPE)
        file_prefix = os.path.basename(os.path.splitext(pdb_file)[0])

    #Remove HETATMS
    pdb_file = os.path.join(pdb_path, "{}.pdb".format(file_prefix))
    with open(pdb_file, "w") as pdb:
        subprocess.Popen(["pdb_striphet.py"], stdin=cleaned_pdb, stdout=pdb)

    if chain is None:
        #Split PDB into chains, 1 chain per file
        subprocess.Popen(["pdb_splitchain.py", pdb_file])

        #Process all chains
        return [run_single_chain(chain_file) for chain_file in \
            glob.glob(os.path.join(pdb_path, "{}_*.pdb".format(file_prefix)))]

    else:
        #Split desired chain in PDB into 1 file
        chain_file = os.apth.join(pdb_path, "{}_{}.pdb".format(file_prefix, chain))
        with open(chain_file, "w") as chain:
            subprocess.Popen(["pdb_chain.py", "-{}".format(chain), pdb_file], stdout=chain)
        return run_single_chain(chain_file)

def run_single_chain(pdb_file):
    assert not pdb_file.endswith(".gz"), "Cannot be a gzip archive, try 'run_protein' instead"
    assert subprocess.check_output(["pdb_wc.py", "-c", pdb_file]).split()[0] == "1", "Must contain only one chain"

    file_prefix = os.path.splitext(pdb_file)[0]

    name_format = "^([A-Za-z0-9]{4})_([A-Za-z0-9]{1,3}).pdb$"
    match = re.match(name_format, pdb_file)
    if match:
        pdb, chain = match.group(1), match.group(2)
        pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
    else:
        print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
        pdb_path = os.getcwd()
        pdb = os.path.basename(file_prefix)
        chain = get_first_chain(pdb_file)
        

    #Remove altLocs
    delocc_file = "{}.delocc".format(pdb_file)
    with open(delocc_file, "w") as delocc:
        subprocess.Popen(["pdb_tidy.py", pdb_file], stdout=delocc)

    #Add hydrogens
    pqr_file = "{}.pqr.pdb".format(file_prefix)
    propka_file = "{}.proka".format(pqr_file)
    subprocess.Popen(["pdb2pqr", "--ff=parse", "--ph-calc-method=propka", "--chain", "--drop-water", delocc_file, pqr_file])

    #rename propKa file
    shutil.move(propka_file, "{}.propka".format(file_prefix))

    #Minimize
    score_file = "{}.sc".format(file_prefix)
    minimized_file = "{}.pqr_0001.pdb".format(file_prefix)
    subprocess.Popen(["minimize", 
        "-s", pqr_file,
        "-run:min_type", "lbfgs_armijo_nonmonotone",
        "-run:min_tolerance", "0.001",
        "-ignore_zero_occupancy", "false",
        "-out:file:scorefile", score_file,
        "-out:path:pdb", pdb_path,
        "-out:path:score", pdb_path])

    #Cleanup Rosetta PDB
    cleaned_minimized_file = "{}.min.pdb".format(file_prefix)
    with open(min_file, "w") as cleaned_minimized:
        subprocess.Popen(["pdb_stripheader.py", minimized_file], stdout=cleaned_minimized)

    return cleaned_minimized_file, pdb, chain

if __name__ == "__main__":
    pass