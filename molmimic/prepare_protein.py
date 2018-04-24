import os, sys
sys.path.append("/data/draizene/pdb-tools")
sys.path.append("/data/draizene/molmimic")

import subprocess
import shutil
import gzip
import glob
import re
from itertools import groupby

from molmimic.calculate_features import SwarmJob
from util import get_first_chain, get_all_chains

PDB_PATH = CUSTOM_DATABASE = "/data/draizene/molmimic/pdb2pqr_structures/full_H2B"
#os.environ.get("MOLMIMIC_PDB_PATH", os.path.join(os.path.dirname(os.path.dirname(__file__)), "pdbs"))

def run_protein(pdb_file, chain=None):
    """Prepare a protein structure for use in molmimic. This will
    0) Unzip if gzipped
    1) Cleanup PDB file and add TER lines in between gaps
    2) Remove HETATMS
    3) Split Chains into separate files
    4) Each chain is then protonated and minimized. See 'run_single_chain'
    for more info

    Parameters
    ----------
    pdb_file : str
        Path to PDB file
    chain : str or None
        Chain ID to split out, protonate, and mininize. If None (default),
        all chains will be split out, protonated, and mininized.

    Returns
    -------
    If no chain was specified a list of all chains is returned with the
    path to the prepared file, the pdb code, and chain. If a chain is
    specified, 3 values are returned: the path to the prepared file,
    the pdb code, and chain. See 'run_single_chain' for info.
    """
    print pdb_file
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))


    base = os.path.basename(pdb_file)
    if base.startswith("pdb") and base.endswith(".ent.gz"):
        name_format = "^pdb([A-Za-z0-9]{4}).ent.gz"
    else:
        name_format = "^([A-Za-z0-9]{4}).pdb"

    match = re.match(name_format, base)
    if match:
        pdb = match.group(1)
        pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
    else:
        print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
        pdb_path = os.getcwd()
        pdb = os.path.basename(os.path.splitext(pdb_file.replace(".gz", ""))[0])

    #Unzip PDB file
    if pdb_file.endswith(".gz"):
        unzipped_pdb = os.path.join(pdb_path, "{}.raw.pdb".format(pdb))
        print unzipped_pdb
        with gzip.open(pdb_file, 'rt') as f, open(unzipped_pdb, "wb") as out_file:
            d = f.read()
            #print d
            print >> out_file, d
        pdb_file = unzipped_pdb

    #Cleanup PDB file and add TER lines in between gaps
    cleaned_pdb = subprocess.Popen(["/data/draizene/pdb-tools/pdb_tidy.py", pdb_file], stdout=subprocess.PIPE)

    #Remove HETATMS
    pdb_file = os.path.join(pdb_path, "{}.pdb".format(pdb))
    with open(pdb_file, "w") as pdbf:
        strip_hetatm = subprocess.Popen(["/data/draizene/pdb-tools/pdb_striphet.py"], stdin=subprocess.PIPE, stdout=pdbf)
        strip_hetatm.communicate(cleaned_pdb.stdout.read())

    if chain is None:
        print "run all chains"
        #Split PDB into chains, 1 chain per file
        subprocess.call(["/data/draizene/pdb-tools/pdb_splitchain.py", pdb_file])

        #Process all chains
        if SwarmJob.slurm_enabled():
            for chain in get_all_chains(pdb_file):
                chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
                yield chain_file
        else:
            for chain in get_all_chains(pdb_file):
                chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
                yield run_single_chain(chain_file)


    else:
        #Split desired chain in PDB into 1 file
        chain_file = os.path.join(pdb_path, "{}_{}.pdb".format(pdb, chain))
        with open(chain_file, "w") as chain:
            subprocess.call(["/data/draizene/pdb-tools/pdb_chain.py", "-{}".format(chain), pdb_file], stdout=chain)
        yield run_single_chain(chain_file)

def run_single_chain(pdb_file, calculate_features=True):
    """Prepare a single chain for use in molmimic (called during 'run_protein').
    This prepares the chains by:
    (1) Removing all altLoc's except the first
    (2) Add hydrogens using pdb2pqr (ff=parse, ph=propka)
    (3) Minimzed usinf rosetta (lbfgs_armijo_nonmonotone with tolerance 0.001)
    (4) Cleaned so that can be parsed by simple PDB parsers

    Note: raises an AssertError if pdb_file is gzipped or contains more than
    one chain

    Parameters
    ----------
    pdb_file : str
        Path to PDB file of chain to prepare. Must not be gzipped or contain
        more than one chain.

    Returns
    -------
    cleaned_minimized_file : str
        Path to prepared PDB file
    pdb : str
        PDB ID
    chain : str
        Chain ID
    """
    print "Running", pdb_file
    if not os.path.isfile(pdb_file):
        raise RuntimeError("Invalid PDB File, cannot find {}".format(pdb_file))

    assert not pdb_file.endswith(".gz"), "Cannot be a gzip archive, try 'run_protein' instead"
    assert len(get_all_chains(pdb_file)) == 1, "Must contain only one chain"

    file_prefix = os.path.splitext(pdb_file)[0]
    base = os.path.basename(pdb_file)
    name_format = "^([A-Za-z0-9]{4})_([A-Za-z0-9]{1,3}).pdb$"
    match = re.match(name_format, base)
    if match:
        pdb, chain = match.group(1), match.group(2)
        pdb_path = os.path.join(PDB_PATH, pdb[1:3].lower())
    else:
        print >> sys.stderr, "Invalid PDB Name, results saved to current working directory"
        pdb_path = os.getcwd()
        pdb = os.path.basename(file_prefix)
        chain = get_first_chain(pdb_file)


    #Remove altLocs
    delocc_file = os.path.join(pdb_path, "{}.delocc".format(pdb_file))
    with open(delocc_file, "w") as delocc:
        subprocess.Popen(["/data/draizene/pdb-tools/pdb_tidy.py", pdb_file], stdout=delocc)

    #Add hydrogens
    pqr_file = os.path.join(pdb_path, "{}.pqr.pdb".format(file_prefix))
    propka_file = os.path.join(pdb_path, "{}.propka".format(pqr_file))
    subprocess.call(["/data/draizene/3dcnn-torch-py2", "shell", "pdb2pqr", "--ff=parse", "--ph-calc-method=propka", "--chain", "--drop-water", delocc_file, pqr_file])

    #rename propKa file
    shutil.move(propka_file, "{}.propka".format(file_prefix))

    if not os.path.isfile(pqr_file):
        raise RuntimeError("Unable to protonate {} using pdb2pqr. Please check pdb2pqr error logs. Most likeley reason for failing is that the structure is missing too many heavy atoms.".format(pdb_file))

    #Minimize
    score_file = "{}.sc".format(file_prefix)
    minimized_file = "{}.pqr_0001.pdb".format(file_prefix)
    subprocess.call(["minimize",
        "-s", pqr_file,
        "-run:min_type", "lbfgs_armijo_nonmonotone",
        "-run:min_tolerance", "0.001",
        "-ignore_zero_occupancy", "false",
        "-out:file:scorefile", score_file,
        "-out:path:pdb", pdb_path,
        "-out:path:score", pdb_path])

    #Cleanup Rosetta PDB
    cleaned_minimized_file = "{}.min.pdb".format(file_prefix)
    with open(cleaned_minimized_file, "w") as cleaned_minimized:
        subprocess.Popen(["/data/draizene/pdb-tools/pdb_stripheader.py", minimized_file], stdout=cleaned_minimized)

    if len(get_all_chains(cleaned_minimized_file)) > 1:
        #Sometimes the chains are split in pdb2pqr, fixes it to one chain
        one_chain_clean_file = "{}.min.one_chain.pdb".format(file_prefix)
        with open(one_chain_clean_file, "w") as one_chain_clean:
            subprocess.Popen(["/data/draizene/pdb-tools/pdb_chain.py", "-{}".format(chain), cleaned_minimized_file], stdout=one_chain_clean)

        return_file = one_chain_clean_file

    else:
        return_file = cleaned_minimized_file

    if calculate_features:
        subprocess.call(["/data/draizene/3dcnn-torch-py2",
            "python", os.path.realpath(__file__), pdb, chain, "None", "False"])

    return return_file, pdb, chain

def load_ibis(ibis_data, minimize=True, cellular_organisms=False):
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data, cellular_organisms=cellular_organisms)
    data = dataset.data

    key = lambda p: p[1:3]
    pdbs = sorted(set(data["pdb"].values.tolist()), key=key)

    job = SwarmJob("prepare_chains", modules="rosetta")

    i = 0
    for pdb_divided, pdb_group in groupby(pdbs, key=key):
        pdb_path = "/data/draizene/molmimic/pdb2pqr_structures/full_H2B/{}".format(pdb_divided.lower())
        feature_path = "/data/draizene/molmimic/features_full_H2B/atom/{}".format(pdb_divided.lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
        if not os.path.exists(feature_path):
            os.makedirs(feature_path)
        for pdb in pdb_group:
            i+=1
            #if len(glob.glob(os.path.join(pdb_path, "{}_*.min.pdb".format(pdb.lower())))) > 0:
            #    continue
            pdb_file = "/pdb/pdb/{}/pdb{}.ent.gz".format(pdb_divided.lower(), pdb.lower())
            try:
                for chain_file in run_protein(pdb_file):
                    job += "python {} single-chain {}\n".format(__file__, chain_file)
            except RuntimeError:
                continue

    job.run()

def submit_ibis(ibis_data):
    job = SwarmJob("prepare_proteins", walltime="96:00:00", individual=True)
    job.add_individual_parameters()
    job += "python {} load {}\n".format(__file__, ibis_data)
    print job.submit_individual()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        submit_ibis(sys.argv[1])
    elif len(sys.argv) == 3 and "load" in sys.argv:
        load_ibis(sys.argv[2])
    elif len(sys.argv) in [3,4] and "protein" in sys.argv:
        chain = sys.argv[3] if len(sys.argv)==4 else None
        run_protein(sys.argv[2], chain)
    elif len(sys.argv) == 3 and "single-chain" in sys.argv:
        run_single_chain(sys.argv[2])
