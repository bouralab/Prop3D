import os

import prepare_protein
from util import get_first_chain

CUSTOM_DATABASE = "/data/draizene/molmimic/pdb2pqr_structures/pdbs/"

def get_pdb(pdb, chain=None, prepare=True):
    """Get the correct PDB file from custom database of protonated and minimized 
    structures, find the raw PDB file in a local spapshot and prepare the file, 
    or download the file from wwwPDB if there is no snapshot and prepare it.

    Parameters
    ----------
    pdb : str
        Either the 4 letter PDB code or path to PDB file to prepare
    chain : str or None
        Chain ID to extract since molmimic can only one chain at a time. Default 
        (None) is to prepare all chains, but only use first.
    prepare : bool
        Prepare the structure by protonating with pdb2pqr and minimizing using
        rosetta

    Returns
    -------
    path : str
        Path to prepared pdb file
    pdb : str
        PDB ID
    chain : str
        Chain ID
    input_format : str
        The format of PDB file
    """
    if not os.path.isfile(pdb) and len(pdb) == 4:
        #First try custom database of split chains, protonated, and minimized

        path = "{}_{}".format(pdb.lower(), chain)
        custom_db_path = os.path.join(CUSTOM_DATABASE, pdb[1:3].lower(), path)
        
        if os.path.isfile(custom_db_path+".min.pdb"):
            #Get chain split, protinated, and minimized structure
            custom_db_path += ".min.pdb"
            return custom_db_path, pdb, chain, "pdb"

        elif os.path.isfile(custom_db_path+".pqr.pdb"):
            #Get chain split, and protonated structure; shouldn't be here
            custom_db_path = ".pqr.pdb"
            return custom_db_path, pdb, chain, "pqr"

        elif os.path.isfile(custom_db_path+".pdb"):
            #Get chain split structure; shouldn't be here
            custom_db_path += ".pdb"
            return custom_db_path, pdb, chain, "pdb"

        elif os.path.isfile(path+".min.pdb"):
            #Get chain split, protinated, and minimized structure
            path += ".min.pdb"
            return path, pdb, chain, "pdb"

        elif os.path.isfile(path+".pqr.pdb"):
            #Get chain split, and protinated structure; shouldn't be here
            path += ".pqr.pdb"
            return path, pdb, chain, "pqr"

        elif os.path.isfile(path+".pdb"):
            #Get chain split structure; shouldn't be here
            path += ".pdb"
            return path, pdb, chain, "pdb"

        else:
            #Structure not seen before, will split chains, protonate, and minimized (only once)
            path = "{}/pdb/{}/pdb{}.ent.gz".format(os.environ.get("PDB_SNAPSHOT", "/pdb"), pdb[1:3].lower(), pdb.lower())

            if not os.path.isfile(path):
                path, input_format = download_pdb(pdb)
            else:
                input_format = "pdb"

            if prepare:
                path = prepare_protein.run_protein(path, chain=chain)

            return path, pdb.lower(), input_format
    else:
        path = pdb

        if input_format not in ["pdb", "pqr", "mmcif", "mmtf"]:
            input_format = "pdb"

        if prepare:
            if chain is None:
                path, pdb, chain = prepare_protein.run_protein(path)[0]
            else:
                path, pdb, chain = prepare_protein.run_protein(path, chain=chain)
        else:
            pdb, input_format = os.path.splitext(path)
            chain = get_first_chain(path)
        
        return path, pdb, chain, input_format
