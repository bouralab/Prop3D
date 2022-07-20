import os

import prepare_protein
from util import InvalidPDB, download_pdb, get_first_chain

CUSTOM_DATABASE = "{}/../data/structures/Ig2".format(os.path.dirname(__file__))

def get_pdb(pdb, chain=None, sdi=None, domain=None, prepare=True, updating=True, allow_download=False):
    """Get the correct PDB file from custom database of protonated and minimized
    structures, find the raw PDB file in a local spapshot and prepare the file,
    or download the file from wwwPDB if there is no snapshot and prepare it.

    Parameters
    ----------
    pdb : str
        Either the 4 letter PDB code or path to PDB file to prepare
    chain : str or None
        Chain ID to extract since Prop3D can only one chain at a time. Default
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
    domain : 2-tuple
        Domain SDI, Domain Num
    input_format : str
        The format of PDB file
    """
    if not os.path.isfile(pdb) and len(pdb) == 4:
        #First try custom database of split chains, protonated, and minimized
        pdb = pdb.lower()
        path = "{}_{}".format(pdb.lower(), chain)

        if "MOLMIMIC_STRUCTURES" in os.environ and os.path.isdir(os.environ["MOLMIMIC_STRUCTURES"]):
            custom_db_path = os.path.join(os.environ["MOLMIMIC_STRUCTURES"], pdb[1:3].lower(), path)

            if sdi is not None and domain is not None:
                custom_domain_path = "{}_sdi{}_d{}.pdb".format(custom_db_path, sdi, domain)
                if os.path.isfile(custom_domain_path):
                    return custom_domain_path, pdb, chain, (sdi, domain), "pdb"

            if os.path.isfile(custom_db_path+".min.pdb"):
                #Get chain split, protinated, and minimized structure
                custom_min_path = "{}.min.pdb".format(custom_db_path)
                return custom_min_path, pdb, chain, (sdi, domain), "pdb"

            elif os.path.isfile(custom_db_path+".pqr.pdb"):
                #Get chain split, and protonated structure; shouldn't be here
                custom_pqr_path = "{}.pqr.pdb".format(custom_db_path)
                return custom_pqr_path, pdb, chain, (sdi, domain), "pqr"

            elif os.path.isfile(custom_db_path+".pdb"):
                #Get chain split structure; shouldn't be here
                custom_pdb_path = "{}.pdb".format(custom_db_path)
                return custom_pdb_path, pdb, chain, (sdi, domain), "pdb"

        elif os.path.isfile(path+".min.pdb"):
            #Get chain split, protinated, and minimized structure
            path += ".min.pdb"
            return path, pdb, chain, (sdi, domain), "pdb"

        elif os.path.isfile(path+".pqr.pdb"):
            #Get chain split, and protinated structure; shouldn't be here
            path += ".pqr.pdb"
            return path, pdb, chain, (sdi, domain), "pqr"

        elif os.path.isfile(path+".pdb"):
            #Get chain split structure; shouldn't be here
            path += ".pdb"
            return path, pdb, chain, (sdi, domain), "pdb"

        else:
            #Structure not seen before, will split chains, protonate, and minimized (only once)
            path = "{}/pdb/{}/pdb{}.ent.gz".format(os.environ.get("PDB_SNAPSHOT", "/pdb"), pdb[1:3].lower(), pdb.lower())

            if not os.path.isfile(path):
                if allow_download:
                    path, input_format = download_pdb(pdb)
                else:
                    raise InvalidPDB("PDB not found: {}".format(pdb))
            else:
                input_format = "pdb"

            if prepare:
                for x in prepare_protein.run_protein(path, chain=chain):
                    try:
                        path, pdb, chain, (sdi, domain) = x
                    except ValueError:
                        print("XXX", x)
                        raise
                    break
                else:
                    raise RuntimeErorr("Invalid chians")

            return path, pdb.lower(), chain, (sdi, domain), input_format
    else:
        path = pdb

        if input_format not in ["pdb", "pqr", "mmcif", "mmtf"]:
            input_format = "pdb"

        if prepare:
            if chain is None:
                path, pdb, chain, (sdi, domain) = next(prepare_protein.run_protein(path))
            else:
                path, pdb, chain, (sdi, domain) = next(prepare_protein.run_protein(path, chain=chain))
        else:
            pdb, input_format = os.path.splitext(path)
            chain = get_first_chain(path)

        return path, pdb, chain, (sdi, domain), input_format
