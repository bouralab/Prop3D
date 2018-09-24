import os
from tempfile import mkdtemp
from shutil import copyfile

import requests
from joblib import Memory

from molmimic.util import data_path_prefix

cachedir = mkdtemp()
memory = Memory(cachedir=cachedir, verbose=0)

CONSURF_PATH = os.path.join(data_path_prefix, "ConSurf")

def download_consurf_scores(pdb, chain, n_tries=3, consurf_path=CONSURF_PATH):
    url = "http://bental.tau.ac.il/new_ConSurfDB/DB/{}/{}/consurf.grades".format(pdb, chain)

    pdb_id = pdb.upper()
    if chain != " ": pdb_id += "_"+chain
    consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)

    for _ in xrange(n_tries):
        try:
            r = requests.get(url)
            break
        except requests.exceptions.ConnectionError:
            pass
        except requests.exceptions.HTTPError:
            return None

    if not os.path.exists(os.path.dirname(consurf_db_file)):
        os.makedirs(os.path.dirname(consurf_db_file))
    with open(consurf_db_file, "w") as f:
        print >> f, r.content
    return consurf_db_file

def download_consurf(pdb=None, chain=None, consurf_path=CONSURF_PATH):
    if (pdb, chain).count(None) == 0:
        pdb_id = pdb.upper()
        if chain != " ": pdb_id += "_"+chain
        consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)
        download_all = False
    elif (pdb, chain).count(None) == 2:
        download_all = True
        pdb_id = None

    with open(os.path.join(consurf_path, "pdbaa_list.nr")) as f:
        for line in f:
            if download_all or pdb_id in line:
                line = line.strip()
                cluster_rep_id, cluster = line.split(":", 1)
                cluster_rep_pdb, cluster_rep_chain = cluster_rep_id.split("_")

                if download_all:
                    consurf_pdb_id = cluster_rep_pdb
                    if chain != " ": consurf_pdb_id += "_"+cluster_rep_chain
                    consurf_db_file = os.path.join(consurf_path, cluster_rep_pdb[1:3].upper(), consurf_pdb_id)

                if os.path.isfile(consurf_db_file):
                    return True

                consurf_f = download_consurf_scores(cluster_rep_pdb, cluster_rep_chain)
                if not consurf_f:
                    if not download_all:
                        return False
                    continue

                if download_all:
                    dir = os.path.dirname(os.path.dirname(consurf_f))
                    for pc in cluster.split(", "):
                        if not pc: continue
                        pdb, chain = pc.split("_")
                        copyfile(consurf_f, os.path.join(dir, pdb[1:3].upper(), pc))

                return True
    return False

@memory.cache
def run_consurf(struct, pdb, chain):
    pdb_id = pdb.upper()
    if chain != " ": pdb_id += "_"+chain
    consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)
    print consurf_db_file

    consurf_result = {}

    if not os.path.isfile(consurf_db_file):
        rt = download_consurf(pdb, chain)
        if not rt:
            #Error no ConSurf
            return consurf_result

    try:
        with open(consurf_db_file) as f:
            pass
    except IOError:
        #Might be empty
        return consurf_result

    parsing = False
    with open(consurf_db_file) as f:
        for line in f:
            if not parsing and line.strip().startswith("(normalized)"):
                parsing = True
                continue

            if parsing:
                fields = line.rstrip().split()

                if len(fields) == 0:
                    break
                elif fields[2] == "-":
                    continue

                try:
                    resseq_parts = natural_keys(fields[2].split(":")[0][3:])
                    resseq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
                    score = (float(fields[3]), int(fields[4].replace("*", "")))
                except IndexError:
                    break

                residue = struct.get_residue_from_resseq(resseq)

                if residue is not None:
                    for a in residue:
                        consurf_result[a.get_id()] = score
    return consurf_result
