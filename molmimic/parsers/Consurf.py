import os
from tempfile import mkdtemp
import zipfile
from shutil import copyfileobj

import requests

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import natural_keys

from toil.realtimeLogger import RealtimeLogger

from joblib import Parallel, delayed

def download_consurf_scores(pdb, chain, n_tries=3, consurf_path=None):
    if consurf_path is None:
        consurf_path = os.getcwd()

    url = "http://bental.tau.ac.il/new_ConSurfDB/DB/{}/{}/consurf.grades".format(pdb, chain)

    pdb_id = pdb.upper()
    if chain != " ":
        pdb_id += "_"+chain

    consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)

    for _ in range(n_tries):
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

def parse_consurf_line(line, pdb_id=None, consurf_path=None, download_all=False, done_consurf=None):
    if consurf_path is None:
        consurf_path = os.getcwd()
    store = IOStore.get("aws:us-east-1:molmimic-consurf")
    if download_all or (pdb_id is not None and pdb_id in line):
        line = line.strip()

        cluster_rep_id, cluster = line.split(":", 1)
        cluster_rep_pdb, cluster_rep_chain = cluster_rep_id.split("_")
        cluster_pdb_ids = cluster[:-1].split(", ")

        consurf_rep_pdb_id = cluster_rep_pdb
        if cluster_rep_chain != " ":
            consurf_rep_pdb_id += "_"+cluster_rep_chain

        if (done_consurf is not None and consurf_rep_pdb_id in done_consurf):
            return

        consurf_rep_db_file = os.path.join(cluster_rep_pdb[1:3].upper(), consurf_rep_pdb_id+".scores")

        if (done_consurf is not None and consurf_rep_pdb_id in done_consurf) or \
          (done_consurf is None and store.exists(consurf_rep_db_file)):
            print("--Exists")
            consurf_f = os.path.join(consurf_path, consurf_rep_db_file+".scores")
            consurf_exists = True
            if not download_all:
                store.read_input_file(consurf_rep_db_file, consurf_f)
                return consurf_f
        else:
            consurf_exists = False

        consurf_f = download_consurf_scores(cluster_rep_pdb, cluster_rep_chain, consurf_path=consurf_path)
        if not consurf_f:
            return False

        if download_all:
            print("--Saving", consurf_rep_db_file)
            if not consurf_exists:
                store.write_output_file(consurf_f, consurf_rep_db_file)
                set_consurf = True
            else:
                set_consurf = False
            for pc in cluster[:-1].split(", "):
                if not pc or (done_consurf is not None and pc in done_consurf):
                    continue
                if consurf_exists and not set_consurf:
                    consurf_f = store.read_input_file(consurf_rep_db_file, consurf_f)
                pdb, chain = pc.split("_")
                print("-----Saving", os.path.join(pdb[1:3].upper(), pc))
                store.write_output_file(consurf_f, os.path.join(pdb[1:3].upper(), pc+".scores"))

        os.remove(consurf_f)

def download_consurf(pdb=None, chain=None, consurf_path=None):
    if consurf_path is None:
        consurf_path = os.getcwd()

    if (pdb, chain).count(None) == 0:
        pdb_id = pdb.upper()
        if chain != " ": pdb_id += "_"+chain
        consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)
        download_all = False
    elif (pdb, chain).count(None) == 2:
        download_all = True
        pdb_id = None

    store = IOStore.get("aws:us-east-1:molmimic-consurf")
    pdb_list = os.path.join(consurf_path, "pdbaa_list.nr")

    done_consurf = [os.path.splitext(os.path.basename(k))[0] for k in \
        store.list_input_directory() if k != "pdbaa_list.nr"]

    if not store.exists("pdbaa_list.nr"):
        r = requests.get("http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip")
        zip_path = os.path.join(consurf_path, "ConSurfDB_list_feature.zip")
        with open(zip_path, "w") as zip:
            zip.write(r.content)
        with zipfile.ZipFile(zip_path) as z:
            with z.open("pdbaa_list.nr") as zf, open(pdb_list, 'wb') as f:
                copyfileobj(zf, f)
        store.write_output_file(pdb_list, "pdbaa_list.nr")
        os.remove(zip_path)
    else:
        store.read_input_file("pdbaa_list.nr", pdb_list)

    with open(pdb_list) as f:
        Parallel(n_jobs=-1)(delayed(parse_consurf_line)(line, pdb_id=pdb_id, \
            consurf_path=consurf_path, download_all=download_all, \
            done_consurf=done_consurf) for line in f)

    os.remove(pdb_list)

def run_consurf(struct, pdb, chain, work_dir=None, download=False):
    if work_dir is None:
        work_dir = os.getcwd()

    pdb_id = pdb.upper()
    if chain != " ":
        pdb_id += "_"+chain
    consurf_db_key = os.path.join(pdb[1:3].upper(), pdb_id+".scores")
    consurf_db_file = os.path.join(work_dir, pdb_id+".scores")

    store = IOStore.get("aws:us-east-1:molmimic-consurf")

    consurf_result = {}
    if store.exists(consurf_db_key):
        store.read_input_file(consurf_db_key, consurf_db_file)
    else:
        if download:
            consurf_db_file = download_consurf(pdb, chain)
            if not consurf_db_file:
                #Error no ConSurf
                return consurf_result
        else:
            return consurf_result

    try:
        with open(consurf_db_file) as f:
            pass
    except IOError as e:
        #Might be empty
        RealtimeLogger.info("Failed reading, {} bc {}".format(consurf_db_file, e))
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
                except IndexError as e:
                    RealtimeLogger.info("Error parsing line {}".format(e))
                    break

                consurf_result[resseq] = score

                # residue = struct.get_residue_from_resseq(resseq)
                #
                # if residue is not None:
                #     for a in residue:
                #         consurf_result[a.get_id()] = score
    return consurf_result
