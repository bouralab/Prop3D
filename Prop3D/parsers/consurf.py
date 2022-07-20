from __future__ import print_function
import os
import zipfile
import shutil

import requests

import numpy as np

from Prop3D.util.iostore import IOStore
from Prop3D.util import natural_keys
from Prop3D.parsers.json import WebService
from Prop3D.parsers.pdbe import PDBEApi

from toil.realtimeLogger import RealtimeLogger

from joblib import Parallel, delayed

class ConSurfApi(WebService):
    def __init__(self, consurf_store, work_dir=None, download=True, max_attempts=2):
        super(ConSurfApi, self).__init__("https://consurfdb.tau.ac.il/scripts/downloadPDB.php?",
            consurf_store, work_dir=work_dir, download=download, max_attempts=max_attempts)

    def _fix_key(self, key):
        return "/".join([kv.split("=")[1] for kv in key.split("&")])

    def extension(self, key):
        return ".consurf"

    def get_conservation_score(self, pdb, chain, bioPDB, raw_values=False):
        consurf = self.get("pdb_ID={}&view_chain={}".format(pdb, chain))
        try:
            scores = [int(s) if s!="." else np.nan for s in consurf["seq3d_grades"]]
        except KeyError:
            print(consurf)
            raise
        scores = {r.get_id():score for r, score in zip(bioPDB[0][chain].get_residues(), scores)}

        if not raw_values:
            scores = {k:(v-1)/8 for k, v in scores.items()}

        return scores

    def should_add_to_store(self, key, fname):
        """Modify file to only include consurf values, remove PDB lines"""
        newf = fname+".tmp"
        with open(newf, "w") as new:
            for line in self.parse_raw_line(fname):
                print(line, file=new)

        try:
            os.remove(fname)
        except OSError:
            pass

        shutil.move(newf, fname)

        return True

    def parse(self, file_path, key, raw_lines=False):
        def _get_value(value, force_string=False):
            if value.startswith("Array"):
                return list(map(int, value[6:-2].split(",")))

            if value.startswith('"') and value.endswith('"'):
                #String
                value = value[1:-1]
                if force_string:
                    return value

            try:
                return int(value)
            except ValueError:
                try:
                    return float(value)
                except ValueError:
                    return value

        consurf = {}
        raw_lines = iter(self.parse_raw_line(file_path))
        for line in raw_lines:
            fields = line[2:-1].split("=")
            if len(fields) ==  2:
                #Key/Value Pairs
                key, value = fields[0].strip(), fields[1].strip()
                value = _get_value(value)
            else:
                #Key, values on next lines
                key = fields[0].strip()
                value = None

                stop = False
                for loop_line in raw_lines:
                    loop_line = loop_line.rstrip()
                    if loop_line.endswith(";"):
                        stop = True
                    if loop_line.startswith("!"):
                        loop_line = loop_line[2:-2]

                    if value is None:
                        value = _get_value(loop_line, force_string=True)
                    else:
                        value += _get_value(loop_line, force_string=True)

                    if stop:
                        break

            consurf[key] = value

        return consurf

    def parse_raw_line(self, file_path):
        parse = False
        with open(file_path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith("! "):
                    parse = True
                    yield line
                if parse and line.startswith("!!"):
                    break

class ConSurfApiOld(WebService):
    def __init__(self, consurf_store, pdbe_store, work_dir=None, download=True,
      use_representative_chains=True, max_attempts=2):
        super(ConSurfApi, self).__init__("http://bental.tau.ac.il/new_ConSurfDB/",
            consurf_store, work_dir=work_dir, download=download, max_attempts=max_attempts)
        self.use_representative_chains = use_representative_chains
        self.pdbe_api = PDBEApi(pdbe_store, work_dir=work_dir, download=download,
            max_attempts=max_attempts)
        self.clusters, self.cluster_members = self._get_clusters()

    def should_add_to_store(self, key, fname):
        if key == "ConSurfDB_list_feature.zip":
            return False
        return True

    def get_conservation_score(self, pdb, chain):
        member_id = (pdb.upper(), chain)
        try:
            #Get representative
            rep_pdb_id = self.cluster_members[member_id]
        except KetError:
            raise KeyError("Cannot get scores for {}_{}".format(pdb, chain))

        rep_scores = self.get("{}/{}/consurf.grades".format(pdb, chain))

        if self.use_representative_chains or rep_pdb_id == member_id:
            return rep_scores

        member_residues = self.pdbe_api.get_pdb_residues(pdb, chain)

        member_scores = pd.merge(rep_scores, member_residues, on="residueNumber",
            suffixes=["_rep", "_mem"])
        member_scores = member_scores.drop(columns=["pdbResidueNumber_rep"]).rename(
            columns={"pdbResidueNumber_mem":"pdbResidueNumber"})

        return member_scores

    def _get_clusters(self):
        for _ in range(3):
            try:
                clusters = self.get("pdbaa_list.nr")
                break
            except KeyError:
                try:
                    self.get("ConSurfDB_list_feature.zip")
                    break
                except KeyError as e:
                    raise
        else:
            raise RuntimeError("Cannot download pdbaa_list.nr")

        cluster_members = {member:cluster for cluster, members in clusters.items() \
            for member in members}
        cluster_members.update({cluster:cluster for cluster in clusters.keys()})

        return clusters, cluster_members

    def parse(self, file_name, key):
        if key == "pdbaa_list.nr":
            return self._parse_clusters(file_name)

        elif key == "ConSurfDB_list_feature.zip":
            return self._decompress_clusters(file_name)

        elif key.endswith(".scores"):
            return self._parse_scores_file(file_name)

    def download_all(self, n_jobs=-1):
        pdbs = cluster_members.keys() if self.use_representative_chains else \
            self.clusters.keys()
        Parallel(n_jobs=n_jobs)(delayed(get_conservation_score)(pdb, chain) \
            for pdb, chain in pdbs)

    def _parse_clusters(self, file_name):
        with open(file_name) as fh:
            return dict(self._parse_consurf_line(line) for line in fh)

    def _decompress_clusters(self, zip_path):
        pdb_list = os.path.join(self.work_dir, "pdbaa_list.nr")
        with open(zip_path, "w") as zip:
            zip.write(r.content)
        with zipfile.ZipFile(zip_path) as z:
            with z.open("pdbaa_list.nr") as zf, open(pdb_list, 'wb') as f:
                copyfileobj(zf, f)
        self.store.write_output_file(pdb_list, "pdbaa_list.nr")

        clusters = self._parse_clusters(pdb_list, "pdbaa_list.nr")

        try:
            os.remove(zip_path)
            del self.files[zip_path]
        except OSError:
            pass

        return clusters

    def _parse_scores_file(self, file_name):
        parsing = False

        consurf_result = {"residueNum":[], "pdbResidueNum":[]}

        with open(file_name) as fh:
            for i, line in enumerate(fh):
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
                        resNum = int(fields[0])
                    except ValueError as e:
                        RealtimeLogger.info("Error parsing line {}".format(e))
                        continue

                    try:
                        pdbResNum = fields[2].split(":")[0][3:]
                    except IndexError as e:
                        RealtimeLogger.info("Error parsing line {}".format(e))
                        continue

                    consurf_result["residueNum"].append(resNum)
                    consurf_result["pdbResidueNum"].append(pdbResNum)

                    # try:
                    #     resseq_parts = natural_keys(fields[2].split(":")[0][3:])
                    #     resseq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
                    #     score = (float(fields[3]), int(fields[4].replace("*", "")))
                    # except IndexError as e:
                    #     RealtimeLogger.info("Error parsing line {}".format(e))
                    #     continue

                    # try:
                    #     resseq_parts = natural_keys(fields[2].split(":")[0][3:])
                    #     resseq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
                    #     score = (float(fields[3]), int(fields[4].replace("*", "")))
                    # except IndexError as e:
                    #     RealtimeLogger.info("Error parsing line {}".format(e))
                    #     continue

        return pd.DataFrame(consurf_result)

    def _parse_consurf_line(self, line):
        print(line)
        cluster_rep_id, cluster_members = line.rstrip().split(":", 1)
        cluster_rep_id = tuple(cluster_rep_id.split("_"))
        cluster_member_ids = [tuple(member) for member in \
            cluster.replace(".", "").split(", ")]

        #Don't think this is necesarsy?
        # consurf_rep_pdb_id = cluster_rep_pdb
        # if cluster_rep_chain != " ":
        #     consurf_rep_pdb_id += "_"+cluster_rep_chain

        return consurf_rep_pdb_id, cluster_member_ids

def get_consurf_feature():
    if consurf and not hasattr(self, "_consurf"):
        #self._consurf = run_consurf(self.pdb, self.chain, work_dir=self.work_dir)
        consurf_store = IOStore.get("aws:us-east-1:Prop3D-consurf-service")
        consurf_api = ConSurfApi(consurf_store, work_dir=self.work_dir)
        self._consurf = consurf_api.get_conservation_score(self.pdb,
            self.chain, self.path)

    if consurf:
        result["consurf_normalized_score"] = self._consurf.get(residue.get_id(), np.nan)
        result["consurf_is_conserved"] = float(consurf_normalized_score>=.6)



# def download_consurf_scores(pdb, chain, n_tries=3, consurf_path=None):
#     if consurf_path is None:
#         consurf_path = os.getcwd()
#
#     url = "http://bental.tau.ac.il/new_ConSurfDB/DB/{}/{}/consurf.grades".format(pdb, chain)
#
#     pdb_id = pdb.upper()
#     if chain != " ":
#         pdb_id += "_"+chain
#
#     consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)
#
#     for _ in range(n_tries):
#         try:
#             r = requests.get(url)
#             break
#         except requests.exceptions.ConnectionError as e:
#             RealtimeLogger.info("ConSurf Error {}: {}".format(type(e), e))
#         except requests.exceptions.HTTPError as e:
#             RealtimeLogger.info("ConSurf Error {}: {}".format(type(e), e))
#             return None
#
#     if not os.path.exists(os.path.dirname(consurf_db_file)):
#         os.makedirs(os.path.dirname(consurf_db_file))
#     with open(consurf_db_file, "w") as f:
#         print(r.content, file=f)
#     return consurf_db_file
#
# def parse_consurf_line(line, pdb_id=None, consurf_path=None, download_all=False, done_consurf=None):
#     if consurf_path is None:
#         consurf_path = os.getcwd()
#     store = IOStore.get("aws:us-east-1:Prop3D-consurf")
#     if download_all or (pdb_id is not None and pdb_id in line):
#         line = line.strip()
#
#         cluster_rep_id, cluster = line.split(":", 1)
#         cluster_rep_pdb, cluster_rep_chain = cluster_rep_id.split("_")
#         cluster_pdb_ids = cluster[:-1].split(", ")
#
#         consurf_rep_pdb_id = cluster_rep_pdb
#         if cluster_rep_chain != " ":
#             consurf_rep_pdb_id += "_"+cluster_rep_chain
#
#         if (done_consurf is not None and consurf_rep_pdb_id in done_consurf):
#             return
#
#         consurf_rep_db_file = os.path.join(cluster_rep_pdb[1:3].upper(), consurf_rep_pdb_id+".scores")
#
#         if (done_consurf is not None and consurf_rep_pdb_id in done_consurf) or \
#           (done_consurf is None and store.exists(consurf_rep_db_file)):
#             print("--Exists")
#             consurf_f = os.path.join(consurf_path, consurf_rep_db_file+".scores")
#             consurf_exists = True
#             if not download_all:
#                 store.read_input_file(consurf_rep_db_file, consurf_f)
#                 return consurf_f
#         else:
#             consurf_exists = False
#
#         consurf_f = download_consurf_scores(cluster_rep_pdb, cluster_rep_chain, consurf_path=consurf_path)
#         if not consurf_f:
#             return False
#
#         if download_all:
#             print("--Saving", consurf_rep_db_file)
#             if not consurf_exists:
#                 store.write_output_file(consurf_f, consurf_rep_db_file)
#                 set_consurf = True
#             else:
#                 set_consurf = False
#             for pc in cluster[:-1].split(", "):
#                 if not pc or (done_consurf is not None and pc in done_consurf):
#                     continue
#                 if consurf_exists and not set_consurf:
#                     consurf_f = store.read_input_file(consurf_rep_db_file, consurf_f)
#                 pdb, chain = pc.split("_")
#                 print("-----Saving", os.path.join(pdb[1:3].upper(), pc))
#                 store.write_output_file(consurf_f, os.path.join(pdb[1:3].upper(), pc+".scores"))
#
#         return consurf_f

# def download_consurf(pdb=None, chain=None, consurf_path=None, download_all=False):
#     if consurf_path is None:
#         consurf_path = os.getcwd()
#
#     if (pdb, chain).count(None) == 0:
#         pdb_id = pdb.upper()
#         if chain != " ": pdb_id += "_"+chain
#         consurf_db_file = os.path.join(consurf_path, pdb[1:3].upper(), pdb_id)
#         download_all = False
#     else:
#         #download_all = True
#         pdb_id = None
#
#     store = IOStore.get("aws:us-east-1:Prop3D-consurf")
#     pdb_list = os.path.join(consurf_path, "pdbaa_list.nr")
#
#     # done_consurf = [os.path.splitext(os.path.basename(k))[0] for k in \
#     #     store.list_input_directory() if k != "pdbaa_list.nr"]
#
#     if download_all:
#         done_consurf = [os.path.splitext(os.path.basename(k))[0] for k in \
#             store.list_input_directory() if k != "pdbaa_list.nr"]
#     else:
#         done_consurf = []
#
#     if not store.exists("pdbaa_list.nr"):
#         r = requests.get("http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip")
#         zip_path = os.path.join(consurf_path, "ConSurfDB_list_feature.zip")
#         with open(zip_path, "w") as zip:
#             zip.write(r.content)
#         with zipfile.ZipFile(zip_path) as z:
#             with z.open("pdbaa_list.nr") as zf, open(pdb_list, 'wb') as f:
#                 copyfileobj(zf, f)
#         store.write_output_file(pdb_list, "pdbaa_list.nr")
#         os.remove(zip_path)
#     else:
#         store.read_input_file("pdbaa_list.nr", pdb_list)
#
#     with open(pdb_list) as f:
#         if download_all:
#             consurf_f = None
#             Parallel(n_jobs=-1)(delayed(parse_consurf_line)(line, pdb_id=pdb_id, \
#                 consurf_path=consurf_path, download_all=download_all, \
#                 done_consurf=done_consurf) for line in f)
#         else:
#             for line in f:
#                 if pdb_id in line:
#                     consurf_f = parse_consurf_line(line, pdb_id=pdb_id, consurf_path=consurf_path,
#                         download_all=False, done_consurf=done_consurf)
#                     break
#             else:
#                 return
#
#     try:
#         os.remove(pdb_list)
#     except OSError:
#         pass
#
#     return consurf_f
#
# def run_consurf(pdb, chain, work_dir=None, download=True):
#     if work_dir is None:
#         work_dir = os.getcwd()
#
#     pdb_id = pdb.upper()
#     if chain != " ":
#         pdb_id += "_"+chain
#     consurf_db_key = os.path.join(pdb[1:3].upper(), pdb_id+".scores")
#     consurf_db_file = os.path.join(work_dir, pdb_id+".scores")
#
#     store = IOStore.get("aws:us-east-1:Prop3D-consurf")
#
#     consurf_result = {}
#
#     if os.path.isfile(consurf_db_file):
#         pass
#     elif store.exists(consurf_db_key):
#         store.read_input_file(consurf_db_key, consurf_db_file)
#     else:
#         if download:
#             consurf_db_file = download_consurf(pdb, chain)
#             if not consurf_db_file:
#                 #Error no ConSurf
#                 return consurf_result
#         else:
#             return consurf_result
#
#     try:
#         with open(consurf_db_file) as f:
#             pass
#     except IOError as e:
#         #Might be empty
#         RealtimeLogger.info("Failed reading, {} bc {}".format(consurf_db_file, e))
#         return consurf_result
