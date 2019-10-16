from __future__ import print_function
import sys
import os
import argparse
import time
from itertools import groupby
import glob

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import InvalidPDB, get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS
from molmimic.generate_data.parse_cath import download_cath_domain

from toil.realtimeLogger import RealtimeLogger
from botocore.exceptions import ClientError

def redownload_cath_if_needed(job, pdb_path, sdi, sfam_id, error, work_dir):
    RealtimeLogger.info("Handle Error {} {} {}".format(type(error), error,
        isinstance(error, InvalidPDB)))
    if "Error get chains" in str(error):
        RealtimeLogger.info("MUST RE-DOWNLOAD PDB FILE {}".format(sdi))
        with open(pdb_path) as f:
            for line in f:
                if "<!DOCTYPE html>" in line:
                    break
            else:
                return False
        sfam_id = map(int, sfam_id)
        hierachy = ["class", "architechture", "topology", "homology"]
        cath_domain = dict(zip(hierachy, sfam_id))
        cath_domain["cath_domain"] = sdi
        download_cath_domain(job, cath_domain, work_dir=work_dir)
        return True
    else:
        return False

def calculate_features(job, pdb_or_key, update_features=None, work_dir=None):
    from molmimic.common.featurizer import ProteinFeaturizer

    if work_dir is None and job is not None:
        work_dir = job.fileStore.getLocalTempDir()

    if work_dir is None or not os.path.isdir(work_dir):
        work_dir = os.getcwd()

    # if [sfam_id, chain, sdi, domNo].count(None) == 0:
    #     #pdb_or_key is pdb
    #     in_store = IOStore.get("aws:us-east-1:molmimic-full-structure")
    #     out_store = IOStore.get("aws:us-east-1:molmimic-features")
    #
    #     pdb = pdb_or_key
    #     key = "{}/{}/{}_{}_sdi{}_d{}".format(int(sfam_id), pdb.lower()[1:3],
    #         pdb.upper(), chain, sdi, domNo)
    # el
    if pdb_or_key.count("_") == 3:
        #pdb_or_key is mmdb sdi key
        in_store = IOStore.get("aws:us-east-1:molmimic-full-structure")
        out_store = IOStore.get("aws:us-east-1:molmimic-features")

        use_cath = False
        key = os.path.splitext(pdb_or_key)[0]
        pdb, chain, sdi, domNo = os.path.basename(key).split("_")
        sdi, domNo = sdi[3:], domNo[1:]
        sfam_id = pdb_or_key.rsplit("/", 1)[0]

    else:
        in_store = IOStore.get("aws:us-east-1:molmimic-cath-structure")
        out_store = IOStore.get("aws:us-east-1:molmimic-cath-features")

        use_cath = True
        key = os.path.splitext(pdb_or_key)[0]
        sdi = os.path.basename(key)
        pdb, chain, domNo = sdi[:4].lower(), sdi[4], sdi[5:]
        sfam_id = pdb_or_key.rsplit("/", 1)[0]

        if update_features is not None:
            for ext in ("atom.h5", "residue.h5", "edges.h5"):
                try:
                    store.read_input_file("{}_{}.h5".format(key, ext),
                        "{}_{}.h5".format(sdi, ext))
                except ClientError:
                    pass

    RealtimeLogger.info("Running {} {} {}".format(key, sdi, sfam_id))

    save_error = None

    for _ in range(3):
        try:
            pdb_path = os.path.join(work_dir, os.path.basename(key)+".pdb")
            in_store.read_input_file(key+".pdb", pdb_path)

            if update_features is None:
                s = ProteinFeaturizer(pdb_path, pdb, chain, sdi, domNo, job, work_dir,
                                      force_feature_calculation=True)
            else:
                s = ProteinFeaturizer(pdb_path, pdb, chain, sdi, domNo, job, work_dir,
                                      update_features=update_features)

            _, atom_features = s.calculate_flat_features()
            RealtimeLogger.info("Finished atom features")
            _, residue_features = s.calculate_flat_features(course_grained=True)
            RealtimeLogger.info("Finished residue features")
            graph_features = s.calculate_graph()
            RealtimeLogger.info("Finished edge features")

            out_store.write_output_file(atom_features, key+"_atom.h5")
            out_store.write_output_file(residue_features, key+"_residue.h5")
            out_store.write_output_file(graph_features, key+"_edges.h5")

            for f in (pdb_path, atom_features, residue_features, graph_features):
                try:
                    os.remove(f)
                except OSError:
                    pass
        except (SystemExit, KeyboardInterrupt):
            raise
        except InvalidPDB as e:
            rc = redownload_cath_if_needed(job, pdb_path, sdi, sfam_id.split("/"), e, work_dir)
            if not rc:
                save_error = e
        except Exception as e:
            rc = redownload_cath_if_needed(job, pdb_path, sdi, sfam_id.split("/"), e, work_dir)
            if not rc:
                save_error = e
                RealtimeLogger.info("Error not caught PDB FILE {}???".format(sdi))
                import traceback
                RealtimeLogger.info(traceback.format_exc())

        break

    if save_error is not None:
        RealtimeLogger.info("Failed to get features for {}: {} - {}".format(os.path.basename(key), type(e), e))
        fail_key = "{}_error.fail".format(key)
        fail_file = os.path.join(work_dir, os.path.basename(key))
        with open(fail_file, "w") as f:
            f.write("{}\n{}\n".format(type(e), e))
        out_store.write_output_file(fail_file, fail_key)
        os.remove(fail_file)

def calculate_features_for_sfam(job, sfam_id, update_features, further_parallelize=True, use_cath=True):
    work_dir = job.fileStore.getLocalTempDir()

    RealtimeLogger.info("Running SFAM {}".format(sfam_id))

    extensions = set([u'edges.h5', u'residue.h5', u'atom.h5']) #set(["atom.h5", "residue.h5", "edges.gz"])

    # done_files = lambda k, o: set([f.rsplit("_", 1)[1] for f in \
    #     o.list_input_directory(k)])

    def done_files(k, o):
        file_types = set(f.rsplit("_", 1)[1] for f in \
            o.list_input_directory(k) if "fail" not in f)
        #RealtimeLogger.info("File types for {}: {}".format(k, file_types))
        return file_types

    if use_cath:
        pdb_store = IOStore.get("aws:us-east-1:molmimic-cath-structure")
        out_store = IOStore.get("aws:us-east-1:molmimic-cath-features")
        pdb_keys_full = set(k for k in pdb_store.list_input_directory(sfam_id) \
            if k.endswith(".pdb"))

        if update_features is None:
            pdb_keys_done = set(k for k in pdb_keys_full if \
                done_files(os.path.splitext(k)[0], out_store)==extensions)
            pdb_keys = list(pdb_keys_full-pdb_keys_done)
        else:
            pdb_keys = list(pdb_keys_full)

        RealtimeLogger.info("RUNNING {}/{} DOMAINS from {}: {}".format(len(pdb_keys),
            len(pdb_keys_full), sfam_id, pdb_keys))

        # pdb_keys = [k for k in pdb_store.list_input_directory(sfam_id) if \
        #     k.endswith(".pdb") and \
        #     extensions != done_files(os.path.splitext(k)[0], out_store)]
    else:
        pdb_store = IOStore.get("aws:us-east-1:molmimic-full-structure")
        out_store = IOStore.get("aws:us-east-1:molmimic-features")
        pdb_keys = [k for k in pdb_store.list_input_directory(
            str(int(sfam_id))) if k.endswith(".pdb") and \
            extensions != done_files(os.path.splitext(k)[0], out_store)]

        RealtimeLogger.info("RUNNING {} DOMAINS from {}: {}".format(len(pdb_keys), sfam_id, pdb_keys))

    if further_parallelize:
        map_job(job, calculate_features, pdb_keys, update_features)
    else:
        for pdb_key in pdb_keys: #pdb_store.list_input_directory(int(sfam_id)):
            try:
                calculate_features(job, pdb_key, update_features, work_dir=work_dir)
            except (SystemExit, KeyboardInterrupt):
                raise
            except Exception as e:
                fail_key = "{}_error.fail".format(os.path.splitext(pdb_key)[0])
                fail_file = os.path.join(work_dir, os.path.basename(fail_key))
                with open(fail_file, "w") as f:
                    f.write("{}\n".format(e))
                out_store.write_output_file(fail_file, fail_key)
                os.remove(fail_file)

def run_cath_hierarchy(job, cathcode, cathFileStoreID, update_features=None, further_parallelize=True):
    work_dir = job.fileStore.getLocalTempDir()
    cath_path = get_file(job, "cath.h5", cathFileStoreID, work_dir=work_dir)

    cath_names = ["class", "architechture", "topology", "homology"]
    cathcode = dict(zip(cath_names, cathcode))
    cath_names = cath_names[:len(cathcode)+1]

    cathcodes = filter_hdf_chunks(
        cath_path,
        "table",
        columns=cath_names,
        drop_duplicates=True,
        **cathcode)[cath_names]

    RealtimeLogger.info("cathcode {} {} {}".format(cathcode, cath_names, cathcodes.columns))

    if len(cathcodes.columns) < 4:
        map_job(job, run_cath_hierarchy, cathcodes.values,
            cathFileStoreID, update_features, further_parallelize)
        RealtimeLogger.info("Running {} {}s".format(len(cathcodes), cathcodes.columns[-1]))
    else:
        sfams = (cathcodes.astype(int).astype(str)+"/").sum(axis=1).str[:-1].tolist()
        RealtimeLogger.info("Running sfam {}".format(cathcode))
        map_job(job, calculate_features_for_sfam, sfams, update_features,
            further_parallelize, True)

    try:
        os.remove(cath_path)
    except (FileNotFoundError, OSError):
        pass

def start_toil(job, further_parallelize=False, use_cath=True, update_features=None):
    import pandas as pd
    work_dir = job.fileStore.getLocalTempDir()


    if use_cath:
        in_store = IOStore.get("aws:us-east-1:molmimic-cath")
        sfam_file = os.path.join(work_dir, "cath.h5")
        in_store.read_input_file("cath-domain-description-file-small.h5", sfam_file)

        cathFileStoreID = job.fileStore.writeGlobalFile(sfam_file)

        # run_cath_hierarchy(job, (2,130,10), cathFileStoreID)
        #
        # map_job(job, calculate_features, ["1/20/80/40/4jk7A02.pdb"])

        classes = filter_hdf_chunks(sfam_file, "table", columns=["class"],
           drop_duplicates=True).sort_values("class")["class"].values[:, None]
        map_job(job, run_cath_hierarchy, classes, cathFileStoreID, update_features)
        #sfams = ["2/60/40/10"]
    else:
        in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
        sfam_file = os.path.join(work_dir, "PDB.h5")
        in_store.read_input_file("PDB.h5", sfam_file)

        sfams = pd.read_hdf(sfam_file, "Superfamilies", columns=
            ["sfam_id"]).drop_duplicates().dropna()["sfam_id"].sort_values()

        RealtimeLogger.info("Running {} SFAMs".format(len(sfams)))
        RealtimeLogger.info("{}".format(sfams))

        #sfams = [299845.0]

        map_job(job, calculate_features_for_sfam, sfams, further_parallelize, use_cath)

    try:
        os.remove(sfam_file)
    except (FileNotFoundError, OSError):
        pass

    #os.remove(pdb_file)

    #job.addChildJobFn(calculate_features, "301320/yc/1YCS_A_sdi225433_d0.pdb")

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as workflow:
        workflow.start(job)
