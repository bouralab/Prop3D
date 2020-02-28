from __future__ import print_function
import sys
import os
import argparse
import time
import traceback
from itertools import groupby
import glob

from molmimic.common.featurizer import ProteinFeaturizer

from molmimic.util import safe_remove
from molmimic.util.iostore import IOStore
from molmimic.util.pdb import InvalidPDB, get_atom_lines
from molmimic.util.hdf import get_file, filter_hdf, filter_hdf_chunks
from molmimic.util.toil import map_job
from molmimic.util.cath import run_cath_hierarchy, download_cath_domain

from molmimic.generate_data import data_stores

from toil.realtimeLogger import RealtimeLogger
from botocore.exceptions import ClientError

class CalculateFeaturesError(RuntimeError):
    def __init__(self, cath_domain, stage, message, errors=None, *args, **kwds):
        super().__init__(*args, **kwds)
        self.cath_domain = cath_domain
        self.stage = stage
        self.message = message
        self.errors = errors if isinstance(errors, list) else []

    def __str__(self):
        return "Error during {}: {}\nErrors:\n".format(self.stage, self.message,
            "\n".join(map(str, self.errors)))

    def save(self, store=None):
        if store is None:
            store = data_stores.cath_features
        fail_file = "{}.{}".format(self.cath_domain, self.stage)
        with open(fail_file, "w") as f:
            print(self.message, file=f)
            print(self.errors, file=f)
        store.write_output_file(fail_file,
            "errors/"+os.path.basename(fail_file))
        safe_remove(fail_file)

def calculate_features(job, cath_domain, cathcode, update_features=None, work_dir=None):
    if work_dir is None:
        if job is not None and hasattr(job, "fileStore"):
            work_dir = job.fileStore.getLocalTempDir()
        else:
            work_dir = os.getcwd()

    cath_key = "{}/{}".format(cathcode, cath_domain)

    if update_features is not None:
        #Download existing features
        for ext in ("atom.h5", "residue.h5", "edges.txt.gz"):
            try:
                data_stores.cath_features.read_input_file(
                    "{}_{}.h5".format(cath_key, ext),
                    os.path.join(work_dir, "{}_{}.h5".format(cath_domain, ext)))
            except ClientError:
                #Ignore, just recalculate
                update_features = None
                break

    domain_file = os.path.join(work_dir, "{}.pdb".format(cath_domain))

    try:
        data_stores.prepared_cath_structures.read_input_file(
            cath_key+".pdb", domain_file)
    except ClientError:
        raise
        RealtimeLogger.info("Failed to download prapred cath file {}".format(
            cath_key+".pdb"))

    structure = ProteinFeaturizer(
        domain_file, cath_domain, job, work_dir,
        force_feature_calculation=update_features is None,
        update_features=update_features)

    for ext, calculate in (("atom.h5", structure.calculate_flat_features),
                           ("residue.h5", structure.calculate_flat_residue_features),
                           ("edges.txt.gz", structure.calculate_graph)):
        try:
            _, feature_file = calculate()
        except (SystemExit, KeyboardInterrupt):
            raise
        except Exception as e:
            tb = traceback.format_exc()
            CalculateFeaturesError(cath_domain, ext.split(".",1)[0], tb).save()
            return
        data_stores.cath_features.write_output_file(
            feature_file, "{}_{}".format(cath_key, ext))
        safe_remove(feature_file)
        RealtimeLogger.info("Finished features for: {}".format(feature_file))

def calculate_features_for_sfam(job, sfam_id, update_features, further_parallelize=True, use_cath=True):
    work_dir = job.fileStore.getLocalTempDir()

    RealtimeLogger.info("Running SFAM {}".format(sfam_id))

    extensions = set([u'atom.h5', u'residue.h5', u'edges.h5'])

    def done_files(k, o):
        return set(f.rsplit("_", 1)[1] for f in \
            o.list_input_directory(k) if "fail" not in f)

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
