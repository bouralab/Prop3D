from __future__ import print_function
import sys
import os
import argparse
import time
from itertools import groupby
import glob

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS

from toil.realtimeLogger import RealtimeLogger

def calculate_features(job, pdb_or_key, sfam_id=None, chain=None, sdi=None, domNo=None, work_dir=None):
    from molmimic.common.featurizer import ProteinFeaturizer

    if work_dir is None and job is not None:
        work_dir = job.fileStore.getLocalTempDir()

    if work_dir is None or not os.path.isdir(work_dir):
        work_dir = os.getcwd()

    if [sfam_id, chain, sdi, domNo].count(None) == 0:
        #pdb_or_key is pdb
        in_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
        out_store = IOStore.get("aws:us-east-1:molmimic-features")

        pdb = pdb_or_key
        key = "{}/{}/{}_{}_sdi{}_d{}".format(int(sfam_id), pdb.lower()[1:3],
            pdb.upper(), chain, sdi, domNo)
    elif pdb_or_key.count("_") == 3:
        #pdb_or_key is mmdb sdi key
        in_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
        out_store = IOStore.get("aws:us-east-1:molmimic-features")

        use_cath = False
        key = os.path.splitext(pdb_or_key)[0]
        pdb, chain, sdi, domNo = os.path.basename(key).split("_")
        sdi, domNo = sdi[3:], domNo[1:]
        sfam_id = pdb_or_key.rsplit("/", 1)[0]

    else:
        in_store = IOStore.get("aws:us-east-1:molmimic-cath-structures")
        out_store = IOStore.get("aws:us-east-1:molmimic-cath-features")

        use_cath = True
        key = os.path.splitext(pdb_or_key)[0]
        sdi = os.path.basename(key)
        pdb, chain, domNo = sdi[:4].lower(), sdi[4], sdi[5:]
        sfam_id = pdb_or_key.rsplit("/", 1)[0]

    RealtimeLogger.info("Running {}".format(key))

    try:
        pdb_path = os.path.join(work_dir, os.path.basename(key)+".pdb")
        in_store.read_input_file(key+".pdb", pdb_path)

        s = ProteinFeaturizer(pdb_path, pdb, chain, sdi, domNo, job, work_dir)

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
    except Exception as e:
        RealtimeLogger.info("Failed to get features for {}: {} - {}".format(os.path.basename(key), type(e), e))
        fail_key = "{}_error.fail".format(key)
        fail_file = os.path.join(work_dir, os.path.basename(key))
        with open(fail_file, "w") as f:
            f.write("{}\n{}\n".format(type(e), w))
        out_store.write_output_file(fail_file, fail_key)
        os.remove(fail_file)


def calculate_features_for_sfam(job, sfam_id, further_parallelize=False, use_cath=True):
    work_dir = job.fileStore.getLocalTempDir()
    pdb_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
    out_store = IOStore.get("aws:us-east-1:molmimic-features")

    extensions = set(["atom.h5", "residue.h5", "edges.gz"])
    done_files = lambda k: set([f.rsplit("_", 1)[1] for f in \
        out_store.list_input_directory(k)])

    if use_cath:
        pdb_store = IOStore.get("aws:us-east-1:molmimic-cath-structures")
        out_store = IOStore.get("aws:us-east-1:molmimic-features-cath")
        pdb_keys = [k for k in pdb_store.list_input_directory(sfam_id) if \
            k.endswith(".pdb") and \
            extensions != done_files(os.path.splitext(k)[0])][:1]
    else:
        pdb_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
        out_store = IOStore.get("aws:us-east-1:molmimic-features")
        pdb_keys = [k for k in pdb_store.list_input_directory(
            str(int(sfam_id))) if k.endswith(".pdb") and \
            extensions != done_files(os.path.splitext(k)[0])]

    if further_parallelize:
        map_job(job, calculate_features, pdb_keys)
    else:
        for pdb_key in pdb_keys: #pdb_store.list_input_directory(int(sfam_id)):
            calculate_features(job, pdb_key, work_dir=work_dir)
    #     except (SystemExit, KeyboardInterrupt):
    #         raise
    #     except Exception as e:
    #         fail_key = "{}_error.fail".format(os.path.splitext(pdb_key)[0])
    #         fail_file = os.path.join(work_dir, os.path.basename(fail_key))
    #         with open(fail_file, "w") as f:
    #             f.write("{}\n".format(e))
    #         out_store.write_output_file(fail_file, fail_key)
    #         os.remove(fail_file)


def start_toil(job, further_parallelize=False, use_cath=True):
    import pandas as pd
    work_dir = job.fileStore.getLocalTempDir()


    if use_cath:
        in_store = IOStore.get("aws:us-east-1:molmimic-cath")
        sfam_file = os.path.join(work_dir, "cath.h5")
        # in_store.read_input_file("cath-domain-description-file-small.h5", sfam_file)
        # sfams = pd.read_hdf(sfam_file, "table", columns=["cathcode"])["cathcode"]
        # sfams = sfams.str.replace(".", "/")
        sfams = ["2/60/40/10"]
    else:
        in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
        sfam_file = os.path.join(work_dir, "PDB.h5")
        in_store.read_input_file("PDB.h5", sfam_file)

        sfams = pd.read_hdf(sfam_file, "Superfamilies", columns=
            ["sfam_id"]).drop_duplicates().dropna()["sfam_id"].sort_values()

    RealtimeLogger.info("Running {} SFAMs".format(len(sfams)))
    RealtimeLogger.info("{}".format(sfams))

    # sfams = [299845.0]

    map_job(job, calculate_features_for_sfam, sfams, further_parallelize, use_cath)

    # try:
    #     os.remove(sfam_file)
    # except (FileNotFoundError, OSError):
    #     pass

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
