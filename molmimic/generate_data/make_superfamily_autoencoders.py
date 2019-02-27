import sys
import os
import argparse
import time
from itertools import groupby
import glob

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job
from molmimic.torch_model.torch_train import train

def start_sae(job, sfam_id):
    work_dir = job.fileStore.getLocalTempDir()
    pdb_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
    feature_store = IOStore.get("aws:us-east-1:molmimic-structure-features")
    sae_store = IOStore.get("aws:us-east-1:molmimic-structure-features")

    extensions = set(["atom.npy", "residue.npy", "edges.gz"])
    done_files = lambda k: set([f.rsplit("_", 1)[1] for f in \
        out_store.list_input_directory(k)])
    pdb_keys = [k for k in pdb_store.list_input_directory(str(int(sfam_id))) if \
        k.endswith(".pdb") and extensions != done_files(os.path.splitext(k)[0])]

    structure_path = os.path.join(work_dir, "structures")
    features_path = os.path.join(work_dir, "features")

    data = []

    sfam_id = str(int(sfam_id))
    for feature_key in feature_store.list_input_directory(sfam_id):
        if feature_key.endswith("atom.npy"):
            data_key = feature_key[:-8]

            feature_file = os.path.join(work_dir, "features", feature_key[len(sfam_id):])
            if not os.path.isdir(os.path.dirname(feature_file)):
                os.makedirs(os.path.dirname(feature_file))
            feature_store.read_input_file(feature_key, feature_file)

            pdb_key = data_key+".pdb"
            pdb_file = os.path.join(work_dir, "structures", pdb_key[len(sfam_id):])
            if not os.path.isdir(os.path.dirname(pdb_file)):
                os.makedirs(os.path.dirname(pdb_file))
            feature_store.read_input_file(pdb_key, pdb_file)

            pdb, chain, sdi, domNo = os.path.basename(pdb_key).split("_")
            sdi, domNo = sdi[3:], domNo[1:]

            data.append([data_key, pdb, chain, sdi, domNo, pdb_key, feature_key])

    data_key = os.path.join(sfam_id, "autoencoder_keys.h5")
    data_file = os.path.join(work_dir, "autoencoder_keys.h5" )
    df = pd.DataFrame(data, names=["key", "pdb", "chain", "sdi", "domNo",
        "structure", "features"])
    df.to_hdf(data_file, "table")
    sae_store.write_output_file(data_file, data_key)
    del df
    del data

    def callback(epoch, statsfile, graphs, model_file, epoch_file):
        key = os.path.join(sfam_id, epoch)
        sae_store.write_output_file(statsfile, os.path.join(key, statsfile))
        sae_store.write_output_file(model_file, os.path.join(key, model_file))
        sae_store.write_output_file(epoch_file, os.path.join(key, epoch_file))
        for graph in graphs:
            sae_store.write_output_file(graph, os.path.join(key, graph))


    train(data_file, model_prefix="{}_sae".format(sfam_id), num_epochs=100,
        use_resnet_unet=True, nFeatures=None, dropout_depth=True,
        dropout_width=True, dropout_p=0.6, autoencoder=True,
        checkpoint_callback=callback)

def start_toil(job):
    import pandas as pd
    work_dir = job.fileStore.getLocalTempDir()
    # in_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    #
    # pdb_file = os.path.join(work_dir, "PDB.h5")
    # in_store.read_input_file("PDB.h5", pdb_file)
    #
    # sfams = pd.read_hdf(pdb_file, "Superfamilies", columns=
    #     ["sfam_id"]).drop_duplicates().dropna()["sfam_id"].sort_values()

    sfams = [299845.0]

    map_job(job, calculate_features_for_sfam, sfams)

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
    with Toil(options) as toil:
        toil.start(job)
