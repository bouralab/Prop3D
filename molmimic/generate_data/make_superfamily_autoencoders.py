import sys
import os
import argparse
import time
from itertools import groupby
import glob

import pandas as pd

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job
from molmimic.torch_model.torch_train import train

def normalize(job, features, work_dir=None):
    """
    one-hot: 0:18,20:23,24:26,27:30,31:34,35:37,38:70,72
    raw: 18,19,23,26,30,34,37,71
    normalized: 70
    """
    import numpy as np
    from sklearn import preprocessing

    if work_dir is None and job is not None:
        work_dir = job.fileStore.getLocalTempDir()
    elif work_dir is None:
        work_dir = os.getcwd()

    one_hot = ((0,18), (20,23), (24, 26), (27, 30), (31, 34), (35,37), (38,70), (72,73))
    one_hot = [n for start, stop in one_hot for n in range(start, stop)]
    raw = [18,19,23,26,30,34,37,71]
    normalized = [70]

    standard_scaler = preprocessing.StandardScaler()
    minmax_scaler = preprocessing.MinMaxScaler()

    for feature_f in features:
        print(feature_f)
        features = np.load(feature_f, mode="r")
        print(features.shape)
        features = features.reshape((features.shape[0], 73))
        standard_scaler.partial_fit(features[:,raw])
        minmax_scaler.partial_fit(features[:,raw+normalized])
        features[:,one_hot].save(os.path.splitext(feature_f)[0]+".one_hot.npy")
        del features

    for feature_f in features:
        base = os.path.splitext(feature_f)[0]
        features = np.memmap(feature_f, mode="c")
        std_scaled = standard_scaler.transform(features[:,raw])
        features[:, raw] = std_scaled
        np.save(base+".stdscaled.npy", features)
        minmax_scaled = minmax_scaler.transform(features[:,raw])
        features[:, raw] = minmax_scaled
        np.save(base+".minmaxscaled.npy", features)
        del features

def start_sae(job, sfam_id, work_dir=None):
    if work_dir is None and job is not None:
        work_dir = job.fileStore.getLocalTempDir()
    elif work_dir is None:
        work_dir = os.getcwd()
    # pdb_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
    # feature_store = IOStore.get("aws:us-east-1:molmimic-structure-features")
    sae_store = IOStore.get("aws:us-east-1:molmimic-structure-features")

    # extensions = set(["atom.npy", "residue.npy", "edges.gz"])
    # done_files = lambda k: set([f.rsplit("_", 1)[1] for f in \
    #     feature_store.list_input_directory(k)])
    # pdb_keys = [k for k in pdb_store.list_input_directory(str(int(sfam_id))) if \
    #     k.endswith(".pdb") and extensions != done_files(os.path.splitext(k)[0])]

    # structure_path = os.path.join(work_dir, "structures")
    # features_path = os.path.join(work_dir, "features")

    data = []

    sfam_id = str(int(sfam_id))
    #for feature_key in feature_store.list_input_directory(sfam_id):
        #if feature_key.endswith("atom.npy"):
    for feature_file in glob.glob(os.path.join(work_dir, "features", sfam_id, "*", "*_atom.npy")):
        if True:
            #data_key = feature_key[:-8]
            #feature_file = feature_key #os.path.join(work_dir, "features", feature_key[len(sfam_id):])

            # if not os.path.isdir(os.path.dirname(feature_file)):
            #     os.makedirs(os.path.dirname(feature_file))
            # feature_store.read_input_file(feature_key, feature_file)

            # pdb_file =
            #
            # pdb_key = data_key+".pdb"
            # #pdb_file = os.path.join(work_dir, "structures", pdb_key[len(sfam_id):])
            # if not os.path.isdir(os.path.dirname(pdb_file)):
            #     os.makedirs(os.path.dirname(pdb_file))
            # feature_store.read_input_file(pdb_key, pdb_file)
            print(feature_file,feature_file[:-9])
            pdb, chain, sdi, domNo = os.path.basename(feature_file[:-9]).split("_")
            sdi, domNo = sdi[3:], domNo[1:]

            key = os.path.join(sfam_id, pdb[1:3].lower(), os.path.basename(feature_file[:-9]))

            pdb_file = os.path.join(work_dir, "structures", sfam_id, pdb[1:3].lower(),
                "{}.pdb".format(os.path.basename(feature_file[:-9])))

            if not os.path.isfile(pdb_file):
                print("Skipping "+os.path.basename(feature_file[:-8]))

            data.append([key, pdb, chain, sdi, domNo, pdb_file, feature_file])

    data_key = os.path.join(sfam_id, "autoencoder_keys.csv")
    data_file = os.path.join(work_dir, "autoencoder_keys.csv" )
    df = pd.DataFrame(data, columns=["key", "pdb", "chain", "sdi", "domNo",
        "structure", "features"])
    df.to_csv(data_file)
    #normalize(job, df["features"], work_dir=work_dir)
    #sae_store.write_output_file(data_file, data_key)
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
