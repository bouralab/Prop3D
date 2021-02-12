import sys
import os
import argparse
import time
from itertools import groupby
import glob

import pandas as pd

from molmimic.util.iostore import IOStore
from molmimic.util.hdf import get_file, filter_hdf, filter_hdf_chunks
from molmimic.util.toil import map_job
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

def cluster_by_structure(job, sfam, structures, work_dir=None):
    from molmimic.parsers.MaxCluster import run_maxcluster, get_centroid
    from molmimic.generate_data.util import PDB_TOOLS, SubprocessChain

    cluster_store = IOStore.get("aws:us-east-1:molmimic-clusters")

    if work_dir is None:
        work_dir = os.getcwd()

    cluster_dir = os.path.join(work_dir, "{}_cluster".format(int(sfam)))
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    file_list = os.path.join(cluster_dir, "{}_file_list.txt".format(int(sfam)))

    with open(file_list, "w") as f:
        for pdb_file in structures:
            if not pdb_file: continue
            cmds = [
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-A", pdb_file],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_reres.py")],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_reatom.py")],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
            ]

            pdb_to_cluster = os.path.join(cluster_dir, os.path.basename(pdb_file))
            with open(pdb_to_cluster, "w") as pdb:
                SubprocessChain(cmds, pdb)

            print >> f, os.path.basename(pdb_to_cluster)

    log_file = os.path.join(cluster_dir, "{}.max_cluster_logs".format(int(sfam)))
    distance_file = os.path.join(cluster_dir, "{}_distances.txt".format(int(sfam)))

    key_base = "{}/structures".format(int(mol_sfam), int(sfam))

    logs = run_maxcluster("rmsd", file_list=file_list, log=True, R=distance_file,
        work_dir=cluster_dir, C=1, P=10, job=job)

    RealtimeLogger.info("CLUSTER_DIR: {}".format(os.listdir(cluster_dir)))
    RealtimeLogger.info("LOG FILE: {}".format(logs))

    with open(logs) as f:
        RealtimeLogger.info("LOG IS {}".format(next(f)))

    cluster_store.write_output_file(file_list, "{}/{}".format(key_base, os.path.basename(file_list)))
    cluster_store.write_output_file(logs, "{}/{}".format(key_base, os.path.basename(log_file)))
    cluster_store.write_output_file(distance_file, "{}/{}".format(key_base, os.path.basename(distance_file)))

    return get_clusters(logs)

def cluster_by_sequence(job, sfam, structure_data, work_dir=None, id=0.90, preemptable=True):
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir() if job is not None else os.getcwd()
    out_store = IOStore.get("aws:us-east-1:molmimic-clusters")

    sdoms_file = copy_pdb_h5(job, pdbFileStoreID)

    sdoms = pd.read_hdf(str(sdoms_file), "merged")
    sdoms = sdoms[sdoms["sfam_id"]==sfam_id]
    sdoms = sdoms[["pdbId", "chnLett", "sdi", "domNo"]].drop_duplicates().dropna()

    #Save all domains to fasta
    domain_fasta = os.path.join(work_dir, "{}.fasta".format(int(sfam_id)))
    domain_ids = {}
    with open(domain_fasta, "w") as fasta:
        for pdb_fname in structure_data["structure"]:
            try:
                seq = subprocess.check_output([sys.executable, os.path.join(PDB_TOOLS,
                    "pdb_toseq.py"), pdb_fname])
                fasta.write(">{}\n{}\n".format(pdb_fname, "\n".join(seq.splitlines()[1:])))
                domain_ids[pdb_fname] = jobStoreID
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as e:
                RealtimeLogger.info("Error getting fasta for : {} {}".format(sfam_id, fname))
                pass

    #Order domains by resolution so the ones with the highest resolutions are centroids
    resolutions = pd.read_hdf(str(sdoms_file), "resolu")

    try:
        pdbs, ids, sequences = list(zip(*[(s.id.split("_", 1)[0].upper(), s.id, str(s.seq)) \
            for s in SeqIO.parse(domain_fasta, "fasta")]))
    except ValueError:
        RealtimeLogger.info("Unable to cluster {}".format(sfam_id))
        return

    domains = pd.DataFrame({"pdbId":pdbs, "domainId":ids, "sequence":sequences})
    domains = pd.merge(domains, resolutions, how="left", on="pdbId")
    xray = domains[domains["resolution"] >= 0.].sort_values("resolution")
    nmr = domains[domains["resolution"] < 0.]

    with open(domain_fasta, "w") as f:
        for row in it.chain(xray.itertuples(index=False), nmr.itertuples(index=False)):
            print(">{} [resolution={}]\n{}".format(row.domainId, row.resolution, row.sequence))

    sfam_key = "{0}/{0}.fasta".format(int(sfam_id))
    out_store.write_output_file(domain_fasta, sfam_key)

    clusters_file, uclust_file = run_usearch(["-cluster_fast",
        "{i}"+domain_fasta, "-id", str(id),
        "-centroids", "{out}"+"{}_clusters.uc".format(int(sfam_id)),
        "-uc", "{o}"+"{}_clusters.uc".format(int(sfam_id))])

    #Convert uclust to h5
    uclust = pd.read_table(str(uclust_file), comment="#", header=None, names=[
        "record_type",
        "cluster",
        "length",
        "pctId",
        "strand",
        "unk1",
        "unk2",
        "alignment",
        "label_query",
        "label_target"
    ])
    del uclust["unk1"]
    del uclust["unk2"]
    hdf_base = "{}_clusters.h5".format(int(sfam_id))
    hdf_file = os.path.join(work_dir, hdf_base)
    uclust.to_hdf(str(hdf_file), "table", complevel=9, complib="bzip2")
    out_store.write_output_file(hdf_file, "{}/{}".format(int(sfam_id), hdf_base))
    os.remove(uclust_file)

    return clustered_pdbs

def compare_clusters(job, seq_clusters, struct_clusters):
    seq_clusters = seq_clusters[["label_query", "cluster"]].rename(columns={
        "label_query":"pdb_file", "cluster":"sequence_cluster"})
    clusters = pd.merge(seq_clusters, struct_clusters, on="pdb_file")


def cluster(job, structures_data):
    seq_job = job.addChildJobFn(cluster_by_sequence, structures_data)
    struct_job = job.addChildJobFn(cluster_by_structre, structures_data["structures"])
    job.addFollowOnJobFn(compare_clusters, seq_job.rv(), struct_job.rv())

def start_sae(job, sfam_id, work_dir=None, sae=True):
    if work_dir is None and job is not None:
        work_dir = job.fileStore.getLocalTempDir()
    elif work_dir is None:
        work_dir = os.getcwd()
    # pdb_store = IOStore.get("aws:us-east-1:molmimic-full-structures")
    # feature_store = IOStore.get("aws:us-east-1:molmimic-structure-features")
    #sae_store = IOStore.get("aws:us-east-1:molmimic-structure-features")

    #if not sae:
    #    ibis_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")
    #
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

    int_file = os.path.join(work_dir, "interfaces", "{}.observed_interactome".format(sfam_id))
    ibis_df = pd.read_hdf(int_file, "table")

    for feature_file in glob.glob(os.path.join(work_dir, "features", sfam_id, "*", "*_atom.npy")):
        if True:
            print(feature_file)
            #data_key = feature_key[:-9]
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
            #print(feature_file,feature_file[:-9])
            pdb, chain, sdi, domNo = os.path.basename(feature_file[:-9]).split("_")
            sdi, domNo = sdi[3:], domNo[1:]

            key = os.path.join(sfam_id, pdb[1:3].lower(), os.path.basename(feature_file[:-9]))

            pdb_file = os.path.join(work_dir, "structures", sfam_id, pdb[1:3].lower(),
                "{}.pdb".format(os.path.basename(feature_file[:-9])))

            if not os.path.isfile(pdb_file):
                print("Skipping "+os.path.basename(feature_file[:-8]))

            datum = [key, pdb, chain, sdi, domNo, pdb_file, feature_file]

            if not sae:
                curr_int = ibis_df[ibis_df["mol_sdi_id"]==int(sdi)]
                res = set([r for res in curr_int["mol_res"] for r in res.split(",")])
                if len(res) == 0:
                    continue

                datum.append(",".join(map(str, sorted(res))))

            data.append(datum)

    data_key = os.path.join(sfam_id, "autoencoder_keys.csv")
    data_file = os.path.join(work_dir, "autoencoder_keys.csv" )
    cols = ["key", "pdb", "chain", "sdi", "domNo", "structure", "features"]
    if not sae:
        cols +=["mol_res"]
    df = pd.DataFrame(data, columns=cols)
    df.to_csv(data_file)
    #normalize(job, df["features"], work_dir=work_dir)
    #sae_store.write_output_file(data_file, data_key)
    del df
    del data

    if not sae:
        return data_file

    def callback(epoch, statsfile, graphs, model_file, epoch_file):
        pass
#         key = os.path.join(sfam_id, epoch)
#         sae_store.write_output_file(statsfile, os.path.join(key, statsfile))
#         sae_store.write_output_file(model_file, os.path.join(key, model_file))
#         sae_store.write_output_file(epoch_file, os.path.join(key, epoch_file))
#         for graph in graphs:
#             sae_store.write_output_file(graph, os.path.join(key, graph))


    train(data_file, model_prefix="{}_sae".format(sfam_id), num_epochs=100,
        use_resnet_unet=True, dropout_depth=True,
        dropout_width=True, dropout_p=0.6, autoencoder=sae,
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
