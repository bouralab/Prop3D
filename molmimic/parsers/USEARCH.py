import os, sys
import glob

import pandas as pd

from toil.realtimeLogger import RealtimeLogger

from molmimic.parsers.container import Container

class USEARCH(Container):
    IMAGE = 'docker://edraizen/usearch:latest'
    LOCAL = ["usearch"]
    RETURN_FILES = True

    @staticmethod
    def parse_uc_file(uclust_file):
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

        return uclust

class ClusterFast(USEARCH):
    ARG_START = "-"
    PARAMETERS = [#"-cluster_fast",
        ("fasta_file", "path:in", ["-cluster_fast", "{}"]),
        (":id", "str"),
        (":centroids", "path:out"),
        (":uc", "path:out"),
        (":consout", "path:out"),
        (":msaout", "path:out:ignore"),
        (":sort", "str"),
        (":threads", "str")
    ]

    def __call__(self, *args, **kwds):
        output = super().__call__(*args, **kwds)

        RealtimeLogger.info("FINISHED UCLUST {}".format(output))
        if "msaout" in kwds and isinstance(output, dict):
            name = output["msaout"]
            RealtimeLogger.info("Similar files: {}".format([f for f in os.listdir(self.work_dir) if "all_binding_sites" in f][:5]))
            output["msaout"] = list(glob.glob(os.path.join(self.work_dir, output["msaout"]+"*")))
            RealtimeLogger.info("Number of MSAOUT: {} {}".format(name, len(output["msaout"])))
        return output

    def cluster(self, fasta_file, id=0.9, centroids=True, uc=True, msa=True, cons=False, parse_uc=True, sort=None, threads=None):
        prefix = os.path.join(self.work_dir, os.path.splitext(os.path.basename(fasta_file))[0])

        RealtimeLogger.info("CLUSTER: msa={}; centroids={};, uc={}; cons={}".format(msa, centroids, uc, cons))

        parameters = {"fasta_file":fasta_file, "id":id}
        if centroids:
            parameters["centroids"] = os.path.join(self.work_dir, "{}_clusterfast.centroids".format(prefix))

        if uc:
            parameters["uc"] = os.path.join(self.work_dir, "{}_clusterfast.uc".format(prefix))

        if msa:
            parameters["msaout"] = os.path.join(self.work_dir, "{}_clusterfast_msa.fasta.".format(prefix))

        if cons:
            parameters["consout"] = os.path.join(self.work_dir, "{}_clusterfast_cons.fasta".format(prefix))

        if isinstance(threads, int):
            parameters["threads"] = str(threads)

        RealtimeLogger.info("RUNNING UCLUST with {}".format(parameters))

        cluster_output = self(**parameters)

        if parse_uc and "uc" in cluster_output:
            cluster_output["uc"] = self.parse_uc_file(cluster_output["uc"])

        return cluster_output

class AllPairsLocal(USEARCH):
    ARG_START = "-"
    PARAMETERS = [#"-cluster_fast",
        ("fasta_file", "path:in", ["-allpairs_local", "{}"]),
        (":id", "str"),
        (":centroids", "path:out"),
        (":uc", "path:out"),
        (":consout", "path:out"),
        (":alnout", "path:out:ignore"),
        (":sort", "str"),
        (":threads", "str"),
        (":acceptall", "store_true")
    ]

class AllPairsGlobal(USEARCH):
    ARG_START = "-"
    PARAMETERS = [#"-cluster_fast",
        ("fasta_file", "path:in", ["-allpairs_global", "{}"]),
        (":id", "str"),
        (":centroids", "path:out"),
        (":uc", "path:out"),
        (":consout", "path:out"),
        (":alnout", "path:out:ignore"),
        (":sort", "str"),
        (":threads", "str"),
        (":acceptall", "store_true")
    ]
