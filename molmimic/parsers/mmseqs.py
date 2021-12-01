import os, sys
import re
import uuid
import subprocess

import pandas as pd
from joblib import Parallel, delayed

from toil.realtimeLogger import RealtimeLogger

from molmimic.parsers.container import Container

class MMSeqs(Container):
    IMAGE = 'docker://edraizen/mmseqs:latest'
    LOCAL = ["mmseqs2"]
    RETURN_FILES = True
    ARG_START = "--"
    DETACH=True

    def all_vs_all(self, fasta_file, output=None, db=None):
        output = os.path.join(self.work_dir, os.path.splitext(os.path.basename(fasta_file))[0])

        if db is None:
            db = f"{output}.mmseqs_db"
            CreateDB(work_dir=self.work_dir)(fasta_file, db)


        allvsallpref = self.fake_pref(db, db, output+".all_vs_all_pref")

        new_params = ["align",
            ("qdb", "path:in", ["{}"]),
            ("tdb", "path:in", ["{}"]),
            ("clu", "path:in", ["{}"]),
            ("aln_out", "path:out:ignore", ["{}"]),
            ("a", "store_true", ["-a"]),
            (":min_covered", "str", ["-c", "{}"]),
            (":max_evalue", "str", ["-e", "{}"]),
            (":min_seq_id", "str", "min-seq-id"),
            (":add_self_matches", "store_true", ["--add-self-matches"]),
            (":cov_mode", "str", "cov-mode"),
            (":seq_id_mode", "str", "seq-id-mode")
            #(":min_covered", "str", ["-c", "{}"])
            #(":sensitivity", "str", ["-s", "{}"]),
            #(":max_seqs:2147483647", "str", "max-seqs")
        ]

        aln_result_file = f"{output}.mmseqs_aln_results"
        with self.custom_parameters(new_params):
            self(qdb=db, tdb=db, clu=allvsallpref, aln_out=aln_result_file, max_evalue=float("inf"),
                min_seq_id=0.0, seq_id_mode=1, cov_mode=2, min_covered=0.0, a=True, add_self_matches=True)

        aln_result_tsv_file = f"{output}.mmseqs_aln_results.tsv"
        CreateTSV(work_dir=self.work_dir)(db, db, aln_result_file, aln_result_tsv_file)

        return aln_result_file, aln_result_tsv_file

    @classmethod
    def fake_pref(cls, qdb, tdb, res):
        """Create a fake prefiltering result for all-vs-all-alignments"""

        # create link to data file which contains a list of all targets that should be aligned
        #subprocess.call("ln -s \"{TDB}.index\" \"{RES}\"".format(TDB=tdb, RES=res), shell=True)
        subprocess.call("cp \"{TDB}.index\" \"{RES}\"".format(TDB=tdb, RES=res), shell=True)

        # create new index repeatedly pointing to same entry
        index_size = subprocess.check_output("wc -c \"{TDB}.index\"".format(TDB=tdb), shell=True)
        index_size = index_size.decode("utf-8").split()[0]
        print(index_size)
        with open(f"{res}.index", "w") as f:
            subprocess.call("awk -v size={INDEX_SIZE} '{{ print $1\"\t0\t\"size; }}' \"{QDB}.index\"".format(INDEX_SIZE=index_size, QDB=qdb), stdout=f, shell=True)

        # create dbtype (7)
        subprocess.call("awk 'BEGIN {{ printf(\"%c%c%c%c\",7,0,0,0); exit; }}' > \"{RES}.dbtype\"".format(RES=res), shell=True)

        return res

class CreateDB(MMSeqs):
    PARAMETERS = ["createdb",
        ("sequence_file", "path:in", ["{}"]),
        ("output_db", "path:out", ["{}"])
    ]

class CreateTSV(MMSeqs):
    PARAMETERS = ["createtsv",
        ("queryDB", "path:in", ["{}"]),
        ("targetDB", "path:in", ["{}"]),
        ("resultDB", "path:in", ["{}"]),
        ("tsvFile", "path:out", ["{}"]),
        (":threads", "str")
    ]

class EasyCluster(MMSeqs):
    PARAMETERS = ["easy-cluster",
        ("fastaFile", "path:in", ["{}"]),
        ("clusterPrefix", "path:out:ignore", ["{}"]),
        ("tmp", "path:in", ["{}"]),
        (":sensitivity", "str", ["-s", "{}"]),
        (":max_seqs:2147483647", "str", "max-seqs"),
        (":min_ungapped_score", "str", "min-ungapped-score"),
        (":add_self_matches", "store_true", ["--add-self-matches"]),

        (":min_covered", "str", ["-c", "{}"]),
        (":cov_mode", "str", "cov-mode"),
        ("a", "store_true", ["-a"]),
        (":alignment_mode", "str", "alignment-mode"),
        (":alignment_output_mode", "str", "alignment-output-mode"),
        (":max_evalue", "str", ["-e", "{}"]),
        (":min_seq_id", "str", "min-seq-id"),

        (":cluster_mode", "str", "cluster-mode"),
        (":cluster_reassign", "store_true", ["--cluster-reassign"]),


        (":threads", "str"),

        (":kmer_per_seq", "str", ["--kmer-per-seq", "{}"])
    ]

    def cluster(self, fastaFile, tmp_dir=None, aggressive=False, prefix="", **kwds):
        for k in ["fastaFile", "tmp", "clusterPrefix"]:
            try:
                del kwds[k]
            except KeyError:
                pass


        if tmp_dir is None or not isinstance(tmp_dir, str) or \
          not (os.path.isdir(tmp_dir) and str(os.path.abspath(tmp_dir)).startswith(
          str(os.path.abspath(self.work_dir)))):
            tmp_dir = os.path.join(self.work_dir, "tmp4", f"tmp-{uuid.uuid4()}")
            os.makedirs(tmp_dir)

        if prefix is None:
            prefix = os.path.join(self.work_dir, os.path.splitext(os.path.basename(fastaFile))[0])

        aggressive_values = dict(sensitivity=8, min_ungapped_score=15,
            min_covered=0.0, max_evalue=1000, #float("inf"),
            min_seq_id=0.0, cluster_mode=0)

        if aggressive:
            kwds.update(aggressive_values)

        print("run", kwds, tmp_dir, self.work_dir)
        self(fastaFile=fastaFile, clusterPrefix=prefix, tmp=tmp_dir, **kwds)

        return f"{prefix}_cluster.tsv"


    def try_all(self, fastaFile, representative_domains, **kwds):
        aggressive_values = dict(
            sensitivity=[6, 7.5, 8, 9.5, 9, 9.5, 10],
            min_ungapped_score=range(31),
            min_covered=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            max_evalue=[10, 100, 1000, float("inf")],
            min_seq_id=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
            cluster_mode=[0,1,2,3])

        import itertools
        keys, values = zip(*aggressive_values.items())
        permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]
        print("Trying", len(permutations_dicts), "hyperparameter combinations")



        prefix = os.path.join(self.work_dir, os.path.splitext(os.path.basename(fastaFile))[0])
        with open(f"{prefix}.hyperparameter_combos", "w") as f:
            print("#n_clusters\tn_sfams\thyperparameters", file=f)

        os.makedirs(prefix, exist_ok=True)

        tmp_dir = os.path.join(self.work_dir, "tmp4")
        os.makedirs(tmp_dir, exist_ok=True)

        hp_dir = os.path.join(self.work_dir, "hyper_params")
        os.makedirs(hp_dir, exist_ok=True)

        def run(i, fasta, hparams, other):
            if i%500==0:
                print(f"Running hyperparameter set: {i}/{len(permutations_dicts)}")
            id_parser = re.compile(r"cath\|4_3_0\|([a-zA-Z0-9]+)\/")
            hparams.update(other)
            os.makedirs(os.path.join(self.work_dir, f"hyper_params", str(i)), exist_ok=True)
            hparams["prefix"] = os.path.join(f"hyper_params", str(i), os.path.splitext(os.path.basename(fastaFile))[0])
            hparams["tmp_dir"] = os.path.join(f"hyper_params", str(i), f"tmp-{uuid.uuid4()}")
            os.makedirs(hparams["tmp_dir"], exist_ok=True)
            results_file = EasyCluster(work_dir=prefix).cluster(fasta, **hparams)
            results = pd.read_csv(results_file, sep="\t", header=None, names="source target targetID alnScore seqIdentity eVal qStart qEnd qLen tStart tEnd tLen cigar".split())
            results["source"] = results["source"].apply(lambda x: id_parser.search(x).groups()[0] if id_parser.search(x) else None)
            results["target"] = results["source"].apply(lambda x: id_parser.search(x).groups()[0] if id_parser.search(x) else None)
            results = pd.merge(results, representative_domains, left_on="source", right_on="cathDomain")

            n_clusters = len(results[["source", "superfamily"]].drop_duplicates())
            n_sfams = len(results["superfamily"].drop_duplicates())

            with open(f"{prefix}.hyperparameter_combos_{i}.tsv", "w") as f:
                print(f"{n_clusters}\t{n_sfams}\t{hparams}", file=f)
                print(f"{n_clusters}\t{n_sfams}\t{hparams}")

        Parallel(n_jobs=20)(delayed(run)(i, fastaFile, hyperparameters, kwds) for i, hyperparameters in enumerate(permutations_dicts))

class Cluster(MMSeqs):
    RETURN_FILES = False
    PARAMETERS = ["cluster",
        ("database_file", "path:in", ["{}"]),
        ("result_file", "path:out:ignore", ["{}"]),
        ("tmp", "path:in", ["{}"]),
        (":min_covered", "str", ["-c", "{}"]),
        (":max_evalue", "str", ["-e", "{}"]),
        (":min_seq_id", "str", "min-seq-id"),
        (":max_seqs:2147483647", "str", "max-seqs"),
        (":cov_mode", "str", "cov-mode"),
        (":alignment_mode", "str", "alignment-mode"),
        (":alignment_output_mode", "str", "alignment-output-mode"),
        (":cluster_reassign", "store_true", ["--cluster-reassign"]),
        (":cluster_mode", "str", "cluster-mode"),
        (":cov_mode", "str", "cov-mode"),
        (":threads", "str"),
        (":sensitivity", "str", ["-s", "{}"]),
        (":kmer_per_seq", "str", ["--kmer-per-seq", "{}"])
    ]



    def cluster(self, fasta_file, tmp_dir=None, **kwds):
        try:
            del kwds["database_file"]
        except KeyError:
            pass

        try:
            del kwds["tmp"]
        except KeyError:
            pass

        if tmp_dir is None or not isinstance(tmp_dir, str) or \
          not (os.path.isdir(tmp_dir) and str(os.path.abspath(tmp_dir)).startswith(
          str(os.path.abspath(self.work_dir)))):
            tmp_dir = os.path.join(self.work_dir, f"tmp-{uuid.uuid4()}")
            os.makedirs(tmp_dir)

        prefix = os.path.join(self.work_dir, os.path.splitext(os.path.basename(fasta_file))[0])
        db_file = f"{prefix}.mmseqs_db"
        CreateDB(work_dir=self.work_dir)(fasta_file, db_file)

        result_file = f"{prefix}.mmseqs_results"
        #kwds["alignment_output_mode"] = 5
        kwds["database_file"] = db_file
        kwds["result_file"] = result_file
        kwds["tmp"] = tmp_dir
        #kwds["min_covered"] = 1.0
        #kwds["max_evalue"] = 1e10
        print(kwds)
        #self(**kwds)
        return self.all_clust(fasta_file, prefix, db=db_file)

        # result_tsv_file = f"{prefix}.mmseqs_results.tsv"
        # CreateTSV(work_dir=self.work_dir)(db_file, db_file, result_file, result_tsv_file)
        #
        # self.all_vs_all(db_file, prefix)
        #
        # return result_tsv_file

    def all_clust(self, db, output):
        aln_result_file, _ = self.all_vs_all(db, output)

        new_params = ["clust",
            ("sequenceDB", "path:in", ["{}"]),
            ("resultDB", "path:in", ["{}"]),
            ("clusterDB", "path:out:ignore", ["{}"]),
            (":cluster_mode", "str", "cluster-mode"),
            (":max_iterations", "str", "max-iterations"),
            (":similarity_type", "str", "similarity-type"),
            (":threads", "store_true"),
        ]

        cluster_result_file = f"{output}.mmseqs_cluster_results"
        with self.custom_parameters(new_params):
            self(sequenceDB=db, resultDB=aln_result_file, clusterDB=cluster_result_file, cluster_mode=1, similarity_type=1)

        cluster_result_tsv_file = f"{output}.mmseqs_cluster_results.tsv"
        CreateTSV(work_dir=self.work_dir)(db, db, cluster_result_file, cluster_result_tsv_file)

        return cluster_result_tsv_file
