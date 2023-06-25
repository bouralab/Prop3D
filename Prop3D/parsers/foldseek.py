import os, sys
import re
import uuid
import subprocess
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

from toil.realtimeLogger import RealtimeLogger

from Prop3D.parsers.container import Container

class FoldSeek(Container):
    IMAGE = 'docker://edraizen/foldseek:latest'
    LOCAL = ["foldseek"]
    ARG_START = "--"
    DETACH=True

class EasyCluster(FoldSeek):
    PARAMETERS = ["easy-cluster",
        ("queryDB", "path:in", ["{}"]),
        ("prefix", "path:out:ignore", ["{}"]),
        ("tmpDir", "path:out:ignore", ["{}"]),
        (":e", str, ["-e", "{}"]), #Sensitivity	List matches below this E-value (range 0.0-inf, default: 0.001); increasing it reports more distant structures
        (":alignment_type", str, ["--alignment-type", "{}"]), #Alignment	0: 3Di Gotoh-Smith-Waterman (local, not recommended), 1: TMalign (global, slow), 2: 3Di+AA Gotoh-Smith-Waterman (local, default)
        (":c", "str", ["-c", "{}"]), #Alignment	List matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0); higher coverage = more global alignment
        (":cov_mode", str, ["--cov-mode", "{}"]), #Alignment	0: coverage of query and target, 1: coverage of target, 2: coverage of query
        (":min_seq_id", str, ["--min-seq-id", "{}"]), 
        (":threads:16", "str", ["--threads", "{}"])
    ]

    def cluster(self, queryDB, prefix=None, tmp_dir=None, e=0.001, alignment_type=2, c=0.9, min_seq_id=0.0, parse=True, **kwds):
        if tmp_dir is None or not isinstance(tmp_dir, str) or \
            not (os.path.isdir(tmp_dir) and str(os.path.abspath(tmp_dir)).startswith(
            str(os.path.abspath(self.work_dir)))):
                tmp_dir = os.path.join(self.work_dir, "tmp4", f"tmp-{uuid.uuid4()}")
                os.makedirs(tmp_dir)

        if prefix is None:
            prefix = f"{Path(queryDB).stem}-db"
        
        self(queryDB, prefix, tmp_dir, e=e, alignment_type=alignment_type, c=c, min_seq_id=min_seq_id, **kwds)

        if parse:
            clusters = pd.read_csv(f"{prefix}_clusters.tsv", sep="\t", names=["cluster_representative", "structure"])
            return clusters
        
        return prefix

    