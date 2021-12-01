import os
import re
from collections import defaultdict

import pandas as pd
from Bio import SearchIO

from molmimic.parsers.container import Container

class HMMER(Container):
    IMAGE = 'docker://edraizen/hmmer:latest'
    LOCAL = ["hmmsearch"]
    RETURN_FILES=True
    ARG_START="-"
    PARAMETERS = [
        (":o", "path:out"),
        (":A", "path:out"),
        (":tblout", "path:out", ["--tblout", "{}"]),
        (":domtblout", "path:out", ["--domtblout", "{}"]),
        (":acc", "store_true", ["--acc"]),
        (":noali", "store_true", ["--noali"]),
        (":notextw", "store_true", ["--notextw"]),
        (":textw", "str", ["--textw", "{}"]),
        (":E", "str"),
        (":T", "str"),
        (":domE", "str", ["--domE", "{}"]),
        (":domT", "str", ["--incE", "{}"]),
        (":incE", "str", ["--incE", "{}"]),
        (":incT", "str", ["--incT", "{}"]),
        (":incdomE", "str", ["--incdomE", "{}"]),
        (":incdomT", "str", ["--incdomT", "{}"]),
        (":incdomT", "str", ["--incdomT", "{}"]),
        (":cut_ga", "store_true", ["--cut_ga"]),
        (":cut_nc", "store_true", ["--cut_nc"]),
        (":cut_tc", "store_true", ["--cut_tc"]),
        (":max", "store_true", ["--max"]),
        (":F1", "str", ["--F1", "{}"]),
        (":F2", "str", ["--F2", "{}"]),
        (":F3", "str", ["--F3", "{}"]),
        (":nobias", "store_true", ["--nobias"]),
        (":nonull2", "store_true", ["--nonull2"]),
        (":Z", "str", ["--Z", "{}"]),
        (":domZ", "str", ["--domZ", "{}"]),
        (":seed", "str", ["--seed", "{}"]),
        (":tformat", "str", ["--tformat", "{}"]),
        (":cpu:4", "str", ["--cpu", "{}"]),
        ("hmmfile", "path:in", ["{}"]),
        ("seqdb", "path:in", ["{}"])
    ]

    def build_combined_hmms(self, *alns, name="combined", press=True):
        """Build HMMs from alignments,
        and outputing them to
        static/browse/hmms/
        to individual dirs as well as combining to pne file combined_hmm_file
        """
        combined_hmm_file = os.path.join(self.work_dir, f"{name}.hmm")
        with open(combined_hmm_file, "w") as combined_hmm:
            for aln in alns:
                print(aln)
                hmm_file = self.build(aln)
                with open(hmm_file) as hmm:
                    print(hmm.read().rstrip(), file=combined_hmm)
        return self.press(combined_hmm_file)

    def build(self, msafile, name=None, out_file_name=None):
        if name is None and out_file_name is None:
            name = os.path.splitext(os.path.basename(msafile))[0]
            out_file_name = os.path.join(os.path.dirname(msafile), f"{name}.hmm")
        elif out_file_name is not None:
            name = os.path.splitext(os.path.basename(msafile))[0]

        build_args = [
            (":n", "str"),
            ("hmmfile_out", "path:out", ["{}"]),
            ("msafile", "path:in", ["{}"])
        ]

        with self.custom_entrypoint("hmmbuild", build_args):
            self(msafile=msafile, hmmfile_out=out_file_name, n=name)

        return out_file_name

    def press(self, combined_hmm):
        """Press the HMMs into a single HMM file, overwriting if present"""
        with self.custom_entrypoint("hmmpress", [("f", "path:in")]):
            self(combined_hmm)
        return combined_hmm

    def search(self, *args, notextw=True, **kwds):
        kwds["notextw"] = notextw
        return self(*args, **kwds)

    def extract_ids(self, sequences):
        with self.custom_entrypoint("esl-sfetch", [("index", "path:in", ["--index", "{}"])]):
            return self(index=sequences)

    def extract_full_sequences(self, sequences, ids_file, output_file):
        with self.custom_entrypoint("esl-sfetch", [("o", "path:out"), ("f", "path:in"), ("ids_file", "path:in", ["{}"])]):
            return self(o=self.full_length_seqs_file, f=sequences, ids_file=self.ids_file)

    @staticmethod
    def parse_hmmsearch_results(results_file, bitscore=True, id_parser=None):
        if id_parser is not None:
            id_parser = re.compile(id_parser)

        results = defaultdict(lambda: defaultdict(float))

        for hmm_query in SearchIO.parse(results_file, "hmmer3-text"):
            hmm_model = os.path.basename(hmm_query.id)
            if hmm_model.endswith(".hmm"):
                hmm_model = hmm_model[:-4]

            for hit in hmm_query:
                #Below we are fetching a list of headers if there are multiple headers for identical sequences
                #Technically HMMER might put the second and on gis in description column.
                #The format should be strictly the genbank format: gi|343434|fafsf gdgfdg gi|65656|534535 fdafaf
                header = "{}{}".format(hit.id, hit.description)
                if id_parser is not None:
                    m = id_parser.search(header)
                    if m:
                        header = m.groups()[0]

                if len(hit)==0:
                    best_score = 0.
                else:
                    if bitscore:
                        best_hsp = max(hit, key=lambda hsp: hsp.bitscore)
                    else:
                        best_hsp = min(hit, key=lambda hsp: hsp.evalue)
                    best_score = best_hsp.bitscore if bitscore else best_hsp.evalue
                    #if best_hsp.bitscore < hit.bitscore:
                        #print(header, "is lower than hit", best_hsp.bitscore, hit.bitscore)
                results[header][hmm_model] = best_score

        return pd.DataFrame.from_dict(results, orient='index')
