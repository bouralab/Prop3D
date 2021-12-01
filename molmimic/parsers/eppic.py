from __future__ import print_function
import os
from shutil import copyfileobj
import json
from collections import defaultdict
from functools import partial

import requests
import numpy as np
import pandas as pd

try:
    from pandas.lib import infer_dtype
except ImportError:
    from pandas.api.types import infer_dtype

from Bio import AlignIO
from Bio.Emboss.Applications import NeedleCommandline

from molmimic.util.iostore import IOStore
from molmimic.util import natural_keys, reset_ip
from molmimic.parsers.json import JSONApi
from molmimic.parsers.pdbe import PDBEApi
from molmimic.parsers.container import Container
from molmimic.util.pdb import s3_download_pdb

from toil.realtimeLogger import RealtimeLogger

from joblib import Parallel, delayed

def stringify(df):
    """Make sure all columns with pdb in name are strings (e.g. pdbCode, ...PdbRes...)"""
    for col in (c for c in df.columns if "pdb" in c.lower() or "type" in c.lower()):
        df[col] = df[col].astype(str)
    return df

def resnum_to_biopdb(resnum):
    try:
        resseq_parts = natural_keys(resseq, use_int=True)
        resseq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
        return resseq
    except Exception:
        return np.nan

class EPPICApi(JSONApi):
    def __init__(self, pdb, eppic_store, pdbe_store, use_representative_chains=True,
      work_dir=None, download=True, clean=True, max_attempts=2):
        self.pdbe_api = PDBEApi(pdbe_store, work_dir=work_dir, download=download,
            max_attempts=max_attempts)
        super(EPPICApi, self).__init__("http://www.eppic-web.org/rest/api/v3/job/",
            eppic_store, work_dir=work_dir, download=download, clean=clean,
            max_attempts=max_attempts)

        self.pdb = pdb.lower()
        self.sequences = self.get_sequences()

        RealtimeLogger.info("Running sequences {}".format(self.sequences))
        # self.chains = {chain:rep.repChain for rep in self.sequences.itertuples() \
        #     for chain in rep.memberChains.split(",")}

        self.chains = {}
        for rep in self.sequences.itertuples():
            RealtimeLogger.info("Running rep {} {}".format(pdb, rep))
            for chain in rep.memberChains.split(","):
                RealtimeLogger.info("Running chain {}".format(chain))
                self.chains[chain] = rep.repChain

        self.use_representative_chains = use_representative_chains

        self._interface_chains = {}
        self._chain_residues = {} if not use_representative_chains else None

    def parse(self, file_path, key):
        result = super().parse(file_path, key)
        if len(result) == 1 and isinstance(result[0], dict) and len(result[0]) == 1 and "uid" in result[0]:
            #Reset to download again
            RealtimeLogger.info("LOAD {} {} {} {}".format(type(result), len(result), type(result[0]), len(result[0])))

            raise ValueError("incomplete data")
        return result

    def _get_chain_residues(self, chain):
        if not self.use_representative_chains:
            if chain not in self._chain_residues:
                self._chain_residues[chain] = self.pdbe_api.get_pdb_residues(
                    self.pdb, chain)
            return self._chain_residues[chain]

    def check_line(self, key, line, attempts):
        RealtimeLogger.info("EPPIC Failed {}".format(line))
        rerun = False
        if "<html>" in line or "<head>" in line or "HTTP 404 Not Found" in line:
            RealtimeLogger.info("Is restarting")
            self.store.remove_file(key)
            rerun = attempts > 0

        if "Too many submissions" in line:
            RealtimeLogger.info("Is changing IP")
            try:
                reset_ip()
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                pass

        try:
            RealtimeLogger.info("Start {} line {}".format(key, line))
            result = json.loads(line)
            RealtimeLogger.info("Result {} {}".format(key, result))
            if len(result) == 1 and isinstance(result[0], dict) and "uid" in result[0]:
                rerun = True
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            import traceback
            tb = traceback.format_exc()
            if "JSON" in tb:
                RealtimeLogger.info("Failed {} Not found JSON error".format(key))
            RealtimeLogger.info("Failed {} {}".format(key, tb))
            rerun = True

        return rerun

    def get_interfaces(self, bio=False):
        result = self.get(("interfaces", self.pdb))

        interfaces = pd.DataFrame(result)

        if interfaces.empty:
            return None

        interfaces = interfaces.assign(interfaceType=
            interfaces["interfaceScores"].apply(lambda x:
                [score["callName"] for score in x["interfaceScore"] \
                    if score["method"] == "eppic"][0]))

        if bio:
            interfaces = interfaces[interfaces["interfaceType"]=="bio"]

        #No residue information so no need to trasnform pdbResNums

        chains = interfaces.set_index("interfaceId")[["chain1", "chain2"]].T.to_dict("list")
        self._interface_chains.update(chains)

        return stringify(interfaces)

    def get_interface_chains(self, interfaceId):
        if interfaceId not in self._interface_chains:
            self.get_interfaces(interfaceId)
        return self._interface_chains[interfaceId]

    def get_interface_residues(self, interfaceId):
        result = self.get(("interfaceResidues", self.pdb, interfaceId))
        result = stringify(pd.DataFrame(result))

        if not self.use_representative_chains:
            result = self.expand_chains_from_interface(result, interfaceId)
            # chain1, chain2 = self.get_interface_chains(interfaceId)
            # chain1 = pd.merge(result[result["side"]==0], self.chain_residues[chain1],
            #     how="left", on="residueNum", suffixes=["_eppic", "_pdbe"])
            # chain2 = pd.merge(result[result["side"]==1], self.chain_residues[chain2],
            #     how="left", on="residueNum", suffixes=["_eppic", "_pdbe"])
            # result = pd.concat((chain1, chain2), axis=0)
            # result = result.drop(columns=["pdbResidueNum_eppic"]).rename( \
            #     columns={"pdbResidueNum_pdbe":"pdbResidueNum"})

        return result

    def get_contacts(self, interfaceId):
        result = self.get(("contacts", self.pdb, interfaceId))

        #No pdb residue information so no need to trasnform pdbResNums

        return stringify(pd.DataFrame(result))

    def get_sequences(self):
        if hasattr(self, "sequences") and self.sequences is not None:
            return self.sequences

        result = self.get(("sequences", self.pdb))
        return stringify(pd.DataFrame(result))

    def expand_chain(self, df, chain_lett):
        chain = pd.merge(df, self._get_chain_residues(chain_lett), how="left",
            on="residueNumber", suffixes=["_eppic", "_pdbe"])
        chain = chain.drop(columns=["pdbResidueNumber_eppic"]).rename(
            columns={"pdbResidueNumber_pdbe":"pdbResidueNumber"})
        return chain

    def expand_chains_from_interface(self, df, interfaceId):
        chain1_lett, chain2_lett = self.get_interface_chains(interfaceId)

        if self.chains[chain1_lett] == chain1_lett:
            chain1 = df[df["side"]==0]
        else:
            chain1 = self.expand_chain(df[df["side"]==0], chain1_lett)

        if self.chains[chain2_lett] == chain2_lett:
            chain2 = df[df["side"]==1]
        else:
            chain2 = self.expand_chain(df[df["side"]==1], chain2_lett)

        chains = pd.concat((chain1, chain2), axis=0)
        return chains

    def expand_chains_from_residues(self, df, chainLett):
        if self.chains[chainLett] == chainLett:
            return df
        else:
            return self.expand_chain(df, chainLett)

    def get_residue_info(self, chain=None):
        seq = self.get_sequences()
        seq = seq.drop(seq.columns.difference(['memberChains','residueInfos']), 1)

        if seq.empty:
            return None

        if chain is not None:
            seq = seq[seq["memberChains"].str.contains(chain)]
            seq["memberChains"] = chain

        result = None
        for seq_chains in seq.itertuples():
            repChain = seq_chains.memberChains.split(",")[0]
            seq_chain_clustered = pd.DataFrame(seq_chains.residueInfos)
            seq_chain_clustered = seq_chain_clustered.dropna(subset=["pdbResidueNumber"])

            for memberChain in seq_chains.memberChains.split(","):
                if self.use_representative_chains:
                    seq_chain = seq_chain_clustered.assign(chain=memberChain)
                else:
                    #Do not use representative and get from PDBE
                    seq_chain = self.expand_chains_from_residues(seq_chain_clustered, memberChain)
                    seq_chain = seq_chain.assign(chain=memberChain)

                if result is None:
                    result = seq_chain
                else:
                    result = pd.concat((result, seq_chain), axis=0)

        return stringify(result)

    def get_entropy_scores(self, chain, maxEntropy=1):
        residueInfo = self.get_residue_info(chain)
        entropy_scores = defaultdict(lambda: maxEntropy, zip(
            residueInfo.pdbResidueNumber, residueInfo.entropyScore))
        return entropy_scores

    # def process_obsolete(self):
    #     parser = PDB.PDBParser()
    #     pdbl = PDBList()
    #
    #     old_pdb_file = pdbl.retrieve_pdb_file(self.pdb.upper(), file_format="pdb", obsolete=True)
    #     old_structure = parser.get_structure(old_pdb_file, "")
    #     old_seq = {chain.id:"".join(PDB.Polypeptide.three_to_one(r.resname) for r in \
    #         chain.get_residues() if r.get_id()[0] == " ") for chain in old_structure.get_chains()}
    #
    #     if not hasattr(self, "obsolete"):
    #         self.obsolete = pd.read_csv("https://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat",
    #             skiprows=2, header=None, delim_whitespace=True, names=
    #                 ["0", "date", "old", "new1", "new2", "new3"], engine="python")
    #
    #     scores = {k:None for k in old_seq.keys()}
    #
    #     superceded_pdbs = self.obsolete[self.obsolete["old"]==self.pdb.upper()]
    #     if len(superceded_pdbs) == 0:
    #         return None
    #
    #     superceded_pdbs = superceded_pdbs.iloc[0][["new1", "new2", "new3"]]
    #     for new_pdb in superceded_pdbs:
    #         if new_pdb is None or new_pdb == float("nan"): continue
    #
    #         new_pdb_file = pdbl.retrieve_pdb_file(new_pdb, file_format="pdb", obsolete=True)
    #         new_structure = parser.get_structure(new_pdb_file, "")
    #         new_seq = {chain.id:"".join(PDB.Polypeptide.three_to_one(r.resname) for r in \
    #             chain.get_residues() if r.get_id()[0] == " ") for chain in new_structure.get_chains()}
    #
    #         for chain1, seq1 in old_seq.items():
    #             for chain2, seq2 in new_seq.items():
    #                 alignments = pairwise2.align.globalms(seq1, seq2, 5, -4, -10, -1,
    #                                   penalize_end_gaps=False, one_alignment_only=True)
    #
    #                 if scores[chain1] is None or alignments[0].score > scores[chain1][0].score:
    #                     scores[chain1] = (alignments[0], new_pdb, chain2)
    #
    #     eppic_for_pdb = {}
    #     for _, new_pdb, _ in scores:
    #         if new_pdb not in eppic_for_pdb:
    #             eppic_for_pdb[new_pdb] = EPPICApi(new_pdb, self.eppic_store, self.pdbe_store,
    #                 use_representative_chains=self.use_representative_chains,
    #                 work_dir=self.work_dir, download=self.download, clean=self.clean,
    #                 max_attempts=self.max_attempts):
    #
    #


class EPPICLocal(Container):
    IMAGE = "docker://edraizen/eppic"
    PARAMETERS = [
        ("input_file", "path:in", "i"),
        (":s", "store_true"),
        (":p", "store_true"),
        (":l", "store_true"),
    ]
    ARG_START = "-"

    def __init__(self, pdb, eppic_local_store=None, job=None, return_files=False,
      force_local=False, fallback_local=False, work_dir=None):
        super().__init__(job=job, return_files=return_files, force_local=force_local,
            fallback_local=fallback_local, work_dir=work_dir)

        if eppic_local_store is None:
            eppic_local_store = data_stores.eppic_local_store

        self.pdb = pdb

        if not os.path.isfile(input) and len(input)==4:
            if eppic_local_store.exists(self.pdb.lower()):
                self.eppic_files = [os.path.basename(f) for f in \
                    eppic_local_store.list_input_directory(self.pdb.lower())]
            else:
                self.eppic_files = []
        else:
            self.eppic_files = None

        self._interface_chains = {}
        self._chain_residues = {} if not use_representative_chains else None

    def run(self, s=True, p=True, l=True):
        if not os.path.isfile(self.pdb) and len(pdb)==4:
            if len(self.eppic_files) == 0:
                pdb_file = s3_download_pdb(self.pdb)
                self(input_file=pdb_file, s=s, p=p, l=l)
                self.eppic_files = [os.path.basename(f) for f in \
                    glob.glob("{}.*".format(os.path.splitext(self.pdb)[0]))]
                for f in self.eppic_files:
                    eppic_local_store.write_output_file(f, "{}/{}".format(self.pdb.lower(), f))
            else:
                for f in self.eppic_files:
                    eppic_local_store.read_input_file(f, os.path.join(self.work_dir, f))
        else:
            #If not pdb id, just run and dont store/check in s3
            self(input_file=self.pdb, s=s, p=p, l=l)

    def get_interfaces(self, bio=True):
        score_file = "{}.scores".format(self.pdb)
        assert os.path.isfile(score_file)
        interfaces = pd.read_csv(score_file, whitespace_delim=True)

        if interfaces.empty:
            return None

        if bio:
            interfaces = interfaces[interfaces["call"]=="bio"]

        #No residue information so no need to trasnform pdbResNums

        chains = interfaces["chains"].str.split("+", expabd=True).T.to_dict("list")
        self._interface_chains.update(chains)

        return stringify(interfaces)



    def get_entropy_scores(self, maxEntropy=1):
        self.run()
        entropy_file = os.path.join(self.work_dir, "{}.entropies".format(
            os.path.splitext(os.path.basename(pdb_file))))
        assert os.path.isfile(entropy_file), "Failed running EPPIC"
        residueInfo = pd.read_csv(entropy_file, delimeter="\t", comment="#",
            names=["residueNumber", "pdbResidueNumber", "uniProtNumber",
                "residueType", "entropyScore"])
        entropy_scores = defaultdict(lambda: maxEntropy, zip(
            residueInfo.pdbResidueNumber, residueInfo.entropyScore))
        return entropy_scores
