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

from molmimic.util.iostore import IOStore
from molmimic.util import natural_keys, reset_ip
from molmimic.parsers.json import JSONApi
from molmimic.parsers.pdbe import PDBEApi

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
        self.chains = {chain:rep.repChain for rep in self.sequences.itertuples() \
            for chain in rep.memberChains.split(",")}
        self.use_representative_chains = use_representative_chains

        self._interface_chains = {}
        self._chain_residues = {} if not use_representative_chains else None

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
            seq_chain_clustered = pd.DataFrame(seq_chains.residueInfos).dropna()
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

class EPPICLocal(Container):
    IMAGE = "docker://edraizen/eppic-cli"
    PARAMETERS = [
        "-i", ("input_file", "path:in"),
        (":s", "store_true", "-s"),
        (":p", "store_true", "-p"),
        (":l", "store_true", "-l"),
    ]

    def get_entropy_scores(self, pdb_file, maxEntropy=1):
        self(input_file=pdb_file, s=True)
        entropy_file = os,path.join(self.work_dir, "{}.entropies".format(
            os.path.splitext(os.path.basename(pdb_file))))
        assert os.path.isfile(entropy_file), "Failed running EPPIC"
        residueInfo = pd.read_csv(entropy_file, delimeter="\t", comment="#",
            names=["residueNumber", "pdbResidueNumber", "uniProtNumber",
                "residueType", "entropyScore"])
        entropy_scores = defaultdict(lambda: maxEntropy, zip(
            residueInfo.pdbResidueNumber, residueInfo.entropyScore))
        return entropy_scores
