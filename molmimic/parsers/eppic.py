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
        print(chain.tail(100))
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

class EPPICLocal(object):
    @staticmethod
    def get_blast_db(job):
        work_dir = job.fileStore.getLocalTempDir()
        store = IOStore.get("aws:us-east-1:molmimic-uniprot")

        blast_db_key = "uniref100.fasta"
        blast_db_prefix = os.path.join(work_dir, "uniref100.fasta")

        blast_db_ext = ["phr", "pin", "pog", "psd", "psi", "psq"]

        blast_db_exists = [store.exists(blast_db_key+"."+ext) for ext in blast_db_ext]

        if all(blast_db_exists):
            blastFileStoreIDs = []
            for ext in blast_db_ext:
                fname = blast_db_prefix+"."+ext
                store.read_input_file(blast_db_key+"."+ext, fname)
                blastFileStoreIDs.append(job.fileStore.writeGlobalFile(fname))

                try:
                    os.remove(fname)
                except OSError:
                    pass
        else:
            #Download Uniref100
            uniref_100_key = "uniref100.fasta.gz"
            uniref_100_file_gz = os.path.join(work_dir, "uniref100.fasta.gz")
            uniref_100_file = os.path.join(work_dir, "uniref100.fasta")

            if store.exists(uniref_100_key):
                store.read_input_file(uniref_100_key, uniref_100_file_gz)
            else:
                url = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"

                try:
                    with requests.get(url, stream=True) as r, open(uniref_100_file_gz, 'wb') as f:
                        r.raw.read = partial(r.raw.read, decode_content=True)
                        copyfileobj(r.raw, f)
                except Exception as e:
                    RealtimeLogger.info("uniprot Error {}: {}".format(type(e), e))
                    return None

                store.write_output_file(uniref_100_file, uniref_100_key)

            #Unzip
            with gzip.open(uniref_100_file_gz, 'rb') as fin, open(uniref_100_file, 'wb') as fout:
                shutil.copyfileobj(fin, fout)

            try:
                os.remove(uniref_100_file_gz)
            except OSError:
                pass

            #Make blast db
            if not os.path.abspath(os.path.dirname(pdb_path)) == os.path.abspath(work_dir):
                shutil.copy(pdb_file, work_dir)

            parameters = ["-in", uniref_100_file, "-dbtype", "prot", "-logfile",
                "makeblastdb.log", "-parse_seqids", "-out", uniref_100_file,
                "-title", uniref_100_key[:-3]]

            try:
                apiDockerCall(job,
                              image='edraizen/blast:latest',
                              working_dir=work_dir,
                              parameters=parameters)
            except ContainerError as e:
                RealtimeLogger.info("Could not create BlastDB")
                raise
                return None

            try:
                os.remove(uniref_100_file_gz)
            except OSError:
                pass

            assert all([os.path.isfile(uniref_100_file+"."+ext) for ext in blast_db_ext])

            blastFileStoreIDs = []
            for ext in blast_db_ext:
                fname = uniref_100_file+"."+ext
                store.write_output_file(fname, uniref_100_key+"."+ext)
                blastFileStoreIDs.append(job.fileStore.writeGlobalFile(fname))

                try:
                    os.remove(fname)
                except OSError:
                    pass

        return blastFileStoreIDs

def download_eppic(*key, **kwds):
    service = kwds.get("service", "sequences")
    work_dir = kwds.get("work_dir")

    if work_dir is None:
        work_dir = os.getcwd()

    service_key = "/".join(k.lower() if isinstance(k, str) else str(k) for k in key)
    url = "http://www.eppic-web.org/rest/api/v3/job/{}/{}".format(service, service_key)
    RealtimeLogger.info("Downloading {}".format(url))

    key = "{}/{}.json".format(service, service_key)
    fname = os.path.join(work_dir, "{}.json".format(service_key.replace("/", "-")))

    try:
        with requests.get(url, stream=True) as r, open(fname, 'wb') as f:
            r.raw.read = partial(r.raw.read, decode_content=True)
            copyfileobj(r.raw, f)
    except Exception as e:
        RealtimeLogger.info("EPPIC Error {}: {}".format(type(e), e))
        return None

    store = IOStore.get("aws:us-east-1:molmimic-eppic-service")
    store.write_output_file(fname, key)

    with open(fname) as f:
        RealtimeLogger.info(f.read())

    return fname

def run_eppic(*key, **kwds):
    service = kwds.get("service", "sequences")
    work_dir = kwds.get("work_dir")
    download = kwds.get("download", True)
    attempts = kwds.get("attempts", 2)

    if work_dir is None:
        work_dir = os.getcwd()

    service_key = "/".join(k.lower() if isinstance(k, str) else str(k) for k in key)
    url = "http://www.eppic-web.org/rest/api/v3/job/{}/{}".format(service, service_key)
    RealtimeLogger.info("Downloading {}".format(url))
    eppic_key = "{}/{}.json".format(service, service_key)
    eppic_file = os.path.join(work_dir, "{}-{}.json".format(service,
        service_key.replace("/", "-")))

    store = IOStore.get("aws:us-east-1:molmimic-eppic-service")

    if os.path.isfile(eppic_file):
        RealtimeLogger.info("EPPIC read from file")
        should_remove = False
    elif store.exists(eppic_key):
        RealtimeLogger.info("EPPIC get from store")
        store.read_input_file(eppic_key, eppic_file)
        should_remove = True
    else:
        should_remove = True
        if download:
            eppic_file = download_eppic(*key, service=service, work_dir=work_dir)
            if not eppic_file:
                RealtimeLogger.info("DOWNLOAD EPPIC -- FAILED")
                #Error no EPPIC file
                return {"file":None, "key":eppic_key, "should_remove":False, "data":None}
        else:
            RealtimeLogger.info("DOWNLOAD EPPIC -- NO DOWNLOAD")
            return {"file":None, "key":eppic_key, "should_remove":False, "data":None}

    try:
        with open(eppic_file) as f:
            pass
    except IOError as e:
        #Might be empty
        try:
            os.remove(eppic_file)
        except OSError:
            pass
        RealtimeLogger.info("Failed reading, {} bc {}".format(eppic_file, e))
        if attempt > 0:
            return run_eppic(*key, service=service, work_dir=work_dir,
                download=download, attempts=attempts-1)
        else:
            return {"file":None, "key":eppic_key, "should_remove":False, "data":None}

    try:
        with open(eppic_file) as f:
            eppic = json.load(f)
    except (SystemExit, KeyboardInterrupt) as e:
        raise
    except ValueError as e:
        rerun = False
        with open(eppic_file) as f:
            for line in f:
                RealtimeLogger.info("EPPIC Failed {}".format(line))
                if "<html>" in line or "<head>" in line:
                    RealtimeLogger.info("Is restarting")
                    store.remove_file(eppic_key)
                    rerun = attempts > 0
                if "Too many submissions" in line:
                    RealtimeLogger.info("Is chainging IP")
                    try:
                        reset_ip()
                    except (SystemExit, KeyboardInterrupt):
                        raise
                    except:
                        pass

        try:
            os.remove(eppic_file)
        except OSError:
            pass

        if rerun:
            return run_eppic(*key, service=service, work_dir=work_dir,
                download=download, attempts=attempts-1)
        else:
            RealtimeLogger.info("Not restarting")
            return {"file":None, "key":eppic_key, "should_remove":False, "data":None}
    except Exception as e:
        RealtimeLogger.info("EPPIC Failed parsing json ({}): {}".format(type(e), e))
        return {"file":None, "key":eppic_key, "should_remove":False, "data":None}

    if should_remove:
        try:
            os.remove(eppic_file)
        except OSError:
            pass

    return {"file":eppic_file, "key":eppic_key, "should_remove":False, "data":eppic}
