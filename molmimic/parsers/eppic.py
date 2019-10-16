from __future__ import print_function
import os
from shutil import copyfileobj
import json
from collections import defaultdict
from functools import partial

import requests
import numpy as np

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import natural_keys, reset_ip

from toil.realtimeLogger import RealtimeLogger

from joblib import Parallel, delayed

def download_eppic(*key, service="sequences", work_dir=None):
    if work_dir is None:
        work_dir = os.getcwd()

    service_key = "/".join(k.lower() if isinstance(k, str) else str(k) for k in key)
    url = "http://www.eppic-web.org/rest/api/v3/job/{}/{}".format(service_key, key)

    key = "{}/{}.json".format(service, service_key)
    fname = os.path.join(work_dir, "{}.json".format(service_key.replace("/", "-")))

    try:
        with requests.get(url, stream=True) as r, open(fname, 'wb') as f:
            r.raw.read = partial(r.raw.read, decode_content=True)
            copyfileobj(r.raw, f)
    except Exception as e:
        RealtimeLogger.info("EPPIC Error {}: {}".format(type(e), e))
        return None

    store = IOStore.get("aws:us-east-1:molmimic-eppic")
    store.write_output_file(fname, key)

    with open(fname) as f:
        RealtimeLogger.info(f.read())

    return fname

def run_eppic(*key, service="sequence", work_dir=None, download=True, attempts=2):
    if work_dir is None:
        work_dir = os.getcwd()

    service_key = "/".join(k.lower() if isinstance(k, str) else str(k) for k in key)
    url = "http://www.eppic-web.org/rest/api/v3/job/{}/{}".format(service_key, key)
    eppic_key = "{}/{}.json".format(service, service_key)
    eppic_file = os.path.join(work_dir, "{}.json".format(service_key.replace("/", "-")))

    store = IOStore.get("aws:us-east-1:molmimic-eppic")

    if os.path.isfile(eppic_file):
        RealtimeLogger.info("EPPIC read from file")
        should_remove = False
    elif store.exists(eppic_db_key):
        RealtimeLogger.info("EPPIC get from store")
        store.read_input_file(eppic_key, eppic_file)
        should_remove = True
    else:
        should_remove = True
        if download:
            RealtimeLogger.info("DOWNLOAD EPPIC")
            eppic_file = download_eppic(*key, service=service, work_dir=work_dir)
            if not eppic_file:
                RealtimeLogger.info("DOWNLOAD EPPIC -- FAILED")
                #Error no EPPIC file
                return {"file":None, "key":eppic_key, "should_remove":False, "data":None}
        else:
            RealtimeLogger.info("DOWNLOAD EPPIC -- NOT")
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
        with open(eppic_db_file) as f:
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

    eppic["key"] = eppic_key

    if should_remove:
        try:
            os.remove(eppic_file)
        except OSError:
            pass

    return {"file":eppic_file, "key":eppic_key, "should_remove":False, "data":eppic}

def get_entropy_scores(pdb, chain, work_dir=None, download=True, attempts=2):
    result = run_eppic(pdb, "sequence", work_dir=work_dir,
        download=download, attempts=attempts)

    entropy_scores = defaultdict(lambda: np.nan)

    if result["file"] is None or result["data"] is None:
        return entropy_scores
    else:
        eppic = result["data"]

    if "chainClusters" in eppic:
        eppic = eppic["chainClusters"]["chainCluster"]

    for chainCluster in eppic:
        if chain in chainCluster["memberChains"]:
            for residueInfo in chainCluster["residueInfos"]:
                if "pdbResidueNumber" in residueInfo:
                    try:
                        resseq_parts = natural_keys(residueInfo["pdbResidueNumber"], use_int=True)
                        resseq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
                    except IndexError as e:
                        RealtimeLogger.info("Error parsing line {}".format(e))
                        break

                    if "entropyScore" in residueInfo:
                        entropy_scores[resseq] = residueInfo["entropyScore"]
                    else:
                        entropy_scores[resseq] = 1
    del eppic
    return entropy_scores

def get_interfaces(pdb, bio=False, work_dir=None, download=True, attempts=2):
    result = run_eppic(pdb, service="interfaces", work_dir=work_dir,
        download=download, attempts=attempts)

    interfaces = pd.DataFrame(result["data"])
    interfaces = interfaces.assign(interfaceType=
        interfaces["interfaceScores"].apply(lambda x:
            [score["callName"] for score in x["interfaceScore"] \
                if score["method"] == "eppic"][0]))

    if bio:
        interfaces = interfaces[interfaces["interfaceType"]=="bio"]

    return interfaces

def get_interface_residues(pdb, interfaceId, work_dir=None, download=True, attempts=2):
    result = run_eppic(pdb, interfaceId, service="interfaces", work_dir=work_dir,
        download=download, attempts=attempts)

    return pd.DataFrame(result["data"])

def get_contacts(pdb, interfaceId, work_dir=None, download=True, attempts=2):
    result = run_eppic(pdb, interfaceId, service="contacts", work_dir=work_dir,
        download=download, attempts=attempts)

    return result["data"]

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
