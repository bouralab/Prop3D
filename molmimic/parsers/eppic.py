from __future__ import print_function
import os
from shutil import copyfileobj
import json
from collections import defaultdict

import requests
import numpy as np

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import natural_keys

from toil.realtimeLogger import RealtimeLogger

from joblib import Parallel, delayed

def download_eppic(pdb, eppic_path=None):
    if eppic_path is None:
        eppic_path = os.getcwd()
        
    for pdb in pdbs["pdb"].drop_duplicates():
        fname = "{}_eppic.json".format(pdb)
        url = "http://www.eppic-web.org/ewui/ewui/dataDownload?type=json&id={}&withResInfo=t&withSeqInfo=t".format(pdb.lower())
        try:
            with requests.get(url, stream=True) as r, open(fname, 'wb') as f:
                r.raw.read = functools.partial(r.raw.read, decode_content=True)
                shutil.copyfileobj(r.raw, f)
            print(pdb)
        except Exception as e:
            print(pdb, "FAILED", e.__class__.__name__, e)

def run_eppic(pdb, chain, work_dir=None, download=False):
    if work_dir is None:
        work_dir = os.getcwd()

    pdb_id = pdb.upper()
    eppic_db_key = pdb_id+"_eppic.json"
    eppic_db_file = os.path.join(work_dir, eppic_db_key)

    store = IOStore.get("aws:us-east-1:molmimic-eppic")
    
    entropy_scores = defaultdict(lambda: np.nan)
    
    if os.path.isfile(eppic_db_file):
        should_remove = False
    elif store.exists(eppic_db_key):
        store.read_input_file(eppic_db_key, eppic_db_file)
        should_remove = False
    else:
        assert 0, ("No eppic file", eppic_db_key)
        if download:
            eppic_db_file = download_eppic(pdb, chain)
            if not eppic_db_file:
                #Error no EPPIC file
                return entropy_scores
        else:
            return entropy_scores

    try:
        with open(eppic_db_file) as f:
            pass
    except IOError as e:
        #Might be empty
        RealtimeLogger.info("Failed reading, {} bc {}".format(eppic_db_file, e))
        return entropy_scores

    
    with open(eppic_db_file) as f:
        try:
            eppic = json.load(f)
        except json.decoder.JSONDecodeError:
            return entropy_scores
    
    for chainCluster in eppic["chainClusters"]["chainCluster"]:
        if chain in chainCluster["memberChains"]:
            for residueInfo in chainCluster["residueInfos"]:
                if "pdbResidueNumber" in residueInfo and "entropyScore" in residueInfo:
                    try:
                        resseq_parts = natural_keys(residueInfo["pdbResidueNumber"], use_int=True)
                        resseq = (" ", int(resseq_parts[1]), resseq_parts[2].rjust(1))
                        entropy_scores[resseq] = residueInfo["entropyScore"]
                    except IndexError as e:
                        RealtimeLogger.info("Error parsing line {}".format(e))
                        break
    del eppic
    
    if should_remove:
        try:
            os.remove(eppic_db_file)
        except OSError:
            pass
        
    return entropy_scores

