import os

import pandas as pd

from molmimic.util.iostore import IOStore
from molmimic.parsers.cath import CATHApi
from molmimic.util.hdf import get_file, filter_hdf, filter_hdf_chunks
from molmimic.util.toil import map_job
from molmimic.util import safe_remove
from molmimic.generate_data import data_stores

from toil.realtimeLogger import RealtimeLogger

def run_cath_hierarchy(job, cathcode, func, cathFileStoreID, **kwds):
    work_dir = job.fileStore.getLocalTempDir()
    further_parallelize = kwds.get("further_parallelize", True)
    level = kwds.get("level", 4)

    if isinstance(cathcode, (int, float, str)):
        #Start with Class if given int, float, or string
        try:
            cathcode = [int(cathocde)]
        except ValueError:
            raise ValueError("Invaid cathcode: {}".format(cathcode))
    elif isinstance(cathcode, (list, tuple)) and len(cathcode)>0:
        if isinstance(cathcode[0], (list, tuple)):
            #Multiple CATH codes
            map_job(job, run_cath_hierarchy, cathcode, func, cathFileStoreID, **kwds)
            return
    elif cathcode is None:
        #Start from the top
        cathcode = []
    else:
        raise ValueError("Invaid cathcode: {}".format(cathcode))


    cath_names = ["class", "architechture", "topology", "homology"]
    if len(cathcode) < level:
        cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
        cathcode = dict(zip(cath_names, cathcode))
        cath_names = cath_names[:len(cathcode)+1]

        cathcodes = filter_hdf_chunks(
            cath_file,
            "table",
            columns=cath_names,
            drop_duplicates=True,
            **cathcode)[cath_names]
    else:
        cathcodes = pd.DataFrame([cathcode], columns=cath_names)

    if cathcodes.shape[1] < level:
        map_job(job, run_cath_hierarchy, cathcodes.values, func, cathFileStoreID, **kwds)
    else:
        sfams = (cathcodes.astype(int).astype(str)+"/").sum(axis=1).str[:-1].tolist()
        RealtimeLogger.info("Running sfam {}".format(cathcode))
        kwds.pop("further_parallelize", None)
        kwds.pop("level", None)
        map_job(job, func, sfams, cathFileStoreID, **kwds)

def download_cath_domain(cath_domain, sfam_id=None, work_dir=None):
    """Download CATH domain from CATh API. Raises KeyError is cath_domain doesn't
    exist"""
    if work_dir is None:
        work_dir = os.getcwd()

    cath_store = data_stores.cath_api_service
    cath = CATHApi(cath_store, work_dir=work_dir)

    domain_file = cath.get_domain_pdb_file(cath_domain)

    return domain_file
