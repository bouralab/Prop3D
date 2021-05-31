import os

import pandas as pd

from molmimic.util.iostore import IOStore
from molmimic.parsers.cath import CATHApi
from molmimic.util.hdf import get_file, filter_hdf, filter_hdf_chunks
from molmimic.util.toil import map_job
from molmimic.util import safe_remove
from molmimic.generate_data import data_stores

from toil.realtimeLogger import RealtimeLogger

def fix_cathcode(c):
    if isinstance(cathcode, (int, float, str)):
        try:
            return int(c)
        except ValueError:
            if isinstance(cathcode, str) and cathcode in ["None", "all"]:
                return None

def run_cath_hierarchy(job, cathcode, func, cathFileStoreID, *args, **kwds):
    RealtimeLogger.info("Run Hierachy: {}".format(cathcode))
    work_dir = job.fileStore.getLocalTempDir()

    skip_cathcode = kwds.get("skip_cathcode", None)
    further_parallelize = kwds.get("further_parallelize", True)
    level = kwds.get("level", 4)

    if isinstance(cathcode, (int, float, str)):
        #Start with Class if given int, float, or string
        try:
            cathcode = [int(cathcode)]
        except ValueError:
            if isinstance(cathcode, str):
                if cathcode in ["None", "all"]:
                    #Start from the top
                    cathcode = []
                elif "." in cathcode:
                    parts = cathcode.split(".")
                    _cathcode = []
                    for p in parts:
                        try:
                            _cathcode.append(int(p))
                        except ValueError:
                            raise ValueError("Invaid cathcode: {}".format(cathcode))
                    cathcode = _cathcode
                else:
                    raise ValueError("Invaid cathcode: {}".format(cathcode))
            else:
                raise ValueError("Invaid cathcode: {}".format(cathcode))
    elif isinstance(cathcode, (list, tuple)):
        if len(cathcode)>1 and isinstance(cathcode[0], (list, tuple)):
            #Multiple CATH codes
            RealtimeLogger.info("Run Mult Hierachy")
            map_job(job, run_cath_hierarchy, cathcode, func, cathFileStoreID, *args, **kwds)
            return
        if len(cathcode)==1 and isinstance(cathcode[0], (list, tuple)):
            run_cath_hierarchy(job, cathcode[0], func, cathFileStoreID, *args, **kwds)
            return
        elif len(cathcode)==1 and cathcode[0] in [None, "None", "all"]:
            cathcode = []
        else:
            #This is correct, e.g.
            #    Full cath code: [2., 40., 60., 10.]
            #    Jast class: [2.]
            #    or all superfamilies: []
            RealtimeLogger.info("Run correct Hierachy: {}".format(cathcode))
            try:
                cathcode = list(map(int, cathcode))
                RealtimeLogger.info("Run correct fixed Hierachy: {}".format(cathcode))
            except ValueError:
                raise ValueError("Invaid cathcode: {}".format(cathcode))

    elif cathcode is None:
        #Start from the top
        cathcode = []
    else:
        raise ValueError("Invaid cathcode: {}".format(cathcode))

    RealtimeLogger.info("Running cathcodexx: {}".format(cathcode))

    cath_names = ["class", "architechture", "topology", "homology"]
    if len(cathcode) < level:
        RealtimeLogger.info("Loading cathcide 0 {}".format(cathcode))
        cathStoreID = kwds.get("cathCodeStoreID", cathFileStoreID)
        cath_file = job.fileStore.readGlobalFile(cathStoreID, cache=True)

        remove_store = False
        if "cathCodeStoreID" not in kwds:
            cathcodes = filter_hdf_chunks(
                cath_file,
                "table",
                columns=cath_names,
                drop_duplicates=True)
            safe_remove(cath_file, warn=True)
            cath_file = os.path.join(work_dir, "cathCodeStore.h5")
            cathcodes.to_hdf(cath_file, "table", format="table",
                complevel=9, complib="bzip2", min_itemsize=1024)
            kwds["cathCodeStoreID"] = job.fileStore.writeGlobalFile(cath_file)
            remove_store = True
            del cathcodes

        RealtimeLogger.info("Loading cathcide {}".format(cathcode))
        cathcode = dict(zip(cath_names, cathcode))
        RealtimeLogger.info("Loading cathcide 2 {}".format(cathcode))
        RealtimeLogger.info("Using logger {}".format(RealtimeLogger.getLogger()))

        cath_names = cath_names[:len(cathcode)+1]
        cathcodes = filter_hdf_chunks(
            cath_file,
            "table",
            columns=cath_names,
            drop_duplicates=True,
            **cathcode)[cath_names]

        rc = safe_remove(cath_file, warn=True)
        if not rc:
            RealtimeLogger.info("Unable to remove file {}")

    else:
        RealtimeLogger.info("skipped")
        cathcodes = pd.DataFrame([cathcode], columns=cath_names)
        RealtimeLogger.info("skipped {}".format(cathcode))

    if cathcodes.shape[1] < level:
        map_job(job, run_cath_hierarchy, cathcodes.values.tolist(), func,
            cathFileStoreID, *args, **kwds)
    else:
        RealtimeLogger.info("Running sfam {}".format(cathcodes))
        RealtimeLogger.info("Running sfam {}".format(cathcode))

        sfams = (cathcodes.astype(int).astype(str)+"/").sum(axis=1).str[:-1]

        if skip_cathcode is not None and len(skip_cathcode)==0:
            skip_cathcode = sfams[~sfams.isin(skip_cathcode)]

        try:
            del kwds["skip_cathcode"]
        except KeyError:
            pass

        sfams = sfams.tolist()
        RealtimeLogger.info("Running sfam {}".format(sfams))
        RealtimeLogger.info("Running func {}".format(func))
        kwds.pop("further_parallelize", None)
        kwds.pop("level", None)
        if "cathCodeStoreID" in kwds:
            del kwds["cathCodeStoreID"]
        RealtimeLogger.info("Running map job {}".format(map_job))
        RealtimeLogger.info("Using logger {}".format(id(RealtimeLogger.getLogger())))

        RealtimeLogger.info("Successors0: {}".format(list(job.description.allSuccessors())))
        map_job(job, func, sfams, cathFileStoreID, *args, **kwds)
        RealtimeLogger.info("Successors1: {}".format(list(job.description.allSuccessors())))

    del cathcodes

def download_cath_domain(cath_domain, sfam_id=None, work_dir=None):
    """Download CATH domain from CATh API. Raises KeyError is cath_domain doesn't
    exist"""
    if work_dir is None:
        work_dir = os.getcwd()

    cath_store = data_stores.cath_api_service
    cath = CATHApi(cath_store, work_dir=work_dir)

    domain_file = cath.get_domain_pdb_file(cath_domain)

    return domain_file
