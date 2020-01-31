from molmimic.util.iostore import IOStore
from molmimic.parsers.cath import CATHApi
from molmimic.util.hdf import get_file, filter_hdf, filter_hdf_chunks

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
            map_job(job, run_cath_hierarchy, cathcodes, func, cathFileStoreID, **kwds)
            return
    else:
        raise ValueError("Invaid cathcode: {}".format(cathcode))


    cath_path = get_file(job, "cath.h5", cathFileStoreID, work_dir=work_dir)

    cath_names = ["class", "architechture", "topology", "homology"]
    cathcode = dict(zip(cath_names, cathcode))
    cath_names = cath_names[:len(cathcode)+1]

    cathcodes = filter_hdf_chunks(
        cath_path,
        "table",
        columns=cath_names,
        drop_duplicates=True,
        **cathcode)[cath_names]

    RealtimeLogger.info("cathcode {} {} {}".format(cathcode, cath_names, cathcodes.columns))

    if len(cathcodes.columns) < level:
        map_job(job, run_cath_hierarchy, cathcodes.values, func, cathFileStoreID, **kwds)
        RealtimeLogger.info("Running {} {}s".format(len(cathcodes), cathcodes.columns[-1]))
    else:
        sfams = (cathcodes.astype(int).astype(str)+"/").sum(axis=1).str[:-1].tolist()
        RealtimeLogger.info("Running sfam {}".format(cathcode))
        kwds.pop("further_parallelize")
        kwds.pop("level")
        map_job(job, func, sfams, update_features, **kwds)

    try:
        os.remove(cath_path)
    except (FileNotFoundError, OSError):
        pass

def download_cath_domain(job, cath_domain, sfam_id, check=False, work_dir=None):
    if work_dir is None and job is None:
        work_dir = os.getcwd()
    elif job is not None:
        work_dir = job.fileStore.getLocalTempDir()

    RealtimeLogger.info("RUNNING {}/{}".format(sfam_id, cath_domain))

    cath_store = IOStore.get("aws:us-east-1:molmimic-cath-service")
    cath = CATHApi(cath_store, work_dir=work_dir)

    domain_file = cath.get_domain_pdb_file(cath_domain)
