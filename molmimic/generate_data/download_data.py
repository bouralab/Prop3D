import os
import subprocess
from util import data_path_prefix

from toil.common import Toil

def get_jobstore(job, name="raw-files"):
    locator = "{}-{}".format(job._fileStore.jobStore.locator, name)
    return Toil.getJobStore(locator)

def add_directory_to_jobstore(jobStore, direcotry):
    remove_dir = None
    for dirpath, dirnames, filenames in os.walk(direcotry):
        if remove_dir is None:
            remove_dir = dirpath+os.sep
        for f in filenames:
            key = os.path.join(dirpath, f)
            jobStore.importFile(key, sharedFileName=key[len(remove_dir):])

def download_pdb(job, force=False):
    jobStore = get_jobstore(job, "pdb")

    try:
        jobStore.getSize("xy/1xyz.ent.gz")
        if not force:
            return
    except RuntimeError:
        pass

    pdb_path = os.path.join(job._fileStore.getLocalTempDir(), "pdb")
    subprocess.call(["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444",
        "rsync.wwpdb.org::ftp_data/structures/divided/pdb/", pdb_path])

    add_directory_to_jobstore(jobStore, pdb_path)

def download_sifts(job, force=False):
    jobStore = get_jobstore(job, "sifts")

    try:
        jobStore.getSize("xy/1xyz.xml.gz")
        if not force:
            return
    except RuntimeError:
        pass

    sifts_path = os.path.join(job._fileStore.getLocalTempDir(), "sifts")
    subprocess.call(["rsync", "-vrltH", "--delete", "--stats", "-D", "--numeric-ids",
        "rsync.ebi.ac.uk::pub/databases/msd/sifts/split_xml/", sifts_path])

    add_directory_to_jobstore(jobStore, sifts_path)

def download_ibis_obs(job, force=False):
    jobStore = get_jobstore(job, "IBIS")

    try:
        jobStore.getSize("IBIS_observed.h5)
        if not force:
            return
    except RuntimeError:
        pass

    ibis_obs_path = os.path.join(job._fileStore.getLocalTempDir(), "IBIS_observed.h5")
    subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/47agor1gx0qewr0/IBIS_observed.h5?dl=0"])
    jobStore.importFile(ibis_obs_path, sharedFileName="IBIS_observed.h5")


def download_ibis_inf(job, force=False):
    jobStore = get_jobstore(job, "IBIS")

    try:
        jobStore.getSize("IBIS_inferred.h5")
        if not force:
            return
    except RuntimeError:
        pass

    ibis_inf_path = os.path.join(job._fileStore.getLocalTempDir(), "IBIS_inferred.h5")
    subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/0rzk2yyspurqbnw/IBIS_inferred.h5?dl=0"])
    jobStore.importFile(ibis_inf_path, sharedFileName="IBIS_inferred.h5")

def download_mmdb(job, force=False):
    jobStore = get_jobstore(job, "IBIS")

    try:
        jobStore.getSize("MMDB.h5")
        if not force:
            return
    except RuntimeError:
        pass

    mmdb_path = os.path.join(job._fileStore.getLocalTempDir(), "MMDB.h5")
    subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/l42w7qq3kixq4v9/MMDB.h5?dl=0"])
    jobStore.importFile(mmdb_path, sharedFileName="MMDB.h5")

def download_consurf(job):
    jobStore = get_jobstore(job, "ConSurf")

    try:
        jobStore.getSize("XY/1XYZ_A")
        if not force:
            return
    except RuntimeError:
        pass

    job.log("START CONSURF")
    import shutil, requests, zipfile
    from cStringIO import StringIO
    from molmimic.parsers.Consurf import download_consurf as download_consurf_all

    #Download nr mapping
    pdb_nr = os.path.join(job._fileStore.getLocalTempDir(), "pdbaa_list.nr")
    r = requests.get("http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip", stream=True)
    with zipfile.ZipFile(StringIO(r.content)) as z:
        with z.open("pdbaa_list.nr") as zf, open(pdb_nr, "w") as f:
            shutil.copyfileobj(zf, f)
    jobStore.importFile(pdb_nr, sharedFileName="pdbaa_list.nr")

    #Download all pdb consurf files
    consurf_path = os.path.join(data_path_prefix, "ConSurf")
    download_consurf(consurf_path=job._fileStore.getLocalTempDir())
    job.log("DONE CONSURF")
    add_directory_to_jobstore(jobStore, consurf_path)

def start_toil(job, name="download_data"):
    job.addChildJobFn(download_pdb)
    job.addChildJobFn(download_sifts)
    job.addChildJobFn(download_consurf)

    #Create jobStore/bucket before so all of these can be saved into
    jobStore = get_jobstore(job, "IBIS")
    job.addChildJobFn(download_ibis_obs)
    job.addChildJobFn(download_ibis_inf)
    job.addChildJobFn(download_mmdb)
