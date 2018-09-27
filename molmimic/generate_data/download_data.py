import os
import subprocess

from toil.common import Toil
from iostore import IOStore
#from molmimic.generate_data.iostore import IOStore

def get_jobstore_name(job, name="raw-files"):
    return = "{}-{}".format(job._fileStore.jobStore.locator, name)

def add_directory_to_jobstore(jobStore, directory):
    for dirpath, dirnames, filenames in os.walk(direcotry):
        for f in filenames:
            fpath = os.path.join(dirpath, f)
            key = fpath[len(direcotry):]
            jobStore.write_output_file(fpath, key)

def download_pdb(job):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "pdb"))

    pdb_path = os.path.join(work_dir, "pdb")
    subprocess.call(["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444",
        "rsync.wwpdb.org::ftp_data/structures/divided/pdb/", pdb_path])

    add_directory_to_jobstore(jobStore, pdb_path)

def download_sifts(job, force=False):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "pdb"))

    sifts_path = os.path.join(work_dir, "sifts")
    subprocess.call(["rsync", "-vrltH", "--delete", "--stats", "-D", "--numeric-ids",
        "rsync.ebi.ac.uk::pub/databases/msd/sifts/split_xml/", sifts_path])

    add_directory_to_jobstore(jobStore, sifts_path)

def download_ibis_obs(job):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "IBIS"))

    ibis_obs_path = os.path.join(work_dir), "IBIS_observed.h5")
    subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/47agor1gx0qewr0/IBIS_observed.h5?dl=0"])
    jobStore.write_output_file(ibis_obs_path, "IBIS_observed.h5")

def download_ibis_inf(job):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "IBIS"))

    ibis_inf_path = os.path.join(work_dir, "IBIS_inferred.h5")
    subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/0rzk2yyspurqbnw/IBIS_inferred.h5?dl=0"])
    jobStore.write_output_file(ibis_obs_path, "IBIS_inferred.h5")

def download_mmdb(job):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "IBIS"))

    mmdb_path = os.path.join(job._fileStore.getLocalTempDir(), "MMDB.h5")
    subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/l42w7qq3kixq4v9/MMDB.h5?dl=0"])
    jobStore.write_output_file(mmdb_path, "MMDB.h5")

def download_consurf(job):
    import shutil, requests, zipfile
    from cStringIO import StringIO
    from molmimic.parsers.Consurf import download_consurf as download_consurf_all

    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "ConSurf"))

    #Download nr mapping
    pdb_nr = os.path.join(work_dir, "pdbaa_list.nr")
    r = requests.get("http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip", stream=True)
    with zipfile.ZipFile(StringIO(r.content)) as z:
        with z.open("pdbaa_list.nr") as zf, open(pdb_nr, "w") as f:
            shutil.copyfileobj(zf, f)
    jobStore.write_output_file(pdb_nr, "pdbaa_list.nr")

    #Download all pdb consurf files
    consurf_path = os.path.join(work_dir, "ConSurf")
    os.path.mkdirs(consurf_path)
    download_consurf(consurf_path=consurf_path)
    add_directory_to_jobstore(jobStore, consurf_path)

def start_toil(job, name="download_data"):
    job.addChildJobFn(download_pdb)
    job.addChildJobFn(download_sifts)
    #job.addChildJobFn(download_consurf)

    #Create jobStore/bucket before so all of these can be saved into
    IOStore.get(get_jobstore_name(job, "IBIS"))
    job.addChildJobFn(download_ibis_obs)
    job.addChildJobFn(download_ibis_inf)
    job.addChildJobFn(download_mmdb)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "never"

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
