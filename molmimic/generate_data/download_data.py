import os
import subprocess

from toil.common import Toil
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.job_utils import map_job

#from molmimic.generate_data.iostore import IOStore

databases = {
    "pdb": ["rsync.wwpdb.org::ftp_data/structures/divided/pdb/",
            ["-rlpt", "-v", "-z", "--delete", "--port=33444"]],
    "sifts": ["rsync.ebi.ac.uk::pub/databases/msd/sifts/split_xml/",
            ["-vrltH", "--delete", "-D", "--numeric-ids"]]
}

def get_jobstore_name(job, name="raw-files"):
    return "{}-{}".format(job.fileStore.jobStore.config.jobStore,
        name.lower())

def call_rsync(arguments, inplace=True, callback=None, args=None):
    assert isinstance(arguments, list)

    if not callable(callback):
        callback = lambda x: False
        args = None

    if isinstance(args, (list,tuple)):
        args = list(args)
    elif args is not None:
        args = [args]
    else:
        args = []

    if inplace:
        command = ["rsync", "--inplace"]+arguments
    else:
        command = ["rsync"]+arguments

    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        while True:
            line = p.stdout.readline().rstrip().decode().encode('ascii', errors='ignore')
            if line == '' and p.poll() is None:
                break

            if callback(*[line]+args):
                break
    except:
        p.terminate()
        p.kill()
        raise

    final = p.stdout.read()
    p.terminate()
    p.kill()
    return final

def compare_rsync_to_store(job, name, options=None, url=None, ending=".gz", no_retries=False, can_retry=True):
    if name in databases:
        _url, _options = databases[name]
        url = url or _url
        options = options or _options
    elif url is None or options is None or not isinstance(options, (list,tuple)):
        raise RuntimeError("Invalid call to compare_rsync_to_store")

    work_dir = job.fileStore.getLocalTempDir()

    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    jobStore = IOStore.get("{}:molmimic-{}".format(prefix, name))

    outpath = os.path.join(work_dir, name)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    num_saved = jobStore.get_number_of_items()

    arguments = list(options)+["--dry-run", "--stats", url, outpath]

    total_num = None
    def get_total_num(line, total):
        job.log(line)
        if line.startswith("Number of files:"):
            job.log("SAVING LINE")
            total = int(line.split(": ")[2][-5].replace(",", ""))
            return True

    output = call_rsync(arguments, callback=get_total_num, args=[total_num])

    for line in output.split("\n"):
        get_total_num(line, total_num)

    job.log("COMPARE Saved: {}, Total: {}".format(num_saved, total_num))
    job.log("OUTPUT: {}".format(output))
    if can_retry and num_saved<total_num:
        job.log("RERUNNING {}".format(name))
        sync_job = job.addChildJobFn(rsync_to_store, name, options, url, ending)
        sync_job.addFollowOnJobFn(compare_rsync_to_store, name, options, url, ending,
            can_retry=False)

    return num_saved, total_num

def rsync_to_store(job, name, options=None, url=None, file_ending=".gz"):
    if name in databases:
        _url, _options = databases[name]
        url = url or _url
        options = options or _options
    elif url is None or options is None or not isinstance(options, (list,tuple)):
        raise RuntimeError("Invalid call to rync_rstore")

    work_dir = job.fileStore.getLocalTempDir()

    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    jobStore = IOStore.get("{}:molmimic-{}".format(prefix, name))

    outpath = os.path.join(work_dir, name)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    arguments = list(options)+[url, outpath]

    def add_to_jobstore(key):
        if key and key.endswith(file_ending):
            full_path = os.path.join(outpath, key)
            if os.path.isfile(full_path):
                jobStore.write_output_file(full_path, key)

    call_rsync(arguments, callback=add_to_jobstore)

def download_pdb(job):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "pdb"))

    pdb_path = os.path.join(work_dir, "pdb")
    os.makedirs(pdb_path)
    subprocess.check_call(["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444",
        "rsync.wwpdb.org::ftp_data/structures/divided/pdb/", pdb_path])

    add_directory_to_jobstore(jobStore, pdb_path)

def download_sifts(job, force=False):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    jobStore = IOStore.get("{}:molmimic-sifts".format(prefix))

    sifts_path = os.path.join(work_dir, "sifts")
    if not os.path.isdir(sifts_path):
        os.makedirs(sifts_path)
    subprocess.check_call(["rsync", "-vrltH", "--delete", "--stats", "-D", "--numeric-ids",
        "rsync.ebi.ac.uk::pub/databases/msd/sifts/split_xml/", sifts_path])

    add_directory_to_jobstore(jobStore, sifts_path)

def download_ibis_obs(job, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "ibis"))

    ibis_obs_path = os.path.join(work_dir, "IBIS_observed.h5")
    subprocess.check_call(["wget", "-N", "-O", ibis_obs_path,
        "https://www.dropbox.com/s/47agor1gx0qewr0/IBIS_observed.h5?dl=0"])
    jobStore.write_output_file(ibis_obs_path, "IBIS_observed.h5")

def download_ibis_inf(job, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "ibis"))

    ibis_inf_path = os.path.join(work_dir, "IBIS_inferred.h5")
    subprocess.check_call(["wget", "-N", "-O", ibis_inf_path,
        "https://www.dropbox.com/s/0rzk2yyspurqbnw/IBIS_inferred.h5?dl=0"])
    jobStore.write_output_file(ibis_inf_path, "IBIS_inferred.h5")

def download_mmdb(job, preemptable=True):
    work_dir = job.fileStore.getLocalTempDir()
    jobStore = IOStore.get(get_jobstore_name(job, "ibis"))

    mmdb_path = os.path.join(job._fileStore.getLocalTempDir(), "MMDB.h5")
    subprocess.check_call(["wget", "-N", "-O", mmdb_path,
        "https://www.dropbox.com/s/l42w7qq3kixq4v9/MMDB.h5?dl=0"])
    jobStore.write_output_file(mmdb_path, "MMDB.h5")

def download_consurf(job):
    import shutil, requests, zipfile
    from io import StringIO
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
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]

    # for name in ("pdb"): #, "sifts"):
    #     jobStore = IOStore.get("{}:molmimic-{}".format(prefix, name))
    #     if jobStore.get_number_of_items() > 0:
    #         #If exists don't rerun
    #         job.addChildJobFn(compare_rsync_to_store, name, can_retry=name!="pdb")
    #     else:
    #         job.addChildJobFn(compare_rsync_to_store, name)

    job.addChildJobFn(download_sifts)

    #job.addChildJobFn(download_consurf)

    #Create jobStore/bucket before so all of these can be saved into
    # IOStore.get(get_jobstore_name(job, "IBIS")).connect()
    # job.addChildJobFn(download_ibis_obs)
    # job.addChildJobFn(download_ibis_inf)
    # job.addChildJobFn(download_mmdb)

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
