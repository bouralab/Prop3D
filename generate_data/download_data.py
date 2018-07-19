import os
import subprocess
from util import data_path_prefix

def download_pdb(job):
    pdb_path = os.path.join(data_path_prefix, "pdb", "pdb")
    if not os.path.isdir(pdb_path):
        subprocess.call(["rsync", "-rlpt", "-v", "-z", "--delete", "--port=33444",
            "rsync.wwpdb.org::ftp_data/structures/divided/pdb/", pdb_path])

def download_sifts(job):
    sifts_path = os.path.join(data_path_prefix, "pdb", "sifts")
    if not os.path.isdir(sifts_path):
        subprocess.call(["rsync", "-vrltH", "--delete", "--stats", "-D", "--numeric-ids",
            "rsync.ebi.ac.uk::pub/databases/msd/sifts/split_xml/", sifts_path])

def download_ibis_obs(job):
    ibis_obs_path = os.path.join(data_path_prefix, "IBIS_observed.h5")
    if not os.path.isfile(ibis_obs_path):
        subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
            "https://www.dropbox.com/s/47agor1gx0qewr0/IBIS_observed.h5?dl=0"])

def download_ibis_inf(job):
    ibis_inf_path = os.path.join(data_path_prefix, "IBIS_inferred.h5")
    if not os.path.isfile(ibis_inf_path):
        subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
            "https://www.dropbox.com/s/0rzk2yyspurqbnw/IBIS_inferred.h5?dl=0"])

def download_mmdb(job):
    mmdb_path = os.path.join(data_path_prefix, "MMDB.h5")
    if not os.path.isfile(mmdb_path):
        subprocess.call(["wget", "-N", "-nc", "-O", ibis_obs_path,
            "https://www.dropbox.com/s/l42w7qq3kixq4v9/MMDB.h5?dl=0"])

def download_consurf(job):
    job.log("START CONSURF")
    import shutil, requests, zipfile
    from cStringIO import StringIO
    from molmimic.parsers.Consurf import download_consurf as download_consurf_all

    #Downlaod nr mapping
    pdb_nr = os.path.join(data_path_prefix, "ConSurf", "pdbaa_list.nr")
    if not os.path.isfile(pdb_nr):
        r = requests.get("http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip", stream=True)
        with zipfile.ZipFile(StringIO(r.content)) as z:
            with z.open("pdbaa_list.nr") as zf, open(pdb_nr, "w") as f:
                shutil.copyfileobj(zf, f)

    #Download all pdb consurf files
    consurf_path = os.path.join(data_path_prefix, "ConSurf")
    if True: #len(next(os.walk(consurf_path))[1]) == 0:
        # mmdb_path = os.path.join(data_path_prefix, "MMDB.h5")
        # pdbs = pd.read_hdf(mmdb_path, "StrucutralDomains")[["pdbId", "chnLett"]].drop_duplicates().dropna()
        # pdbs.apply(lambda r: download_consurf(r["pdbId"], r["chnLett"]), axis=1)
        download_consurf_all()
    job.log("DONE CONSURF")
    return

def start_toil(job, name="download_data"):
    job.addChildJobFn(download_pdb)
    job.addChildJobFn(download_sifts)
    job.addChildJobFn(download_ibis_obs)
    job.addChildJobFn(download_ibis_inf)
    mmdb = job.addChildJobFn(download_mmdb)
    mmdb.addFollowOnJobFn(download_consurf)
