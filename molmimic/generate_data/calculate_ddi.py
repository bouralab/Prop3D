import os
from molmimic.generate_data.iostore import IOStore
from molmimic.parsers.haddock import dock
from molmimic.parsers.ce import align

def inferred_dock(job, dock_name, superfam_inf, pdb_inf, chain_inf, sdi_inf, \
  domNo_inf, resi_inf, superfam_obs, pdb_obs, chain_obs, sdi_obs, domNo_obs, resi_obs, align_chain):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    inf_prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(superfam_inf, pdb_inf, chain_inf,
        sdi_inf, domNo_inf)
    inf_file = os.path.join(work_dir, inf_prefix)
    in_store.read_input_file(inf_prefix, inf_file)

    obs_prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(superfam_obs, pdb_obs, chain_obs,
        sdi_obs, domNo_obs)
    obs_file = os.path.join(work_dir, obs_prefix)
    in_store.read_input_file(obs_prefix, obs_file)

    aligned_inf_pdb = align(inf_file, chain_inf, obs_file, align_chain)
    dock(dock_name, aligned_inf_pdb, chain_inf, resi_inf, obs_file, chain_obs,
        domNo_inf, structures0=10, structures1=2, anastruc1=2, job=job)

def observed_dock(job, dock_name, superfam1, pdb1, chain1, sdi1, domNo1, resi1, \
  superfam2, pdb2, chain2, sdi2, domNo2, resi2):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    prefix1 = "{}/{}_{}_sdi{}_d{}.pdb".format(superfam1, pdb1, chain1, sdi1, domNo1)
    file1 = os.path.join(work_dir, prefix1)
    in_store.read_input_file(prefix1, file1)

    prefix2 = "{}/{}_{}_sdi{}_d{}.pdb".format(superfam2, pdb2, chain2, sdi2, domNo2)
    file2 = os.path.join(work_dir, prefix2)
    in_store.read_input_file(prefix2, file2)

    dock(dock_name, file1, chain1, resi1, file2, chain2, domNo2, structures0=10,
        structures1=2, anastruc1=2, job=job)

def start_toil(job):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-ddi".format(prefix))

    pdb1 = "1YCS_A_sdi225433_d0.pdb"
    pdb1_file = os.path.join(work_dir, pdb1)
    in_store.read_input_file(pdb1, pdb1_file)

    pdb2 = "1YCS_B_sdi225436_d2.pdb"
    pdb2_file = os.path.join(work_dir, pdb2)
    in_store.read_input_file(pdb2, pdb2_file)

    tbl = "1YCS_B_sdi225436_d2_1YCS_A_sdi225433_d0.tbl"
    tbl_file = os.path.join(work_dir, tbl)
    in_store.read_input_file(tbl, tbl_file)

    settings = "new.html"
    settings_file = os.path.join(work_dir, settings)
    in_store.read_input_file(settings, settings_file)

    dock("p53-sh3", pdb1, "A", None, pdb2, "B", None, tbl_file=tbl_file, \
        structures0=10, structures1=2, anastruc1=2, settings_file=settings_file,  \
        work_dir=work_dir, job=job)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    dataset_name = options.jobStore.split(":")[-1]

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as toil:
        toil.start(job)
