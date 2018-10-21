import os
from molmimic.generate_data.iostore import IOStore
from molmimic.parsers.haddock import dock

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

    dock("p53-sh3", pdb1, "A", None, pdb2, "B", None, tbl_file=tbl_file, settings_file=settings_file, work_dir=work_dir, job=job)

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
