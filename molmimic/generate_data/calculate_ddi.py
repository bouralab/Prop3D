import os
import subprocess

from molmimic.generate_data.iostore import IOStore
from molmimic.parsers.haddock import dock
from molmimic.parsers.superpose import align

def prep(*pdbs, work_dir=None):
    if work_dir is None:
        work_dir = os.getcwd()

    for i, (f, c) in pdbs:
        newf = os.path.join(work_dir, "{}.pdb".format(i))
        with open(newf, "w") as new:
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(i), f], stdout=new)
        yield newf

def binary_dock(job, inferred_row, int_type):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    dock_name = inferred_row["int_id"]

    mol_prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(inferred_row.mol_superfam_id,
        inferred_row.mol_pdb, inferred_row.mol_chain, inferred_row.mol_sdi,
        inferred_row.mol_domNo)
    mol_file = os.path.join(work_dir, mol_prefix)
    in_store.read_input_file(mol_prefix, mol_file)

    int_prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(inferred_row.int_superfam_id,
        inferred_row.int_pdb, inferred_row.int_chain, inferred_row.int_sdi,
        inferred_row.int_domNo)
    int_file = os.path.join(work_dir, obs_prefix)
    in_store.read_input_file(int_prefix, int_file)

    #Set chain names as M=mol, I=int
    mol_file, int_file = list(prep((mol_file, "M"), (int_file, "I")))

    if int_type = "inferred":
        nbr_prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(inferred_row.nbr_superfam_id,
            inferred_row.nbr_pdb, inferred_row.nbr_chain, inferred_row.nbr_sdi,
            inferred_row.nbr_domNo)
        nbr_file = os.path.join(work_dir, nbr_prefix)
        in_store.read_input_file(nbr_prefix, nbr_file)

        #Set chain name: N=nbr
        nbr_file = list(prep((nbr_file, "N")))

        #Reset mol_file with it superposed on top of nbr
        mol_file = align(mol_file, "M", nbr_file, "N")

    #Perform docking
    docked_file = dock(dock_name, mol_nbr_file, "M", inferred_row.mol_resi,
        int_file, "I", inferred_row.int_resi, small_refine=True, job=job)

    #Analyze docking, from PRODIGY

    # Parse structure
    structure, n_chains, n_res = parse_structure(struct_path)
    prodigy = Prodigy(structure, cmd.selection, cmd.temperature)
    prodigy.predict(distance_cutoff=cmd.distance_cutoff, acc_threshold=cmd.acc_threshold)

    #Molmimic bsa


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
