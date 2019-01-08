import os
import subprocess
import time

from molmimic.generate_data.job_utils import cleanup_ids, map_job_rv, map_job

from molmimic.generate_data.iostore import IOStore
from molmimic.parsers.MaxCluster import get_centroid
from molmimic.parsers.haddock import dock
from molmimic.parsers.superpose import align
from molmimic.parsers.CNS import Minimize
from molmimic.generate_data.iRMSD import Complex

def prep(*pdbs, **kwds):
    merge = kwds.get("merge", False)
    work_dir = kwds.get("work_dir", None)

    if work_dir is None:
        work_dir = os.getcwd()

    if merge:
        newf = os.path.join(work_dir, "{}.pdb".format("__".join(
            ["{}_{}".format(f, c) for f, c in pdbs])))
        with open(newf, "w") as new:
            for i, (f, c) in pdbs:
                subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(i), f], stdout=new)
        yield newf

    else:
        for i, (f, c) in pdbs:
            newf = os.path.join(work_dir, "{}.pdb".format(i))
            with open(newf, "w") as new:
                subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(i), f], stdout=new)
            yield newf

def download_pdb(job, sfam_id, pdb, chain, sdi, domNo, in_store=None):
    work_dir = job.fileStore.getLocalTempDir()

    if in_store is None:
        prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
        in_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(sfam_id, pdb, chain, sdi, domNo)
    file = os.path.join(work_dir, os.path.basename(prefix))
    in_store.read_input_file(prefix, file)
    return file

def process_interface(job, inferred_row, int_type):
    mol_file = download_pdb(job, inferred_row.mol_superfam_id, inferred_row.mol_pdb,
        inferred_row.mol_chain, inferred_row.mol_sdi, inferred_row.mol_domNo)

    int_file = download_pdb(job, inferred_row.int_superfam_id, inferred_row.int_pdb,
        inferred_row.int_chain, inferred_row.int_sdi, minferred_row.int_domNo)

    #Set chain names as M=mol, I=int
    mol_file, int_file = list(prep((mol_file, "M"), (int_file, "I")))

    if int_type == "inferred":
        nbr_prefix = "{}/{}_{}_sdi{}_d{}.pdb".format(inferred_row.nbr_superfam_id,
            inferred_row.nbr_pdb, inferred_row.nbr_chain, inferred_row.nbr_sdi,
            inferred_row.nbr_domNo)
        nbr_file = os.path.join(work_dir, nbr_prefix)
        in_store.read_input_file(nbr_prefix, nbr_file)

        #Set chain name: N=nbr
        nbr_file = list(prep((nbr_file, "N")))

        #Reset mol_file with it superposed on top of nbr
        mol_file = align(mol_file, "M", nbr_file, "N")

    return mol_file, int_file

def advanced_dock(job, inferred_row, int_type):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-full-structures".format(prefix))

    dock_name = inferred_row["int_id"]

    mol_file, int_file = process_interface(job, inferred_row, int_type)

    #Perform docking
    docked_file = dock(dock_name, mol_file, "M", inferred_row.mol_resi,
        int_file, "I", inferred_row.int_resi, small_refine=True, job=job)

    #Analyze docking, from PRODIGY
    original_complex = next(prep((mol_file, "M"), (int_file, "I"), merge=True))
    results = analyze(job, docked_file, original_complex)
    return docked_file, results

def analyze(job, complex_path, ref1=None, ref2=None):
    complex = Complex(complex_path)
    results = complex.get_stats(ref=ref1)

    #Molmimic bsa
    #get_bsa()
    job.log(str(results))

    return results

def simple_dock(job, complex_file=None, row=None):
    assert [complex_file, row].count(None)==1
    work_dir = job.fileStore.getLocalTempDir()

    if row is not None:
        mol_file, int_file = process_interface(job, inferred_row, int_type)
        complex_file = next(prep((mol_file, "M"), (int_file, "I"), merge=True))

    minimized_complex_file = Minimize(complex_file, work_dir=work_dir, job=job)
    results = analyze(job, minimized_complex_file, complex_file)
    return minimized_complex_file, results

def cluster(job, mol_sfam, int_sfam, group):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    structure_store = IOStore.get("{}:molmimic-full-structures".format(prefix))
    file_list = os.path.join(work_dir, "{}_{}.txt".format(mol_sfam, int_sfam))
    files = []

    with open(file_list, "w") as f:
        for _, row in group:
            row = row.iloc[0]
            mol_file = download_pdb(job, row.mol_superfam_id, row.mol_pdb,
                row.mol_chain, row.mol_sdi, row.mol_domNo)

            int_file = download_pdb(job, row.int_superfam_id, row.int_pdb,
                row.int_chain, row.int_sdi, row.int_domNo)

            merged_file = next(prep((mol_file, "M"), (int_file, "I"), merge=True))

            print >> f, merged_file
            files.append(merged_file)

    centroid = get_centroid(file_list)
    centroid_index = files.index(centroid)

    return centroid_index, centroid, files

def process_sfam(job, sfam_id, observed=True):
    int_type = "observed" if observed else "inferred"
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    interface_store = IOStore.get("{}:molmimic-interfaces".format(prefix))

    interfaces_key = "{s}/{s}.{o}_interactome".format(
        s=sfam_id, o="observed" if observed else "inferred")
    interfaces_file = os.path.basename(interfaces_key)
    interface_store.read_input_file(interfaces_key, interfaces_file)

    interfaces = pd.read_hdf(interfaces_file, "table")

    if observed:
        for (mol_sfam, int_sfam), group in interfaces.groupby(["mol_supferfam_id", "int_supferfam_id"]):
            #Remove redundant interfaces
            group = group.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
                as_index=False).nth(0).reset_index(drop=True).copy()

            if "mol_sdi" in cdd_interactome:
                key = "mol_sdi"
            elif "mol_sdi_id" in cdd_interactome:
                key = "mol_sdi_id"
            else:
                raise RuntimeError("sdi not in df")

            sdi_group = group.groupby(key, as_index=False)

            centroid_index, centroid, files = cluster(job, mol_sfam, int_sfam, sdi_group)
            for i, (row, file) in enumerate(izip(sdi_group, files)):
                if i == centroid_index:
                    advanced_dock(job, row.iloc[0], int_type)
                else:
                    simple_dock(job, file)
    else:
        for (mol_sfam, int_sfam), group in interfaces.groupby(["mol_supferfam_id", "int_supferfam_id"]):
            best_complex_idx = group["nbr_score"].idxmax()
            best_complex = group.loc[best_complex_idx]
            advanced_dock(job, best_complex, "inferred")

            other_complexes = group[group.index!=best_complex_idx]
            for name, sdi_group in other_complexes.groupby(["mol_sdi", "nbr_obs_int_id"]):
                simple_dock(job, row=sdi_group.iloc[0])

def test_advanced(job, pdb1, pdb2, tbl, settings, face1, face2, work_dir):
    t0 = time.time()
    complex_file = dock("p53-sh3", pdb1, "A", None, pdb2, "B", None, tbl_file=tbl, \
        structures0=10, structures1=2, anastruc1=2, settings_file=settings, \
        work_dir=work_dir, job=job)
    results = analyze(job, complex_file, pdb1, pdb2)
    area = calculate_buried_surface_area(complex_file, "p53-sh3", face1=face1, face2=face2)
    advanced_time = time.now()-t0
    results = results.update(area)
    results["time"] = advanced_time
    results["method"] = "advanced"
    return complex_file, results

def test_simple(job, pdb1, pdb2, face1, face2):
    t0 = time.time()
    complex_file_orig = next(prep((pdb1, "M"), (pdb2, "I"), merge=True))
    complex_file_min = Minimize(complex_file_orig, work_dir=job.fileStore.getLocalTempDir(), job=job)
    results = analyze(job, complex_file_min, pdb1, pdb2, ref1=complex_file_orig)
    area = calculate_buried_surface_area(complex_file, "p53-sh3", face1=face1, face2=face2)
    simple_time = time.now()-t0
    results = results.update(area)
    results["time"] = simple_time
    results["method"] = "simple"
    return complex_file_min, results

def merge(job, *runs):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    store = IOStore.get("{}:molmimic-ddi".format(prefix))

    results = pd.DataFrame(runs)
    runs.to_csv("results.txt")
    store.write_output_file("results.txt", "resulst.txt")

def start_toil_test(job):
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

    face1 = "75 148 149 150 154 155 156 187".split(),
    face2 = "469 472 473 475 476 492 493 494 495 497 498 509 511 513".split()

    #adv = job.addChildJobFn(test_advanced, pdb1_file, pdb2_file, tbl_file, settings_file, face1, face2)
    #simple = job.addChildJobFn(test_simple, pdb1_file, pdb2_file, face1, face2)
    #job.addFollowOnJobFn(merge, (adv, simple))

    job.log(str(os.listdir(work_dir)))

    adv_file, adv_results = test_advanced(job, pdb1_file, pdb2_file, tbl_file, settings_file, face1, face2, work_dir)
    simple_file, simple_results = test_simple(job, pdb1_file, pdb2_file, face1, face2)
    merged1 = analyze(adv_file, ref1=simple_file)
    merged2 = analyze(simple_file, ref1=advanced_file)
    merge(job, (adv_results, simple_results, merged1, merged2))

def start_toil(job):
    work_dir = job.fileStore.getLocalTempDir()
    prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
    in_store = IOStore.get("{}:molmimic-interfaces".format(prefix))
    keys = list(in_store.list_input_directory())
    observed = set(key.split("/")[0] for key in keys if key.endswith("observed_interactome"))
    inferred = set(key.split("/")[0] for key in keys if key.endswith("inferred_interactome"))

    map_job(job, process_sfam, observed)
    map_job(job, process_sfam, inferred, True)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    dataset_name = options.jobStore.split(":")[-1]

    job = Job.wrapJobFn(start_toil_test)
    with Toil(options) as toil:
        toil.start(job)
