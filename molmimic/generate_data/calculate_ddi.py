import os
import sys
import time
import logging

import pandas as pd

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from botocore.exceptions import ClientError

from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.iRMSD import Complex

def prep(*pdbs, **kwds):
    import subprocess
    from molmimic.generate_data.util import PDB_TOOLS
    merge = kwds.get("merge", False)
    work_dir = kwds.get("work_dir", None)

    if work_dir is None:
        work_dir = os.getcwd() #pdbs[0][0][:-1*len(os.path.basename(pdbs[0][0]))]

    if merge:
        newf = os.path.join(work_dir, "{}.pdb".format("__".join(
            ["{}_{}".format(os.path.splitext(os.path.basename(f))[0], c) for f, c in pdbs])))
        with open(newf+".temp", "w") as new:
            for i, (f, c) in enumerate(pdbs):
                print f, c
                subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(c), f], stdout=new)
        with open(newf, "w") as new:
            subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py"), newf+".temp"], stdout=new)
        try:
            os.remove(newf+".temp")
        except OSError:
            pass

        yield newf

    else:
        for i, (f, c) in enumerate(pdbs):
            newf = os.path.join(work_dir, "{}_{}.pdb".format(os.path.splitext(os.path.basename(f))[0], c))
            with open(newf, "w") as new:
                subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-{}".format(c), f], stdout=new)
            yield newf

def download_pdb(job, sfam_id, pdb, chain, sdi, domNo, extracted=False, work_dir=None):
    if work_dir is None:
        work_dir = os.getcwd()

    in_store = IOStore.get("aws:us-east-1:molmimic-full-structures")

    prefix = "{}/{}/{}_{}_sdi{}_d{}.pdb".format(int(sfam_id), pdb[1:3].lower(),
        pdb.upper(), chain, sdi, domNo)
    file = os.path.join(work_dir, os.path.basename(prefix))
    if extracted:
        prefix += ".extracted"
    assert in_store.exists(prefix)
    in_store.read_input_file(prefix, file)
    return file

def process_interface(job, row, int_type, work_dir=None):
    from molmimic.parsers.superpose import align

    if work_dir is None:
        work_dir = os.getcwd()

    RealtimeLogger.info("ROW {}".format(row))

    try:
        mol_file = download_pdb(job, row.mol_superfam_id, row.mol_pdb,
            row.mol_chain, row.mol_sdi_id, row.mol_domNo, extracted=True, work_dir=work_dir)
    except (ClientError, AssertionError):
        RealtimeLogger.info("Could not download mol {} {} {} {}".format(
            row.mol_superfam_id, row.mol_pdb,
            row.mol_chain, row.mol_sdi_id))
        return None, None, None, None

    try:
        int_file = download_pdb(job, row.int_superfam_id, row.int_pdb,
            row.int_chain, row.int_sdi_id, row.int_domNo, extracted=True, work_dir=work_dir)
    except (ClientError, AssertionError):
        RealtimeLogger.info("Could not download int {} {} {} {}".format(
            row.mol_superfam_id, row.mol_pdb,
            row.mol_chain, row.mol_sdi_id))
        return None, None, None, None

    #Set chain names as M=mol, I=int
    mol_file, int_file = list(prep((mol_file, "M"), (int_file, "I"), work_dir=work_dir))
    mol_resi = map(int, row["mol_res"].split(","))
    int_resi = map(int, row["int_res"].split(","))

    if int_type == "inferred":
        try:
            nbr_file = download_pdb(job, row.nbr_superfam_id, row.nbr_pdb,
                row.nbr_chain, row.nbr_sdi_id, row.nbr_domNo)
        except (ClientError, AssertionError):
            RealtimeLogger.info("Could not download nbr {} {} {} {}".format(
                row.nbr_superfam_id, row.nbr_pdb,
                row.nbr_chain, row.nbr_sdi_id))
            return None, None, None, None

        #Set chain name: N=nbr
        nbr_file = next(prep((nbr_file, "N")))

        #Reset mol_file with it superposed on top of nbr
        mol_file = align(mol_file, "M", nbr_file, "N", work_dir=work_dir, job=job)

    RealtimeLogger.info("MOL FILE: "+mol_file)
    RealtimeLogger.info("INT FILE: "+int_file)

    return mol_file, mol_resi, int_file, int_resi

def advanced_dock(job, mol_file, mol_resi, int_file, int_resi, original_complex, int_id=None, work_dir=None, cores=2):
    from molmimic.parsers.haddock import dock

    if work_dir is None:
        work_dir = os.getcwd()

    #dock_name = row["obs_int_id"] if int_type=="observed" else "fix_me"

    #mol_file, mol_resi, int_file, int_resi = process_interface(job, row, int_type, work_dir=work_dir)

    fname = "{}__{}".format(
        os.path.splitext(os.path.basename(mol_file))[0],
        os.path.splitext(os.path.basename(int_file))[0])
    dock_name = "{}_{}".format(int_id if int_id is not None else "inf", fname)

    # #Perform docking
    # docked_file = dock(dock_name, mol_file, "M", inferred_row.mol_resi,
    #     int_file, "I", inferred_row.int_resi, refine=True, small_refine=True, job=job)
    #
    # #Analyze docking, from PRODIGY
    # original_complex = next(prep((mol_file, "M"), (int_file, "I"), merge=True))
    # results = analyze(job, docked_file, original_complex)
    # return docked_file, results

    t0 = time.time()

    #Refine:HADDOCK #Score:Haddock
    orig_file, complex_file, haddock_result, haddock_files = dock(
        dock_name,
        mol_file, "M", mol_resi,
        int_file, "I", int_resi,
        small_refine=True,
        work_dir=work_dir,
        cores=cores,
        job=job)

    refined_complex = Complex(complex_file, face1=mol_resi, face2=int_resi,
        method="advanced", work_dir=work_dir, job=job)
    refined_complex.results.append(haddock_result)
    refined_complex.compare_to_reference(original=original_complex)
    refined_complex["time"] = time.time()-t0
    #Refine:None/Original #Score:CNS
    # energy_begin = calculate_energy(orig_file, work_dir=work_dir, job=job)
    # energy_begin = energy_begin.rename(lambda l: "begin_"+l)
    # results = results.append(energy_begin)

    #Refine:HADDOCK #Score:CNS
    # energy_after = calculate_energy(complex_file, work_dir=work_dir, job=job)
    # energy_after = energy_after.rename(lambda l: "after_"+l)
    # results = results.append(energy_after)

    #Refine:None/Original #Score:Haddock
    # haddock_begin = score_complex(orig_file, "MI", work_dir=work_dir, job=job)
    # haddock_begin = haddock_begin.rename(lambda l: "begin_"+l)
    # results = results.append(haddock_begin)

    # advanced_time = time.time()-t0
    #
    # results["time"] = advanced_time
    # results["method"] = "advanced"
    return refined_complex, haddock_files

# def analyze(job, complex_path, chain1="M", chain2="I", face1=None, face2=None, ref1=None, ref2=None, work_dir=None):
#     from molmimic.generate_data.iRMSD import Complex
#     complex = Complex(complex_path, chain1=chain1, chain2=chain2, face1=face1, face2=face2, work_dir=work_dir, job=job)
#     results = complex.get_stats(ref=)
#     RealtimeLogger.info(str(results))
#     return results

def simple_dock(job, orig_file, mol_resi, int_resi, original_complex, work_dir=None):
    from molmimic.parsers.CNS import Minimize

    if work_dir is None:
        work_dir = os.getcwd()

    t0 = time.time()
    # mol_file, mol_resi, int_file, int_resi = process_interface(job, row, int_type, work_dir=work_dir)
    # orig_file = next(prep((mol_file, "M"), (int_file, "I"), merge=True, work_dir=work_dir))

    # minimized_complex_file = Minimize(complex_file, work_dir=work_dir, job=job)
    # import shutil
    # shutil.copy(minimized_complex_file, "/root/minimized_complex.pdb")
    # assert os.path.isfile("/root/minimized_complex.pdb")
    # results = analyze(job, minimized_complex_file, complex_file, work_dir=work_dir)
    # return minimized_complex_file, results

    #t0 = time.time()
    #complex_file_orig = next(prep((pdb1, "M"), (pdb2, "I"), merge=True))

    #Refine:CNS #Score:CNS
    complex_file, cns_results = Minimize(orig_file, work_dir=work_dir, job=job)
    refined_complex = Complex(complex_file, face1=mol_resi, face2=int_resi,
        method="simple", work_dir=work_dir, job=job)
    refined_complex.results.append(cns_results)
    refined_complex.compare_to_reference(original=original_complex)
    refined_complex["time"] = time.time()-t0

    #Refine:None/Original #Score:CNS
    # energy_begin = calculate_energy(orig_file, work_dir=work_dir, job=job)
    # energy_begin = energy_begin.rename(lambda l: "begin_"+l)
    # results = results.append(energy_begin)

    #Refine:None/Original #Score:Haddock
    # haddock_begin = score_complex(orig_file, "MI", work_dir=work_dir, job=job)
    # haddock_begin = haddock_begin.rename(lambda l: "begin_"+l)
    # results = results.append(haddock_begin)

    #Refine:CNS #Score:Haddock
    # haddock_after = score_complex(complex_file, "MI", work_dir=work_dir, job=job)
    # haddock_after = haddock_after.rename(lambda l: "after_"+l)
    # results = results.append(haddock_after)

    # simple_time = time.time()-t0
    #
    # results["method"] = "simple"
    return refined_complex

def get_original_complexes(job, mol_sfam, int_sfam, group, int_type, work_dir=None):
    if work_dir is None:
        work_dir = os.cwd()

    complex_files = []
    for _, row in group:
        row = row.iloc[0]
        RealtimeLogger.info("ROW: {}".format(row))

        mol_file, mol_resi, int_file, int_resi = process_interface(job,
            row, int_type, work_dir=work_dir)
        # try:
        #     mol_file = download_pdb(job, row.mol_superfam_id, row.mol_pdb,
        #         row.mol_chain, row.mol_sdi_id, row.mol_domNo, work_dir=work_dir)
        #
        #     int_file = download_pdb(job, row.int_superfam_id, row.int_pdb,
        #         row.int_chain, row.int_sdi_id, row.int_domNo, work_dir=work_dir)
        # except (KeyboardInterrupt, SystemExit):
        #     raise
        # except Exception as e:
        if mol_file is None or int_file is None:
            #PDB files not found, skip
            RealtimeLogger.info("Cannot download PDB {}.{}.{} bc it was not found".format(row.mol_pdb, row.mol_chain, row.mol_sdi_id))
            complex_files.append(None)
            continue

        merged_file = next(prep((mol_file, "M"), (int_file, "I"), merge=True,
            work_dir=work_dir))

        complex_files.append(merged_file)

    return complex_files

def cluster(job, mol_sfam, int_sfam, pdb_complexes, work_dir=None):
    from molmimic.parsers.MaxCluster import run_maxcluster, get_centroid
    from molmimic.generate_data.util import PDB_TOOLS, SubprocessChain

    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")

    if work_dir is None:
        work_dir = os.getcwd()

    cluster_dir = os.path.join(work_dir, "{}_{}_cluster".format(int(mol_sfam), int(int_sfam)))
    if not os.path.exists(cluster_dir):
        os.makedirs(cluster_dir)

    file_list = os.path.join(cluster_dir, "{}_{}_file_list.txt".format(int(mol_sfam), int(int_sfam)))

    with open(file_list, "w") as f:
        for pdb_file in pdb_complexes:
            if not pdb_file: continue
            cmds = [
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-A", pdb_file],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_reres.py")],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_reatom.py")],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
            ]

            pdb_to_cluster = os.path.join(cluster_dir, os.path.basename(pdb_file))
            with open(pdb_to_cluster, "w") as pdb:
                SubprocessChain(cmds, pdb)

            print >> f, os.path.basename(pdb_to_cluster)

    log_file = os.path.join(cluster_dir, "{}_{}.max_cluster_logs".format(int(mol_sfam), int(int_sfam)))
    distance_file = os.path.join(cluster_dir, "{}_{}_distances.txt".format(int(mol_sfam), int(int_sfam)))

    ddi_base = "{}_{}".format(int(mol_sfam), int(int_sfam))

    logs = run_maxcluster("rmsd", file_list=file_list, log=True, R=distance_file,
        work_dir=cluster_dir, C=1, P=10, job=job)

    RealtimeLogger.info("CLUSTER_DIR: {}".format(os.listdir(cluster_dir)))
    RealtimeLogger.info("LOG FILE: {}".format(logs))

    with open(logs) as f:
        RealtimeLogger.info("LOG IS {}".format(next(f)))

    ddi_store.write_output_file(file_list, "{}/{}".format(ddi_base, os.path.basename(file_list)))
    ddi_store.write_output_file(logs, "{}/{}".format(ddi_base, os.path.basename(log_file)))
    ddi_store.write_output_file(distance_file, "{}/{}".format(ddi_base, os.path.basename(distance_file)))

    centroid = get_centroid(logs, job=job)
    centroid_index = None
    for i, p in enumerate(pdb_complexes):
        if p and p.endswith(centroid):
            return i, p
    return None, centroid
#

# def test_advanced(job, pdb1, pdb2, face1, face2, work_dir):
#     from molmimic.parsers.haddock import dock
#     t0 = time.time()
#     RealtimeLogger.info("ADV")
#     orig_file, complex_file, haddock_result, haddock_files = dock(
#         "p53-sh3",
#         pdb1, "A", face1,
#         pdb2, "B", face2, \
#         small_refine=True,
#         work_dir=work_dir, job=job)
#     results = analyze(job, complex_file, chain1="A", chain2="B", face1=face1, face2=face2, ref1=orig_file, work_dir=work_dir)
#     advanced_time = time.now()-t0
#     results = results.append(haddock_result)
#     energy_begin = calculate_energy(orig_file, work_dir=work_dir, job=job)
#     energy_begin = energy_begin.rename(lambda l: "begin_"+l)
#     results = results.append(energy_begin)
#     results["time"] = advanced_time
#     results["method"] = "advanced"
#     return complex_file, results
#
# def test_simple(job, pdb1, pdb2, face1, face2, work_dir):
#     from molmimic.parsers.CNS import Minimize, calculate_energy
#     t0 = time.time()
#     complex_file_orig = next(prep((pdb1, "M"), (pdb2, "I"), merge=True))
#     complex_file_min, cns_results = Minimize(complex_file_orig, work_dir=work_dir, job=job)
#     import shutil
#     shutil.copy(complex_file_min, "/root/minimized_complex.pdb")
#     assert os.path.isfile("/root/minimized_complex.pdb")
#     results = analyze(job, complex_file_min, face1=face1, face2=face2, ref1=complex_file_orig, work_dir=work_dir)
#     results = results.append(cns_results)
#     energy_begin = calculate_energy(complex_file_orig, work_dir=work_dir, job=job)
#     energy_begin = energy_begin.rename(lambda l: "begin_"+l)
#     results = results.append(energy_begin)
#     simple_time = time.time()-t0
#     results["time"] = simple_time
#     results["method"] = "simple"
#     return complex_file_orig, complex_file_min, results

# def merge(job, *runs):
#     work_dir = job.fileStore.getLocalTempDir()
#     prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
#     store = IOStore.get("aws:us-east-1:molmimic-ddi")
#
#     assert 0, runs
#
#     results = pd.DataFrame(runs)
#     runs.to_csv("results.txt")
#     store.write_output_file("results.txt", "resulst.txt")

# def start_toil_test(job):
#     RealtimeLogger.info("START")
#
#     work_dir = job.fileStore.getLocalTempDir()
#     #prefix = job.fileStore.jobStore.config.jobStore.rsplit(":", 1)[0]
#     in_store = IOStore.get("aws:us-east-1:molmimic-ddi")
#
#     pdb1 = "1YCS_A_sdi225433_d0.pdb"
#     pdb1_file = os.path.join(work_dir, pdb1)
#     in_store.read_input_file(pdb1, pdb1_file)
#
#     pdb2 = "1YCS_B_sdi225436_d2.pdb"
#     pdb2_file = os.path.join(work_dir, pdb2)
#     in_store.read_input_file(pdb2, pdb2_file)
#
#     tbl = "1YCS_B_sdi225436_d2_1YCS_A_sdi225433_d0.tbl"
#     tbl_file = os.path.join(work_dir, tbl)
#     in_store.read_input_file(tbl, tbl_file)
#
#     settings = "new.html"
#     settings_file = os.path.join(work_dir, settings)
#     in_store.read_input_file(settings, settings_file)
#
#     naccess = "naccess.config"
#     naccess_file = os.path.join(work_dir, naccess)
#     in_store.read_input_file(naccess, naccess_file)
#     os.environ["FREESASA_PAR"] = os.path.join(work_dir, "naccess.config")
#
#     face1 = map(int, "75 148 149 150 154 155 156 187".split())
#     face2 = map(int, "469 472 473 475 476 492 493 494 495 497 498 509 511 513".split())
#     #adv = job.addChildJobFn(test_advanced, pdb1_file, pdb2_file, tbl_file, settings_file, face1, face2)
#     #simple = job.addChildJobFn(test_simple, pdb1_file, pdb2_file, face1, face2)
#     #job.addFollowOnJobFn(merge, (adv, simple))
#
#     RealtimeLogger.info(str(os.listdir(work_dir)))
#
#     adv_file, adv_results = test_advanced(job, pdb1_file, pdb2_file, face1, face2, work_dir)
#     orig_file, simple_file, simple_results = test_simple(job, pdb1_file, pdb2_file, face1, face2, work_dir)
#     merged1 = analyze(adv_file, ref1=simple_file)
#     merged2 = analyze(simple_file, ref1=advanced_file)
#     merge(job, (adv_results, simple_results, merged1, merged2))

def process_sfam(job, sfam_id, int_sfams=None, observed=True, min_cluster=0, cores=2):
    from itertools import izip
    int_type = "observed" if observed else "inferred"
    work_dir = job.fileStore.getLocalTempDir()
    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")



    interfaces_key = "{s}/{s}.{o}_interactome".format(
        s=sfam_id, o="observed" if observed else "inferred")
    interfaces_file = os.path.basename(interfaces_key)
    interface_store.read_input_file(interfaces_key, interfaces_file)

    interfaces = pd.read_hdf(interfaces_file, "table")

    if observed:
        num_groups = 0
        igroups = interfaces.fillna(-1.).groupby(["mol_superfam_id", "int_superfam_id"])

        for (mol_sfam, int_sfam), group in igroups:
            ddi_base = "{}_{}".format(int(mol_sfam), int(int_sfam))
            if any([k for k in ddi_store.list_input_directory(ddi_base) if k.endswith("results.csv")]):
                continue

            if group.shape[0] < min_cluster:
                #Only keep groups with more than threshold
                continue

            if int_sfams != -1. and int_sfam not in int_sfams:
                #Only keep groups where known sfams are specified
                pass

            #Remove redundant interfaces
            group = group.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
                as_index=False).nth(0).reset_index(drop=True).copy()

            if "mol_sdi" in group.columns:
                key = "mol_sdi"
            elif "mol_sdi_id" in group.columns:
                key = "mol_sdi_id"
            else:
                raise RuntimeError("sdi not in df")

            sdi_group = group.groupby(key, as_index=False)
            files = get_original_complexes(job, mol_sfam, int_sfam, sdi_group, "observed", work_dir=work_dir)

            if len(files) == 0 or files.count(None) >= len(files)-1:
                RealtimeLogger.info("No PDBs for {}:{}".format(mol_sfam, int_sfam))
                continue
            elif len(files) == 1:
                centroid_index, centroid = 0, files[0]
                sdi_group = [next(sdi_group)]
            else:
                centroid_index, centroid = cluster(job, mol_sfam, int_sfam, files, work_dir=work_dir)

            for i, ((sdi_key, row), file) in enumerate(izip(sdi_group, files)):
                if i == centroid_index:
                    mol_file, mol_resi, int_file, int_resi = process_interface(job, row.iloc[0], int_type, work_dir=work_dir)
                    orig_file = next(prep((mol_file, "M"), (int_file, "I"), merge=True, work_dir=work_dir))
                    original_complex = Complex(orig_file, face1=mol_resi, face2=int_resi,
                        method="original", work_dir=work_dir, job=job)

                    advanced_complex, adv_files = advanced_dock(job, mol_file,
                        mol_resi, int_file, int_resi, original_complex,
                        int_id=row.iloc[0]["obs_int_id"], work_dir=work_dir, cores=cores)
                    simple_complex = simple_dock(job, orig_file, mol_resi, int_resi,
                        original_complex, work_dir=work_dir)
                    simple_complex.compare_to_reference(adv_ref=advanced_complex)
                    advanced_complex.compare_to_reference(simple_ref=simple_complex)
                    # orig_to_adv = {"ref2_"+k:v for k, v in simple_file.compare(adv_file).iteritems()}
                    # adv_to_orig = {"ref2_"+k:v for k, v in adv_file.compare(simple_file).iteritems()}
                    # simple_results = _simple_results.append(pd.Series(orig_to_adv))
                    # adv_results = _adv_results.append(pd.Series(adv_to_orig))
                    all_results = pd.concat([
                        original_complex.set_prefix(),
                        simple_complex.set_prefix(),
                        advanced_complex.set_prefix()]).to_frame().T
                    all_results_path = os.path.join(work_dir, "{}.csv".format(row.iloc[0]["obs_int_id"]))
                    all_results.to_csv(all_results_path)

                    ddi_key = "{}/{}".format(ddi_base, row.iloc[0]["obs_int_id"])
                    ddi_store.write_output_file(orig_file, "{}/original.pdb".format(ddi_key))
                    ddi_store.write_output_file(advanced_complex.pdb, "{}/haddock.pdb".format(ddi_key))
                    ddi_store.write_output_file(adv_files, "{}/haddock_info.tgz".format(ddi_key))
                    ddi_store.write_output_file(simple_complex.pdb, "{}/cns.pdb".format(ddi_key))
                    ddi_store.write_output_file(all_results_path, "{}/results.csv".format(ddi_key))



    # else:
    #     for (mol_sfam, int_sfam), group in interfaces.groupby(["mol_supferfam_id", "int_supferfam_id"]):
    #         best_complex_idx = group["nbr_score"].idxmax()
    #         best_complex = group.loc[best_complex_idx]
    #         advanced_dock(job, best_complex, "inferred")
    #
    #         other_complexes = group[group.index!=best_complex_idx]
    #         for name, sdi_group in other_complexes.groupby(["mol_sdi", "nbr_obs_int_id"]):
    #             simple_dock(job, row=sdi_group.iloc[0])

def process_sfams(job, max_sfams=300, memory="1G"):
    import json
    from collections import defaultdict

    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-ddi")

    json_file = os.path.join(work_dir, "sorted_sfams.json")
    store.read_input_file("sorted_sfams.json", json_file)

    with open(json_file) as f:
        counts = json.load(f)

    sfams = defaultdict(set)
    num_ddi = 0
    for i, ((mol_sfam, int_sfam), count) in enumerate(counts):
        RealtimeLogger.info("{}: {}-{}".format(i, mol_sfam, int_sfam))
        if -1 in map(int, (mol_sfam, int_sfam)):
            RealtimeLogger.info("SKIPPED")
            continue

        if count<90:
            RealtimeLogger.info("SKIPPED count is less than 100")
            continue

        print (mol_sfam, int_sfam)
        #if i<30: continue
        # and (max_sfams is None or num_ddi<max_sfams) and \
      #mol_sfam not in sfams[int_sfam]:
        if mol_sfam != int_sfam:
            if (mol_sfam not in sfams and int_sfam not in sfams) or (mol_sfam in sfams and int_sfam in sfams):
                sfams[mol_sfam].add(int_sfam)
                RealtimeLogger.info("Added {}: {}-{}".format(i, mol_sfam, int_sfam))
            elif mol_sfam in sfams and int_sfam not in sfams:
                sfams[mol_sfam].add(int_sfam)
                RealtimeLogger.info("Added {}: {}-{}".format(i, mol_sfam, int_sfam))
            elif mol_sfam not in sfams and int_sfam in sfams:
                sfams[int_sfam].add(mol_sfam)
                RealtimeLogger.info("Added {}: {}-{}".format(i, int_sfam, mol_sfam))
            else:
                RealtimeLogger.info("Could not add {}: {}-{}; mol_sfam {} sfams; int_sfam {} sfams".format(
                    i, mol_sfam, int_sfam, "in" if mol_sfam in sfams else "not in",
                    "in" if int_sfam in sfams else "not in"))
            num_ddi += 1
            #break

    RealtimeLogger.info("{} starting domains".format(len(sfams)))
    #mol_sfam = list(sfams.keys())[pd.np.random.randint(len(sfams))]
    #int_sfams = sfams[mol_sfam]
    #job.addChildJobFn(process_sfam, mol_sfam, int_sfams)

    for mol_sfam, int_sfams in sfams.iteritems():
        job.addChildJobFn(process_sfam, mol_sfam, int_sfams)

def best_sfams(job, all_counts, max_sfams=300):
    import json
    work_dir = job.fileStore.getLocalTempDir()
    out_store = IOStore.get("aws:us-east-1:molmimic-ddi")

    #Merge into one dataframe
    counts = pd.concat(all_counts)

    #mol->int should be same as int->mol: remove dupes
    ddi_counts = {}
    for counts in all_counts:
        for row in counts.itertuples():
            ddi = tuple(map(int, sorted((row.mol_superfam_id, row.int_superfam_id))))
            if ddi in ddi_counts:
                RealtimeLogger.info("{} {}, are counts symmetrical? {}".format(ddi[0], ddi[1],
                    "Yes" if ddi_counts[ddi] == row.count else "No"))
                continue
            ddi_counts[ddi] = row.count

    sfams = sorted(ddi_counts.iteritems(), key=lambda x: x[1], reverse=True)
    RealtimeLogger.info("sfams is {}".format(sfams))
    sfam_file = os.path.join(work_dir, "sorted_sfams.json")
    with open(sfam_file, "w") as f:
        json.dump(sfams, f)
    out_store.write_output_file(sfam_file, "sorted_sfams.json")

    return sfams[:max_sfams]

def get_sfam_ddi_sizes(job, sfam_id, observed=True):
    int_type = "observed" if observed else "inferred"
    work_dir = job.fileStore.getLocalTempDir()
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces")

    interfaces_key = "{s}/{s}.{o}_interactome".format(
        s=sfam_id, o="observed" if observed else "inferred")
    interfaces_file = os.path.basename(interfaces_key)
    interface_store.read_input_file(interfaces_key, interfaces_file)

    interfaces = pd.read_hdf(interfaces_file, "table")

    RealtimeLogger.info("COLS: {}".format(interfaces.columns))
    counts = interfaces.fillna(-1.).groupby(["mol_superfam_id","int_superfam_id"]).size().reset_index(name="count")

    RealtimeLogger.info("SIZES :{}".format(counts))

    try:
        os.remove(interfaces_file)
    except OSError:
        pass

    return counts

def start_toil(job, memory="1G"):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-ddi")
    if not store.exists("sorted_sfams.json"):
        interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces")
        observed = set(key.split("/")[0] for key in interface_store.list_input_directory() \
            if key.endswith("observed_interactome"))
        #inferred = set(key.split("/")[0] for key in keys if key.endswith("inferred_interactome"))

        interface_counts = [job.addChildJobFn(get_sfam_ddi_sizes, o).rv() for o in observed]
        merge_job = job.addFollowOnJobFn(best_sfams, interface_counts)
    else:
        merge_job = job

    merge_job.addFollowOnJobFn(process_sfams)
    #map_job(job, process_sfam, observed, 5000)
    #map_job(job, process_sfam, inferred, True)

if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"

    detector = logging.StreamHandler()
    logging.getLogger().addHandler(detector)

    with Toil(options) as workflow:
        print "STARTING"
        workflow.start(Job.wrapJobFn(start_toil))
