from __future__ import print_function
import os
import sys
import time
import logging

import numpy as np
import pandas as pd

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from botocore.exceptions import ClientError

from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.iRMSD import Complex
from molmimic.generate_data.prepare_protein import copy_pdb_h5, process_domain
from molmimic.generate_data.util import filter_hdf, get_all_chains
from molmimic.parsers.Electrostatics import run_pdb2pqr
from molmimic.generate_data.job_utils import map_job
from molmimic.generate_data.BioUnit import check_contacts, build_biounit

from six.moves import zip

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
        pdb.upper(), chain, int(sdi), int(domNo))
    file = os.path.join(work_dir, os.path.basename(prefix))
    if extracted:
        prefix += ".extracted"
    assert in_store.exists(prefix)
    in_store.read_input_file(prefix, file)
    return file

def getOverlap(a0, a1, b0, b1):
    overlap = max(0, min(a1, b1) - max(a0, b0))
    overlap /= float(min(a1-a0, b1-b0))
    return overlap

def get_int_file(job, mol_pdb, int_pdb, mol_chain, int_chain,
  int_sdi_id, int_domNo, mol_sdi_from, mol_sdi_to, int_sdi_from,
  int_sdi_to, mol_superfam_id, int_superfam_id, pdbFileStoreID,
  download, work_dir, nbr=None):

    def to_int(s):
        try:
            return int(s)
        except ValueError:
            return int(float(s))

    int_inside = False
    int_file = None

    if mol_pdb == int_pdb and mol_chain == int_chain:
        int_inside = True

        if nbr is not None:
            obs_file = "{0}.observed_interactome".format(mol_superfam_id)
            if not os.path.isfile(obs_file):
                obs_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")
                try:
                    obs_key = "{0}/{0}.observed_interactome".format(int(mol_superfam_id))
                    obs_store.read_input_file(obs_key, obs_file)
                except ClientError:
                    return None, None
            obs_df = filter_hdf(obs_file, "table", "obs_int_id", nbr)
            obs_df = obs_df.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
                as_index=False).nth(0).reset_index(drop=True).copy()

            mol_sdi_from = map(to_int, obs_df.mol_sdi_from.tolist())
            mol_sdi_to = map(to_int, obs_df.mol_sdi_to.tolist())
            int_sdi_from = map(to_int, obs_df.int_sdi_from.tolist())
            int_sdi_to = map(to_int, obs_df.int_sdi_to.tolist())

            del obs_df
        else:
            mol_sdi_from = map(to_int, mol_sdi_from.split(","))
            mol_sdi_to = map(to_int, mol_sdi_to.split(","))
            int_sdi_from = map(to_int, int_sdi_from.split(","))
            int_sdi_to = map(to_int, int_sdi_to.split(","))


        mol_sdi_ranges = list(zip(mol_sdi_from, mol_sdi_to))
        int_sdi_ranges = list(zip(int_sdi_from, int_sdi_to))

        for mol_from, mol_to in mol_sdi_ranges:
            need_match = False
            found_match = False
            for int_from, int_to in int_sdi_ranges:
                if mol_from<=int_to and int_from<=mol_to:
                    #They overlap, no int is inside mol => maybe wrong mol_id
                    need_match = True
                    pdb_info_file = copy_pdb_h5(job, pdbFileStoreID)
                    pdb = filter_hdf(str(pdb_info_file), "merged", "pdbId", mol_pdb)
                    if pdb[pdb["sdi"]==int_sdi_id].iloc[0]["whole_chn"] == 1.0:
                        chains = pdb.groupby("chnLett", as_index=False)
                        mol_chain_df = chains.get_group(mol_chain)
                        mol_chain_ranges = list(zip(mol_chain_df["from"], mol_chain_df["to"]))

                        if chains.ngroups > 2:
                            #Not sure how to handle this
                            match = False
                            break

                        for chnLett, chain in chains:
                            if chnLett == mol_chain: continue
                            chain_ranges = list(zip(chain["from"], chain["to"]))

                            match = True
                            for a0, a1 in mol_chain_ranges:
                                for b0, b1 in chain_ranges:
                                    if getOverlap(a0, a1, b0, b1) >= 0.9:
                                        break
                                else:
                                    match = False
                                    break

                            if match:
                                if download:
                                    try:
                                        int_file = download_pdb(job, int_superfam_id, int_pdb,
                                            chnLett, int_sdi_id, int_domNo, work_dir=work_dir)
                                    except (ClientError, AssertionError):
                                        #Multiple ranges for one sdi. Might be domain swaped?
                                        rslices = ["{}:{}".format(chain["from"].min(), chain["to"].max())]
                                        try:
                                            int_file, _, _ = process_domain(job,
                                                int_sdi_id,
                                                pdbFileStoreID,
                                                force_chain=chnLett,
                                                force_rslices=rslices,
                                                work_dir=work_dir)
                                        except ClientError:
                                            continue
                                else:
                                    int_file = -1

                                found_match = True
                                break
                            else:
                                continue

                    if found_match:
                        break

                if found_match:
                    break

            if found_match:
                break

        if not need_match and int_file is None:
            #Regular non overlapping domian domian interaction
            if download:
                try:
                    int_file = download_pdb(job, int_superfam_id, int_pdb,
                        int_chain, int_sdi_id, int_domNo, work_dir=work_dir)
                except (ClientError, AssertionError):
                    RealtimeLogger.info("Could not download int {} {} {} {}".format(
                        int_superfam_id, int_pdb,
                        int_chain, int_sdi_id))
                    return None, None
            else:
                int_file = -1

        elif not need_match:
            RealtimeLogger.info("No match for int {} {} {} {}".format(
                int_superfam_id, int_pdb,
                int_chain, int_sdi_id))
            return None, None

    else:
        #Normal, get regualr interacting pdb_file
        if download:
            try:
                int_file = download_pdb(job, int_superfam_id, int_pdb,
                    int_chain, int_sdi_id, int_domNo, work_dir=work_dir)
            except (ClientError, AssertionError):
                RealtimeLogger.info("Could not download int {} {} {} {}".format(
                    int_superfam_id, int_pdb,
                    int_chain, int_sdi_id))
                return None, None
        else:
            int_file = -1

    return int_file, int_inside

def check_full_chain(job, sdi, pdb, res, pdbFileStoreID, work_dir):
    mol = Complex.parser.get_structure("ref", pdb)
    face = [r.id[1] for r in mol.get_residues() if r.id[1] in res]
    del mol
    if len(face)==0:
        try:
            pdbf, _, _ = process_domain(job,
                sdi,
                pdbFileStoreID,
                force=True,
                cleanup=False,
                work_dir=work_dir)
            RealtimeLogger.info("NEW FILE IS {}".format(pdbf))
            assert os.path.isfile(pdbf)
        except (RuntimeError, AssertionError):
            return None
    else:
        pdbf = pdb

    return pdbf

def process_interface(job, row, int_type, pdbFileStoreID, download=True, work_dir=None):
    from molmimic.parsers.superpose import align

    if work_dir is None:
        work_dir = os.getcwd()

    if int_type=="observed":
        mol_superfam_id = row.mol_superfam_id
        mol_sdi_id = row.mol_sdi_id
        mol_domNo = row.mol_domNo
        int_sdi_id = row.int_sdi_id
        int_domNo = row.int_domNo
    else:
        mol_superfam_id = row.nbr_superfam_id
        mol_sdi_id = int(row.mol_sdi)
        mol_domNo = int(row.mol_domNo)
        int_sdi_id = int(row.int_sdi)
        int_domNo = int(row.int_domNo)


    if download:
        try:
            mol_file = download_pdb(job, mol_superfam_id, row.mol_pdb,
                row.mol_chain, mol_sdi_id, mol_domNo, work_dir=work_dir)
        except (ClientError, AssertionError):
            try:
                #Multiple ranges for one sdi. Might be domain swaped?
                mol_file, _, _ = process_domain(job,
                    mol_sdi_id,
                    pdbFileStoreID,
                    work_dir=work_dir)
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                raise
                RealtimeLogger.info("Could not download mol {} {} {} {}".format(
                    mol_superfam_id, row.mol_pdb,
                    row.mol_chain, mol_sdi_id))
                return None, None, None, None, False, {}
    else:
        mol_file = None

    if int_type=="observed":
        int_file, int_inside = get_int_file(job, row.mol_pdb, row.int_pdb,
            row.mol_chain, row.int_chain, row.int_sdi_id,
            row.int_domNo, row.mol_sdi_from, row.mol_sdi_to,
            row.int_sdi_from, row.int_sdi_to, row.mol_superfam_id,
            row.int_superfam_id, pdbFileStoreID, download, work_dir)
    else:
        int_file, int_inside = get_int_file(job, row.nbr_pdb, row.int_pdb,
            row.nbr_chain, row.int_chain, row.int_sdi,
            row.int_domNo, row.mol_sdi_from, row.mol_sdi_to,
            row.int_sdi_from, row.int_sdi_to, row.nbr_superfam_id,
            row.int_superfam_id, pdbFileStoreID, download, work_dir,
            nbr=row.nbr_obs_int_id)

    if int_file is None:
        RealtimeLogger.info("Could not download int {} {} {} {}".format(
            row.int_superfam_id, row.int_pdb,
            row.int_chain, int_sdi_id))
        return None, None, None, None, False, {}

    #Set chain names as M=mol, I=int
    if download:
        mol_file, int_file = list(prep((mol_file, "M"), (int_file, "I"), work_dir=work_dir))

    try:
        mol_resi = list(map(int, str(row["mol_res"]).split(",")))
    except ValueError:
        if isinstance(row["mol_res"], np.float64):
            mol_resi = [int(row["mol_res"])]
        else:
            RealtimeLogger.info("Could not convert mol resi {} {} {} {} {}".format(
                row["mol_res"],
                row.mol_superfam_id, row.mol_pdb,
                row.mol_chain, mol_sdi_id))
            return None, None, None, None, False, {}

    mol_file = check_full_chain(job, mol_sdi_id, mol_file, mol_resi, pdbFileStoreID, work_dir)
    if mol_file is None:
        return None, None, None, None, False, {}

    try:
        int_resi = list(map(int, str(row["int_res"]).split(",")))
    except ValueError:
        if isinstance(row["int_res"], np.float64):
            int_resi = [int(row["int_res"])]
        else:
            RealtimeLogger.info("Could not convert int resi {} {} {} {} {}".format(
                row["int_res"],
                row.int_superfam_id, row.int_pdb,
                row.int_chain, int_sdi_id))
            return None, None, None, None, False, {}

    int_file = check_full_chain(job, int_sdi_id, int_file, int_resi, pdbFileStoreID, work_dir)
    if int_file is None:
        return None, None, None, None, False, {}

    if int_type == "inferred":
        try:
            nbr_file = download_pdb(job, row.nbr_superfam_id, row.nbr_pdb,
                row.nbr_chain, int(row.nbr_sdi_id), int(row.nbr_domNo))
        except (ClientError, AssertionError):
            RealtimeLogger.info("Could not download nbr {} {} {} {}".format(
                row.nbr_superfam_id, row.nbr_pdb,
                row.nbr_chain, row.nbr_sdi_id))
            return None, None, None, None, int_inside, {}

        #Set chain name: N=nbr
        nbr_file = next(prep((nbr_file, "N")))

        #Reset mol_file with it superposed on top of nbr
        mol_file, rmsd, tm_score, _ = align(nbr_file, "N",
            mol_file, "M", work_dir=work_dir, job=job)
        align_info = {"nbr_rmsd":rmsd, "nbr_tm-score":tm_score}
    else:
        align_info = {}

    return mol_file, mol_resi, int_file, int_resi, int_inside, align_info

def advanced_dock(job, mol_file, mol_resi, int_file, int_resi, original_complex, int_id=None, rerun=False, work_dir=None, cores=1):
    from molmimic.parsers.haddock import dock

    if work_dir is None:
        work_dir = os.getcwd()

    #dock_name = row["obs_int_id"] if int_type=="observed" else "fix_me"

    #mol_file, mol_resi, int_file, int_resi = process_interface(job, row, int_type, work_dir=work_dir)



    # #Perform docking
    # docked_file = dock(dock_name, mol_file, "M", inferred_row.mol_resi,
    #     int_file, "I", inferred_row.int_resi, refine=True, small_refine=True, job=job)
    #
    # #Analyze docking, from PRODIGY
    # original_complex = next(prep((mol_file, "M"), (int_file, "I"), merge=True))
    # results = analyze(job, docked_file, original_complex)
    # return docked_file, results

    t0 = time.time()

    if not rerun:
        #Refine:HADDOCK #Score:Haddock
        fname = "{}__{}".format(
            os.path.splitext(os.path.basename(mol_file))[0],
            os.path.splitext(os.path.basename(int_file))[0])
        dock_name = "{}_{}".format(int_id if int_id is not None else "inf", fname)

        orig_file, complex_file, haddock_result, haddock_files = dock(
            dock_name,
            mol_file, "M", mol_resi,
            int_file, "I", int_resi,
            small_refine=True,
            work_dir=work_dir,
            clean_docked_file=False,
            cores=cores,
            job=job)
    else:
        complex_file = mol_file
        haddock_files = None

    refined_complex = Complex(complex_file, face1=mol_resi, face2=int_resi,
        method="advanced", work_dir=work_dir, job=job)
    if not rerun:
        refined_complex.results.append(haddock_result)
    refined_complex.compare_to_reference(original=original_complex)

    if not rerun:
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

def simple_dock(job, orig_file, mol_resi, int_resi, original_complex, rerun=False, work_dir=None):
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
    if not rerun:
        complex_file, cns_results = Minimize(orig_file, work_dir=work_dir, job=job)
    else:
        complex_file = orig_file

    refined_complex = Complex(complex_file, face1=mol_resi, face2=int_resi,
        method="simple", work_dir=work_dir, job=job)
    if not rerun:
        refined_complex.results.append(cns_results)
    refined_complex.compare_to_reference(original=original_complex)

    if not rerun:
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

def get_original_complexes(job, mol_sfam, int_sfam, group, int_type, pdbFileStoreID, work_dir=None):
    if work_dir is None:
        work_dir = os.cwd()

    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")

    files_to_clean = []
    complex_files = []
    interface_info = []
    for _, row in group:
        orow = row.copy()
        row = row.iloc[0]

        if int_type == "observed":
            int_id = row["obs_int_id"]
        else:
            int_id = "{}_{}_{}".format(row["nbr_obs_int_id"], row["nbr_sdi_id"], row["mol_sdi"])

        unrefined_files = list(ddi_store.list_input_directory("{}_{}/unrefined_{}/{}".format(
            int(mol_sfam), int(int_sfam), int_type, int_id)))

        if len(unrefined_files) == 0:
            #Put together interaction from monomers

            mol_file_f, mol_resi, _int_file_f, int_resi, int_inside, _align_info = process_interface(job,
                row, int_type, pdbFileStoreID, work_dir=work_dir)

            files_to_clean += [mol_file_f, _int_file_f]

            if mol_file_f is None:
                #PDB files not found, skip
                RealtimeLogger.info("Cannot download mol PDB {}.{}.{} bc it was not found {}".format(
                    row.mol_pdb, row.mol_chain, row.mol_sdi_id if int_type=="observed" else row.mol_sdi, "INT INSIDE" if int_inside else ""))
                if int_inside:
                    complex_files.append((-1, None))
                else:
                    complex_files.append((None, None))
                continue
            if _int_file_f is None:
                RealtimeLogger.info("Cannot download int PDB {}.{}.{} bc it was not found {}".format(
                    row.int_pdb, row.int_chain, row.int_sdi_id if int_type=="observed" else row.int_sdi, "INT INSIDE" if int_inside else ""))
                if int_inside:
                    complex_files.append((-1, None))
                else:
                    complex_files.append((None, None))
                continue

            #Check for biounit
            try:
                for i, int_file_f in enumerate(build_biounit(
                  mol_file_f, row.mol_pdb, row.mol_chain, _int_file_f, row.int_pdb,
                  row.int_chain, work_dir)):
                    RealtimeLogger.info("Trying Biounit {}".format(i))
                    merged_file_f = next(prep((mol_file_f, "M"), (int_file_f, "I"), merge=True,
                        work_dir=work_dir))
                    files_to_clean += [int_file_f, merged_file_f]

                    if check_contacts(merged_file_f, mol_file_f, mol_resi, int_file_f, int_resi):
                        break
                else:
                    RealtimeLogger.info("Failed to find complex for {} {}: {}".format(mol_file_f, _int_file_f,
                        check_contacts(merged_file_f, mol_file_f, mol_resi, int_file_f, int_resi, return_vars=True)))
                    complex_files.append((None, None))
                    continue
            except RuntimeError:
                RealtimeLogger.info("Failed (RuntimeError) to find complex for {} {}".format(mol_file_f, _int_file_f))
                complex_files.append((None, None))
                continue

        else:
            #Unrefined DDI exists
            merged_file_f = os.path.join(work_dir, os.path.basename(unrefined_files[0]))
            files_to_clean.append(merged_file_f)
            try:
                ddi_store.read_input_file(unrefined_files[0], merged_file_f)

                _, mol_resi, _, int_resi, int_inside, _align_info = process_interface(job,
                    row, int_type, pdbFileStoreID, download=False, work_dir=work_dir)

                #Split chains
                subprocess.call([sys.executable, os.path.join(PDB_TOOLS, "pdb_splitchain.py"), merged_file_f])
                fname_root = merged_file_f[:-4]
                mol_file_f = fname_root+"_M.pdb"
                int_file_f = fname_root+"_I.pdb"
                assert os.path.isfile(mol_file_f)
                assert os.path.isfile(int_file_f)
                files_to_clean += [mol_file_f, int_file_f]
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                complex_files.append((None, None))
                continue


        mol_file = (os.path.basename(mol_file_f), job.fileStore.writeGlobalFile(mol_file_f, cleanup=True))
        int_file = (os.path.basename(int_file_f), job.fileStore.writeGlobalFile(int_file_f, cleanup=True))
        merged_file = (os.path.basename(merged_file_f), job.fileStore.writeGlobalFile(merged_file_f, cleanup=True))

        complex_files.append(merged_file)

        interface_info.append((mol_file, mol_resi, int_file, int_resi, int_inside, _align_info))

        if int_type == "observed":
            int_id = row["obs_int_id"]
        else:
            int_id = "{}_{}_{}".format(row["nbr_obs_int_id"], row["nbr_sdi_id"], row["mol_sdi"])

        ddi_store.write_output_file(merged_file_f, "{}_{}/unrefined_{}/{}_{}".format(
            int(mol_sfam), int(int_sfam), int_type, int_id, merged_file[0]))

        for f in files_to_clean:
            if not f: continue
            try:
                os.remove(f)
            except OSError:
                pass

    return complex_files, interface_info

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
    log_file = os.path.join(cluster_dir, "{}_{}.max_cluster_logs".format(int(mol_sfam), int(int_sfam)))
    distance_file = os.path.join(cluster_dir, "{}_{}_distances.txt".format(int(mol_sfam), int(int_sfam)))

    ddi_base = "{}_{}".format(int(mol_sfam), int(int_sfam))
    file_key = "{}/{}".format(ddi_base, os.path.basename(file_list))
    log_key = "{}/{}".format(ddi_base, os.path.basename(log_file))
    distance_key = "{}/{}".format(ddi_base, os.path.basename(distance_file))

    successful_pdbs = []

    with open(file_list, "w") as f:
        for i, (pdb_file_name, pdb_file_store) in enumerate(pdb_complexes):
            if not pdb_file_store: continue
            pdb_file = job.fileStore.readGlobalFile(pdb_file_store, os.path.join(work_dir, pdb_file_name))
            cmds = [
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-A", pdb_file],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_reres.py")],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_reatom.py")],
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_tidy.py")]
            ]

            pdb_to_cluster = os.path.join(cluster_dir, os.path.basename(pdb_file))
            with open(pdb_to_cluster, "w") as pdb:
                SubprocessChain(cmds, pdb)

            print(os.path.basename(pdb_to_cluster), file=f)
            successful_pdbs.append((i, os.path.basename(pdb_to_cluster)))

    ddi_store.write_output_file(file_list, file_key)

    if ddi_store.exists(log_key):
        try:
            ddi_store.read_input_file(log_key, log_file)
            ddi_store.read_input_file(distance_key, distance_file)
        except ClientError:
            return None, None
        logs = log_file
    else:
        if len(successful_pdbs)>1:
            logs = run_maxcluster("rmsd", file_list=file_list, log=True, R=distance_file,
                work_dir=cluster_dir, C=1, P=10, job=job)
        else:
            return successful_pdbs[0][0], successful_pdbs[0][1]

        if logs is None:
            return None, None

        ddi_store.write_output_file(logs, log_key)
        ddi_store.write_output_file(distance_file, distance_key)

    centroid = get_centroid(logs, job=job)
    centroid_index = None
    for i, p in successful_pdbs:
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

def combine_from_to(df):
    row = df.iloc[0].copy()
    row["mol_sdi_from"] = ",".join(df["mol_sdi_from"].astype(str))
    row["mol_sdi_to"] = ",".join(df["mol_sdi_to"].astype(str))
    row["int_sdi_from"] = ",".join(df["int_sdi_from"].astype(str))
    row["int_sdi_to"] = ",".join(df["int_sdi_to"].astype(str))
    return row

def process_sfam_sfam(job, group, pdbFileStoreID, int_type, rerun=False, cores=1):
    RealtimeLogger.info("STARTING sfam sfam {}".format(group))
    work_dir = job.fileStore.getLocalTempDir()
    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")

    if isinstance(group, (list, tuple)):
        mol_sfam, int_sfam = group[0]
        group = group[1]
    else:
        mol_sfam, int_sfam = group.iloc[0]["mol_superfam_id"], group.iloc[0]["int_superfam_id"]

    if int_type == "observed":
        #Get PPI Bind Structures
        pdb_bind_f = os.path.join(work_dir, "pdb_bind_ppi.csv")
        try:
            ddi_store.read_input_file("pdb_bind_ppi.csv", pdb_bind_f)
            pdb_bind = pd.read_csv(pdb_bind_f, usecols=["pdbId"])["pdbId"]
        except ClientError:
            pdb_bind = pd.Series([], name="pdbId")

        in_pdbbind = group["mol_pdb"].apply(lambda x: x in pdb_bind)
        in_pdbbind.index = group.index
        group = group.assign(in_pdbbind=in_pdbbind)

    ddi_base = "{}_{}".format(int(mol_sfam), int(int_sfam))

    #finished = [k for k in ddi_store.list_input_directory(ddi_base) if k.endswith("results.csv")]
    #if not rerun and any(finished):
    #    continue

    if int_sfam == -1.:
        #Only keep groups where known sfams are specified
        return

    group[["mol_sdi_from", "mol_sdi_to"]] = group[["mol_sdi_from", "mol_sdi_to"]].astype(str)

    if int_type == "observed":
        #Remove redundant interfaces
        group = group.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
            as_index=False).nth(0).reset_index(drop=True).copy()
    else:
        group = group.groupby(["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi", "mol_sdi_from", "mol_sdi_to"],
            as_index=False).nth(0).reset_index(drop=True).copy()

    if int_type == "observed":
        #Merge similar obs ints
        group = group.groupby(["obs_int_id"], as_index=False).apply(combine_from_to)
        if "mol_sdi" in group.columns:
            key = "mol_sdi"
        elif "mol_sdi_id" in group.columns:
            key = "mol_sdi_id"
        else:
            raise RuntimeError("sdi not in df")
        sdi_group = group.groupby(key, as_index=False)
    else:
        #Merge similar inf ints
        group = group.groupby(["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi"],
            as_index=False).apply(combine_from_to).reset_index(drop=True)
        sdi_group = group.groupby(["nbr_obs_int_id", "nbr_sdi_id", "mol_sdi"],
            as_index=False)

    if rerun:
        results_files = [k for k in ddi_store.list_input_directory(ddi_base) \
            if k.endswith("results.csv")]
        if len(results_files) == 0:
            return
        results_file = results_files[0]
        results_dir = os.path.dirname(results_file)
        files = [os.path.join(results_dir, "original.pdb")]
        orig_file = files[0]
        advanced_file = os.path.join(results_dir, "haddock.pdb")
        simple_file = os.path.join(results_dir, "cns.pdb")
        best_sdi = int(os.path.basename(results_dir) )
        sdi_group = [(n,g) for n, g in sdi_group if g.iloc[0].obs_int_id == best_sdi]

        for rk in [results_file, orig_file, advanced_file, simple_file]:
            ddi_store.read_input_file(rk, os.path.join(work_dir, os.path.basename(rk)))

        orig_file = os.path.join(work_dir, os.path.basename(orig_file))
        advanced_file = os.path.join(work_dir, os.path.basename(advanced_file))
        simple_file = os.path.join(work_dir, os.path.basename(simple_file))
        results_file = os.path.join(work_dir, os.path.basename(results_file))

    else:
        files, interface_info = get_original_complexes(job, mol_sfam, int_sfam, sdi_group,
            int_type, pdbFileStoreID, work_dir=work_dir)

    if len(files) == 0 or files.count((None, None))+files.count((-1, None)) == len(files):
        RealtimeLogger.info("No PDBs for {}:{}".format(mol_sfam, int_sfam))
        fail_file = os.path.join(work_dir, "fail")
        with open(fail_file, "w") as f:
            print("No PDBs for {}:{}".format(mol_sfam, int_sfam), file=f)
            if files.count(-1) > 0:
                for file, sdi in zip(files, sdi_group):
                    if not file==-1: continue
                    print("Failed finding match for alternate chain in {}".format(
                        sdi), file=f)
        ddi_store.write_output_file(fail_file, "{}/failed".format(ddi_base))
        try:
            os.remove(fail_file)
        except OSError:
            pass
        return
    elif len(files) == 1:
        centroid_index, centroid = 0, files[0]
        if not isinstance(sdi_group, (list, tuple)):
            if isinstance(sdi_group, pd.core.groupby.generic.DataFrameGroupBy):
                if sdi_group.ngroups > 0:
                    sdi_group = [(list(sdi_group.groups.keys())[0], sdi_group.first())]
                else:
                    #Error getting group, skip
                    return
            else:
                try:
                    sdi_group = [next(sdi_group)]
                except TypeError:
                    #Error getting group, skip
                    return
    else:
        if int_type == "observed":
            centroid_index, centroid = cluster(job, mol_sfam, int_sfam, files, work_dir=work_dir)
            if centroid_index is None and centroid is None:
                return
        else:
            centroid = centroid_index = sdi_group["nbr_score"].max().idxmax()

    def ddi_exists(sdi_key, row):
        if int_type == "observed":
            return ddi_store.exists("{}/{}/results.csv".format(ddi_base, row.iloc[0]["obs_int_id"]))
        else:
            int_id = "_".join(map(lambda x: str(int(float(x))), sdi_key))
            ddi_key = "{}/inferred/{}/results.csv".format(ddi_base, int_id)
            return ddi_store.exists(ddi_key)

    groups = [grouping for grouping in enumerate(zip(sdi_group, files, interface_info)) \
        if not ddi_exists(*grouping[1][0]) and \
        grouping[1][1][1] is not None and \
        (grouping[1][0][1].iloc[0].in_pdbbind or \
            (int_type == "observed" and grouping[0]==centroid_index) or \
            (int_type == "inferred" and row.iloc[0]["nbr_score"]==centroid))
        ]

    RealtimeLogger.info("RUNNING groups {} {}".format(len(groups), groups))



    map_job(job, run_refinement, groups, centroid_index, centroid, ddi_base, int_type, rerun)
        # i = grouping[0]
        # group_info = grouping[1][0]
        # file = grouping[1][1]
        #
        # sdi_key, row = group_info
        # mol_file, mol_resi, int_file, int_resi, int_inside, align_info = interface_info[i]
            # process_interface(job, row.iloc[0], int_type, pdbFileStoreID,
            #     download=not rerun, work_dir=work_dir)
        # if not rerun:
        #     orig_file = file #next(prep((mol_file, "M"), (int_file, "I"), merge=True, work_dir=work_dir))
        # else:
        #     mol_file = advanced_file

def flatten(l):
    if isinstance(l, (list, tuple)):
        return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else [])
    return [l]

def run_refinement(job, run_params, centroid_index, centroid, ddi_base, int_type, rerun, cores=1):
    work_dir = job.fileStore.getLocalTempDir()
    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")

    i = run_params[0]
    sdi_key, row = run_params[1][0]
    orig_file = job.fileStore.readGlobalFile(
        run_params[1][1][1], os.path.join(work_dir, run_params[1][1][0])) #run_params[1][1]
    mol_file, mol_resi, int_file, int_resi, int_inside, align_info = run_params[1][2]

    mol_file = job.fileStore.readGlobalFile(mol_file[1], os.path.join(work_dir, mol_file[0]))
    int_file = job.fileStore.readGlobalFile(int_file[1], os.path.join(work_dir, int_file[0]))

    # flattend_params = flatten(run_params)
    # RealtimeLogger.info("RUN PARAMS {}".format(run_params))
    # RealtimeLogger.info("RUN PARAMS {}".format(len(run_params)))
    # RealtimeLogger.info("FLATTEDN PARAMS {}".format(flattend_params))
    # RealtimeLogger.info("FLATTEDN PARAMS {} 9".format(len(flattend_params), ))
    # sdi_key, row, orig_file, mol_file, mol_resi, int_file, int_resi, int_inside, align_info = \
    #     flatten(run_params)
    # orig_file = file
    # RealtimeLogger.info("MOL RESI {}".format(mol_resi))
    # RealtimeLogger.info("INT RESI {}".format(int_resi))
    if int_type == "observed":
        ddi_key = "{}/{}".format(ddi_base, row.iloc[0]["obs_int_id"])
        all_results_path = os.path.join(work_dir, "{}.csv".format(row.iloc[0]["obs_int_id"]))
    else:
        int_id ="_".join(map(lambda x: str(int(float(x))), sdi_key))
        ddi_key = "{}/inferred/{}".format(ddi_base, int_id)
        all_results_path = os.path.join(work_dir, "{}.csv".format(int_id))

    original_complex = Complex(orig_file, face1=mol_resi, face2=int_resi,
        method="original", work_dir=work_dir, job=job)

    ddi_store.write_output_file(orig_file, "{}/original.pdb".format(ddi_key))
    #original_complex.set_prefix().to_frame().T.to_csv(all_results_path+".tmp")
    #ddi_store.write_output_file(all_results_path+".tmp", "{}/results.csv.tmp".format(ddi_key))

    if int_type == "observed":
        if not rerun and ddi_store.exists("{}/results.csv".format(ddi_key)):
            return

        advanced_complex, adv_files = advanced_dock(job, mol_file,
            mol_resi, int_file, int_resi, original_complex,
            int_id=row.iloc[0]["obs_int_id"], rerun=rerun,
            work_dir=work_dir, cores=cores)

        ddi_store.write_output_file(adv_files, "{}/haddock_info.tgz".format(ddi_key))
        ddi_store.write_output_file(advanced_complex.pdb, "{}/haddock.pdb".format(ddi_key))

        simple_complex = simple_dock(job, orig_file, mol_resi, int_resi,
            original_complex, rerun=rerun, work_dir=work_dir)

        ddi_store.write_output_file(simple_complex.pdb, "{}/cns.pdb".format(ddi_key))

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

        if rerun:
            old_results = pd.read_csv(results_file, index_col=0)
            old_results.update(all_results)
            all_results = old_results

        all_results.to_csv(all_results_path)

        if not rerun:
            ddi_store.write_output_file(adv_files, "{}/haddock_info.tgz".format(ddi_key))

        ddi_store.write_output_file(advanced_complex.pdb, "{}/haddock.pdb".format(ddi_key))


        ddi_store.write_output_file(all_results_path, "{}/results.csv".format(ddi_key))
    if int_type == "inferred":
        redock_inferred(job, sdi_key, row, file, centroid,
            mol_file, mol_resi, int_file, int_resi,
            original_complex, align_info, ddi_base, work_dir,
            rerun=rerun)

def redock_inferred(job, sdi_key, row, file, centroid,
  mol_file, mol_resi, int_file, int_resi, original_complex,
  align_info, ddi_base, work_dir, advanced_inf_dock=False, rerun=False):
    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")
    int_id ="_".join(map(lambda x: str(int(float(x))), sdi_key))
    ddi_key = "{}/inferred/{}".format(ddi_base, int_id)

    for k,v in align_info.items():
        original_complex[k] = v

    simple_complex = simple_dock(job, original_complex.pdb, mol_resi, int_resi,
        original_complex, rerun=rerun, work_dir=work_dir)

    if advanced_inf_dock and sdi_key==centroid:
        run_adv = True
        advanced_complex, adv_files = advanced_dock(job, mol_file,
            mol_resi, int_file, int_resi, original_complex,
            int_id=int_id, rerun=rerun,
            work_dir=work_dir, cores=cores)
        simple_complex.compare_to_reference(adv_ref=advanced_complex)
        advanced_complex.compare_to_reference(simple_ref=simple_complex)

        all_results = pd.concat([
            original_complex.set_prefix(),
            simple_complex.set_prefix(),
            advanced_complex.set_prefix()]).to_frame().T
    else:
        run_adv = False
        all_results = pd.concat([
            original_complex.set_prefix(),
            simple_complex.set_prefix()]).to_frame().T

    if rerun:
        old_results = pd.read_csv(results_file, index_col=0)
        old_results.update(all_results)
        all_results = old_results

    all_results_path = os.path.join(work_dir, "{}.csv".format(int_id))
    all_results.to_csv(all_results_path)

    if not rerun:
        ddi_store.write_output_file(original_complex.pdb, "{}/original.pdb".format(ddi_key))

    if not rerun and run_adv:
        ddi_store.write_output_file(adv_files, "{}/haddock_info.tgz".format(ddi_key))

    if run_adv:
        ddi_store.write_output_file(advanced_complex.pdb, "{}/haddock.pdb".format(ddi_key))

    ddi_store.write_output_file(simple_complex.pdb, "{}/cns.pdb".format(ddi_key))
    ddi_store.write_output_file(all_results_path, "{}/results.csv".format(ddi_key))

def process_sfam(job, sfam_id, pdbFileStoreID, observed=True, int_sfams=None, min_cluster=0, rerun=False, cores=1):
    int_type = "observed" if observed else "inferred"
    work_dir = job.fileStore.getLocalTempDir()
    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")

    if isinstance(sfam_id, (list, tuple)) and len(sfam_id) == 2:
        sfam_id, int_sfams = sfam_id[0], sfam_id[1]

    interfaces_key = "{s}/{s}.{o}_interactome".format(
        s=sfam_id, o="observed" if observed else "inferred")
    interfaces_file = os.path.basename(interfaces_key)

    try:
        interface_store.read_input_file(interfaces_key, interfaces_file)
    except ClientError:
        RealtimeLogger.info("Cannot read SFAM: {}".format(sfam_id))
        return

    interfaces = pd.read_hdf(interfaces_file, "table")

    if int_sfams is not None:
        interfaces = interfaces[interfaces["int_superfam_id"].isin(int_sfams)]

    if observed:
        igroups = interfaces.fillna(-1.).groupby(["mol_superfam_id", "int_superfam_id"], as_index=False)
    else:
        igroups = interfaces.fillna(-1.).groupby(["nbr_superfam_id", "int_superfam_id"], as_index=False)
    igroups = [g for g in igroups]

    map_job(job, process_sfam_sfam, igroups, pdbFileStoreID, int_type, rerun)

def process_sfams(job, pdbFileStoreID, observed=True, max_sfams=300, memory="1G"):
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
        if -1 in list(map(int, (mol_sfam, int_sfam))):
            RealtimeLogger.info("SKIPPED")
            continue

        if count<90:
            RealtimeLogger.info("SKIPPED count is less than 100")
            #continue

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

    sfams = list(sfams.items())

    map_job(job, process_sfam, sfams, pdbFileStoreID, observed)

    # for mol_sfam, int_sfams in list(sfams.items()):
    #     job.addChildJobFn(process_sfam, mol_sfam, pdbFileStoreID, int_sfams=int_sfams)

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

    sfams = sorted(iter(list(ddi_counts.items())), key=lambda x: x[1], reverse=True)
    RealtimeLogger.info("sfams is {}".format(sfams))
    sfam_file = os.path.join(work_dir, "sorted_sfams.json")
    with open(sfam_file, "w") as f:
        json.dump(sfams, f)
    out_store.write_output_file(sfam_file, "sorted_sfams.json")

    return sfams[:max_sfams]

def get_sfam_ddi_sizes(job, sfam_id, observed=True):
    int_type = "observed" if observed else "inferred"
    work_dir = job.fileStore.getLocalTempDir()
    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")

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

def start_toil(job, observed=True, run_haddock=True, memory="1G"):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-ddi")
    pdb_store = IOStore.get("aws:us-east-1:molmimic-ibis")
    ddi_store = IOStore.get("aws:us-east-1:molmimic-ddi")

    #Download PDB info
    sdoms_file = os.path.join(work_dir, "PDB.h5")
    pdb_store.read_input_file("PDB.h5", sdoms_file)

    #Add pdb info into local job store
    pdbFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

    if not store.exists("sorted_sfams.json"):
        interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces")
        observed = set(key.split("/")[0] for key in interface_store.list_input_directory() \
            if key.endswith("observed_interactome"))
        #inferred = set(key.split("/")[0] for key in keys if key.endswith("inferred_interactome"))

        interface_counts = [job.addChildJobFn(get_sfam_ddi_sizes, o).rv() for o in observed]
        merge_job = job.addFollowOnJobFn(best_sfams, interface_counts)
    else:
        merge_job = job

    #merge_job.addFollowOnJobFn(process_sfam, "72144", pdbFileStoreID, int_sfams=["281925"])
    #295086_296406
    #278516_281043
    #304940_305068
    #merge_job.addFollowOnJobFn(process_sfam, "304408", pdbFileStoreID, int_sfams=["304408"])
    #merge_job.addFollowOnJobFn(process_sfams, pdbFileStoreID, observed)

    #map_job(job, process_sfam, observed, 5000)
    #map_job(job, process_sfam, inferred, True)
    if False:
        from itertools import groupby
        import pandas as pd

        pdb_bind_f = os.path.join(work_dir, "pdb_bind_ppi.csv")
        ddi_store.read_input_file("pdb_bind_ppi.csv", pdb_bind_f)
        pdbbind_pdbIds = pd.read_csv(pdb_bind_f, usecols=["pdbId"])["pdbId"].apply(lambda c: c.upper())

        pdb = pd.read_hdf(sdoms_file, "merged")
        pdbbind_pdbs = pdb[pdb["pdbId"].isin(pdbbind_pdbIds)]
        sfams = pdbbind_pdbs["sfam_id"].drop_duplicates().dropna().apply(lambda x: str(int(x)))

        map_job(job, process_sfam, sfams, pdbFileStoreID, observed)

    interface_store = IOStore.get("aws:us-east-1:molmimic-interfaces-1")
    observed_sfam = set(key.split("/")[0] for key in interface_store.list_input_directory() \
        if key.endswith("observed_interactome"))
    done = set([key.split("/")[0] for key in  ddi_store.list_input_directory() \
        if "inferred" not in key and key.endswith("results.csv")])
    sfams_to_run = list(observed_sfam-done)
    RealtimeLogger.info("RUNNING {} SFAMS".format(len(sfams_to_run)))
    map_job(job, process_sfam, sfams_to_run, pdbFileStoreID, observed)


if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--inferred", default=False, action="store_true")
    parser.add_argument("--no-haddock", default=False, action="store_true")
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"

    detector = logging.StreamHandler()
    logging.getLogger().addHandler(detector)

    with Toil(options) as workflow:
        workflow.start(Job.wrapJobFn(start_toil, not options.inferred))
