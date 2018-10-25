from __future__ import print_function
import os, sys

import subprocess
import gzip
import tempfile
import re
import shlex
import json
import shutil
from multiprocessing.pool import ThreadPool

import pandas as pd
import numpy as np
import dask
import dask.dataframe as dd
from joblib import Parallel, delayed

import freesasa

from Bio.PDB import PDBList

from molmimc.util import data_path_prefix, get_interfaces_path, iter_cdd, iter_unique_superfams
from molmimc.generate_data.mmcif2pdb import mmcif2pdb
from molmimic.parsers.FreeSASA import run_freesasa

NUM_WORKERS = 20
dask.config.set(scheduler='multiprocessing', num_workers=NUM_WORKERS)
dask.config.set(pool=ThreadPool(NUM_WORKERS))

RAW_PDB_PATH = os.path.join(data_path_prefix, "pdb", "pdb")
PDB_TOOLS = os.path.join(os.path.dirname(data_path_prefix), "pdb-tools")

cutoffs = {
    "weak transient": (0, 1500),
    "transient": (1500, 2500),
    "permanent": (2500, float("inf"))
}

def get_pdb(pdb, chain1, chain2, sdi1, sdi2):
    prefix = "freesasa_{}_{}_{}_full.pdb".format(pdb, chain1, chain2)
    tmp_pdb = tempfile.NamedTemporaryFile(prefix=prefix, delete=False)
    tmp_pdb_path1 = tmp_pdb_path2 = tmp_pdb.name
    pdb_path = os.path.join(RAW_PDB_PATH, pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    if os.path.isfile(pdb_path):
        with tmp_pdb as tmp, gzip.open(pdb_path, 'rb') as pdbf:
            for line in pdbf:
                tmp.write(line)
    else:
        try:
            tmp_pdb_path1, tmp_pdb_path2 = list(mmcif2pdb(pdb, [chain1, chain2]))
        except IOError:
            raise IOError("Cannot load "+pdb)

            #Try obsolete
            PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=os.getcwd(), file_format="pdb")
            tmp_pdb_path = "pdb{}.ent".format(pdb.lower())

            if not os.path.isfile(tmp_pdb_path) or os.path.getsize(tmp_pdb_path) == 0:
                PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=os.getcwd(), file_format="pdb")
                raise IOError("Cannot load "+pdb)

    tmp_pdb.close()

    prefix = "freesasa_{}_{}_{}_d1d2.pdb".format(pdb, chain1, chain2)
    tmp_pdb2 = tempfile.NamedTemporaryFile(prefix=prefix, delete=False)

    with open(tmp_pdb2.name, "w") as output:
        print("REMARK 300 ORIGINAL PDB: {}".format(pdb), file=output)
        print("REMARK 300 CHAIN {}:1 SEGMENTS: {}".format(chain1, "; ".join(sdi1)), file=output)
        if chain2 is not None and sdi2 is not None:
            print("REMARK 300 CHAIN {}:2 SEGMENTS: {}".format(chain2, "; ".join(sdi2)), file=output)
        print("REMARK 300", file=output)

        splitmodel1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", tmp_pdb_path1],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        splitchain1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain1)],
            stdin=splitmodel1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        delocc1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
            stdin=splitchain1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        splitdomain1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+list(sdi1),
            stdin=delocc1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        renamechain1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-1"],
            stdin=splitdomain1.stdout,
            stdout=output,
            stderr=subprocess.PIPE)
        out1, err1 = renamechain1.communicate()

        if chain2 is not None and sdi2 is not None:
            splitmodel2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_selmodel.py"), "-1", tmp_pdb_path2],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            splitchain2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain2)],
                stdin=splitmodel2.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            delocc2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_delocc.py")],
                stdin=splitchain2.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            splitdomain2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+list(sdi2),
                stdin=delocc2.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

            #Change chain name to fix intra-chain interactions
            renamechain2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_chain.py"), "-2"],
                stdin=splitdomain2.stdout,
                stdout=output,
                stderr=subprocess.PIPE)
            out2, err2 = renamechain2.communicate()

    os.remove(tmp_pdb.name)

    return tmp_pdb2.name

# def run_freesasa(command):
#     FNULL = open(os.devnull, 'w')
#     try:
#         freesasa = subprocess.check_output(command, stderr=FNULL)
#     except subprocess.CalledProcessError:
#         raise FreeSASAException("Failed running free sasa: {}".format(command))
#     freesasa = "{"+freesasa.split("\n",1)[1].strip().replace("-nan", "NaN").replace("nan", "NaN")
#
#     try:
#         return json.loads(freesasa)
#     except:
#         raise FreeSASAException("Unable to freesasa convert to JSON: {}".format(freesasa))

def calculate_buried_surface_area(pdb_file, pdb, sdi_sel1=None, sdi_sel2=None, face1=None, face2=None, job=None):
    """Calculate the burried surface area (BSA) of a complex. Assumes the interacting
    partners are in the same file in two separate chains, named "1" and "2" respectively.

    BSA(12) = ASA(1)+ASA(2)-ASA(12)

    where, 12 is the complex and 1 is the first monomer, 2 is the second monomer

    Paramters
    ---------
    Returns
    -------
    bsa : float
        the calculated burried surface area
    ppi_type,
    c1_asa,
    c2_asa,
    complex_asa
    """
    parameters = ["--chain-groups=12+1+2"]#["freesasa", "--format=json", "--chain-groups=12+1+2"]

    if sdi_sel1 is not None and sdi_sel1 is not None:
        command.append("--select='domain1, chain 1 and resi {}'".format("+".join(sdi_sel1).replace("-", "\-")))
        command.append("--select='domain2, chain 2 and resi {}'".format("+".join(sdi_sel2).replace("-", "\-")))

    if face1 is not None:
        residues1 = face1.replace(",", "+").replace("-", "\-")
        selection1 = "binding-site1, chain 1 and resi {}".format(residues1)
        command.append("--select='{}'".format(selection1))

    if face2 is not None:
        residues2 = face2.replace(",", "+").replace("-", "\-")
        selection2 = "binding-site2, chain 2 and resi {}".format(residues2)
        command.append("--select='{}'".format(selection2))

    #command.append(pdb_file)

    try:
        freesasa = run_freesasa(pdb_file, parameters, job)
    except FreeSASAException as e:
        print(e)
        return pd.Series({
            "c1_asa": np.nan,
            "c2_asa": np.nan,
            "complex_asa": np.nan,
            "face1_asa": np.nan,
            "face2_asa": np.nan,
            "bsa": np.nan,
            "ppi_type": "unknown",
        })

    result = {
        "c1_asa": freesasa["results"][2]["structure"][0]["area"]["total"],
        "c2_asa": freesasa["results"][3]["structure"][0]["area"]["total"],
        "complex_asa": freesasa["results"][1]["structure"][0]["area"]["total"]
    }

    if sdi_sel1 is not None and sdi_sel1 is not None:
        try:
            result["c1_asa"] = freesasa["results"][2]["structure"][0]["selections"][0]["area"]
        except (IndexError, KeyError):
            #Use full chain area, it is assumed the interacting domains are split out
            pass

        try:
            result["c2_asa"] = freesasa["results"][3]["structure"][0]["selections"][0]["area"]
        except (IndexError, KeyError):
            #Use full chain area, it is assumed the interacting domains are split out
            pass

        _complex_asa = sum([sel["area"] for sel in freesasa["results"][1]["structure"][0]["selections"] \
            if "domain" in sel["name"]])

        if _complex_asa > 0.0:
            result["complex_asa"] = _complex_asa
        else:
            print("Error, used total", complex_asa, freesasa["results"][1]["structure"][0]["selections"], command)

    if face1 is not None:
        face1_sel = {s["name"]:s["area"] for s in freesasa["results"][2]["structure"][0]["selections"]}
        result["face1_asa"] = face1_sel["binding-site1"]

    if face2 is not None:
        face2_sel = {s["name"]:s["area"] for s in freesasa["results"][3]["structure"][0]["selections"]}
        result["face2_asa"] = face2_sel["binding-site2"]

    result["bsa"] = result["c1_asa"]+result["c2_asa"]-result["complex_asa"]

    for ppi_type, (low_cut, high_cut) in cutoffs.iteritems():
        if low_cut <= result["bsa"] < high_cut:
            result["ppi_type"] = ppi_type
            break
    else:
        result["ppi_type"] = "unknown"

    return pd.Series(result)

def calculate_surface_area_chain(pdb_file, pdb, sdi=None, face=None, job=None):
    """Calculate the accesable surface area (ASA) of a complex. Assumes the chain
    of interest is named "1".

    Paramters
    ---------
    Returns
    -------
    asa : float
        the calculated accesable surface area
    """
    parameters = ["--chain-groups=1"]

    if sdi is not None:
        command.append("--select=domain, resi {}".format(sdi))

    if face is not None:
        residues = face.replace(",", "+")
        selection = "binding-site, chain 1 and resi {}".format(residues)
        command.append("--select='{}'".format(selection))

    #command.append(pdb_file)

    try:
        freesasa = run_freesasa(pdb_file, parameters, job)
    except FreeSASAException as e:
        print(e)
        if face is not None:
            return np.NaN, np.NaN
        else:
            return np.NaN

    asa = freesasa["results"][1]["structure"][0]["area"]["total"]
    if sdi is not None:
        try:
            asa = freesasa["results"][1]["structure"][0]["selections"][0]["area"]
        except (IndexError, KeyError):
            #Use full area
            pass

    if face is not None:
        selections = {s["name"]:s["area"] for s in freesasa["results"][1]["structure"][0]["selections"]}
        face1_asa = selections["binding-site"]
        return asa, face1_asa
    else:
        return asa

def get_bsa(df):
    r = df.iloc[0]

    if any(r[["mol_chain", "int_chain"]].isna()):
        #Cannot process chain, so bsa is unknown
        return pd.Series({
            "mol_sdi":r.mol_sdi,
            "bsa":np.nan,
            "ppi_type":"unknown",
            "c1_asa":np.nan,
            "c2_asa":np.nan,
            "face1_asa":np.nan,
            "face2_asa":np.nan,
            "complex_asa":np.nan})

    #If interacting domain is not well defined, use the entire chain
    df[["int_sdi_from"]] = df[["int_sdi_from"]].fillna(1)
    df[["int_sdi_to"]] = df[["int_sdi_to"]].fillna(0)

    try:
        sdi1, sdi2 = zip(*[(
          "{}:{}".format(int(row.mol_sdi_from), int(row.mol_sdi_to)),
          "{}:{}".format(int(row.int_sdi_from), int(row.int_sdi_to) if row.int_sdi_to else "")) \
          for row in df.itertuples()])
    except ValueError:
        #There is a NaN in mol sdi or from to: Invalid!
        return pd.Series({
            "mol_sdi":r.mol_sdi,
            "bsa":np.nan,
            "ppi_type":"unknown",
            "c1_asa":np.nan,
            "c2_asa":np.nan,
            "face1_asa":np.nan,
            "face2_asa":np.nan,
            "complex_asa":np.nan})

    try:
        pdb_file = get_pdb(r["mol_pdb"], r["mol_chain"], r["int_chain"], sdi1, sdi2)
    except IOError:
        return pd.Series({
            "mol_sdi":r.mol_sdi,
            "bsa":np.nan,
            "ppi_type":"unknown",
            "c1_asa":np.nan,
            "c2_asa":np.nan,
            "face1_asa":np.nan,
            "face2_asa":np.nan,
            "complex_asa":np.nan})

    area = calculate_buried_surface_area(pdb_file, r["mol_pdb"], face1=r["mol_res"], face2=r["int_res"])

    try:
        os.remove(pdb_file)
    except OSError:
        pass

    return pd.Series({
        "mol_sdi":r.mol_sdi,
        "bsa":area.bsa,
        "ppi_type":area.ppi_type,
        "c1_asa":area.c1_asa,
        "c2_asa":area.c2_asa,
        "face1_asa":area.face1_asa,
        "face2_asa":area.face2_asa,
        "complex_asa":area.complex_asa})

def observed_bsa(job, dataset_name, cdd, cores=NUM_WORKERS):
    job.log("CDD {}".format(cdd))
    prefix = os.path.join(get_interfaces_path(dataset_name), "by_superfamily", str(int(cdd)), str(int(cdd)))

    # if os.path.isfile(prefix+"_bsa.h5"):
    #     store = pd.HDFStore(unicode(prefix+"_bsa.h5"))
    #     if "/observed" in store.keys():
    #         store.close()
    #         return
    #     store.close()

    cdd_interactome_path = prefix+".observed_interactome"

    cdd_interactome = pd.read_hdf(unicode(cdd_interactome_path), "table")

    if cdd_interactome.shape[0] == 0:
        job.log("CDD observed interactome is empty -- FIX!!!")
        return

    if cdd_interactome.shape[0] == 0:
        job.log("CDD observed interactome contains intra-chain PPI, skipped -- FIX!!!")
        return

    #Remove redundant interfaces
    cdd_interactome = cdd_interactome.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
        as_index=False).nth(0).reset_index(drop=True).copy()

    if "mol_sdi" in cdd_interactome:
        key = "mol_sdi"
    elif "mol_sdi_id" in cdd_interactome:
        key = "mol_sdi_id"
    else:
        raise RuntimeError("sdi not in df")

    bsa = Parallel(n_jobs=NUM_WORKERS)(delayed(get_bsa)(group) for _, group in \
        cdd_interactome.groupby(key, as_index=False))
    bsa = pd.concat(bsa, axis=1).T
    bsa[key] = bsa[key].astype(int)
    bsa = bsa.astype({
        "bsa": np.float64,
        "c1_asa": np.float64,
        "c2_asa": np.float64,
        "complex_asa": np.float64,
        "face1_asa": np.float64,
        "face2_asa": np.float64,
        "ppi_type": str})

    cdd_interactome = pd.merge(cdd_interactome, bsa, how="left", on=key)

    cdd_interactome.to_hdf(unicode(prefix+"_bsa.h5"), "observed", table=True, format='table', complevel=9, complib="bzip2")
    print(unicode(prefix+"_bsa.h5"))

def get_asa(df):
    df = df.reset_index(drop=True)
    r = df.iloc[0]

    try:
        sdi1 = ["{}:{}".format(int(row.mol_sdi_from), int(row.mol_sdi_to)) for row in df.itertuples()]
    except:
        print("Failed due to being Series? {}".format(type(df)))
        print("DF is {}".format(df))
        raise

    try:
        pdb_file = get_pdb(r["mol_pdb"], r["mol_chain"], None, sdi1, None)
    except IOError:
        return pd.Series({
            "mol_sdi":r.mol_sdi,
            "nbr_obs_int_id":r.nbr_obs_int_id,
            "c1_asa":np.nan,
            "face1_asa": np.nan,
            "bsa":np.nan,
            "complex_asa":np.nan,
            "pred_ratio":np.nan,
            "ppi_type":"unknown"})

    asa, face1_asa = calculate_surface_area_chain(
        pdb_file, r.mol_pdb, face=r.mol_resi)

    try:
        os.remove(pdb_file)
    except OSError:
        pass

    obs_interface_asa = r.face1_asa_obs+r.face2_asa
    ratio = (face1_asa+r.face2_asa)/obs_interface_asa if obs_interface_asa > 0. else 0.
    
    predicted_bsa = r.obs_bsa*ratio
    complex_asa = asa+r.c2_asa-predicted_bsa

    for ppi_type, (low_cut, high_cut) in cutoffs.iteritems():
        if low_cut <= predicted_bsa < high_cut:
            break
    else:
        ppi_type = "unknown"

    return pd.Series({
        "mol_sdi":float(r.mol_sdi),
        "nbr_obs_int_id":float(r.nbr_obs_int_id),
        "c1_asa":asa,
        "face1_asa": face1_asa,
        "bsa":predicted_bsa,
        "complex_asa":complex_asa,
        "pred_ratio":ratio,
        "ppi_type":ppi_type})

def inferred_bsa(job, dataset_name, cdd, cores=NUM_WORKERS):
    job.log("INF CDD {}".format(cdd))
    cdd_bsa_path = os.path.join(get_interfaces_path(dataset_name), "by_superfamily", str(int(cdd)), str(int(cdd)))

    if not os.path.isfile(cdd_bsa_path+"_bsa.h5"):
        job.log("observed bsa must exist")
        return

    print("Reading obs bsa")
    store = pd.HDFStore(unicode(cdd_bsa_path+"_bsa.h5"))

    # if "/inferred" in store.keys():
    #     return

    try:
        cdd_obs_bsa = store.get("/observed")
    except KeyError:
        raise RuntimeError("Must calculate observed BSAs first")

    try:
        cdd_obs_bsa = cdd_obs_bsa[["obs_int_id", "bsa", "c1_asa", "c2_asa", "face1_asa", "face2_asa", "complex_asa", "ppi_type"]]
    except KeyError:
        job.log("Failed due to column select {}".format(cdd_obs_bsa.columns))
        raise

    cdd_obs_bsa = cdd_obs_bsa.rename(columns={
        "obs_int_id": "nbr_obs_int_id",
        "bsa": "obs_bsa",
        "c1_asa": "c1_asa_obs",
        "face1_asa": "face1_asa_obs",
        "complex_asa": "complex_asa_obs",
        "ppi_type": "ppi_type_obs"
    })

    inf_interactome_path = unicode(cdd_bsa_path+".inferred_interactome")
    try:
        print("Reading  inf interactome")
        int_store = pd.HDFStore(unicode(cdd_bsa_path+".inferred_interactome"))
        if "/table" not in int_store.keys():
            return
        m = re.search("nrows->(\d+)", int_store.info())
        if not m:
            int_store.close()
            job.log("Unable to read inferred interactome")
            return

        if int(m.group(1))> 1000000:
            int_store.close()
            return inferred_bsa_dask(cdd_obs_bsa, cdd_bsa_path)
        inf_interactome = int_store.get("/table") #pd.read_hdf(unicode(cdd_bsa_path+".inferred_interactome"), "table").reset_index()
    except MemoryError:
        return inferred_bsa_dask(cdd_obs_bsa, cdd_bsa_path)

    if inf_interactome.shape[0] > 1000000:
        int_store.close()
        del inf_interactome
        return inferred_bsa_dask(cdd_obs_bsa, cdd_bsa_path)

    inf_interactome = pd.merge(inf_interactome, cdd_obs_bsa, how="left", on="nbr_obs_int_id")

    #Remove redundant interfaces
    inf_interactome = inf_interactome.groupby(["mol_sdi", "nbr_obs_int_id", "mol_sdi_from", "mol_sdi_to"],
        as_index=False).nth(0).reset_index(drop=True).copy()

    bsa = Parallel(n_jobs=NUM_WORKERS)(delayed(get_asa)(group) for _, group in \
        inf_interactome.groupby(["mol_sdi", "nbr_obs_int_id"], as_index=False))
    bsa = pd.concat(bsa, axis=1).T
    bsa = bsa.astype({
        "mol_sdi":np.float64,
        "nbr_obs_int_id":np.float64,
        "c1_asa":np.float64,
        "face1_asa": np.float64,
        "bsa":np.float64,
        "complex_asa":np.float64,
        "pred_ratio":np.float64,
        "ppi_type":str})

    inf_interactome = pd.merge(inf_interactome, bsa, how="left", on=["mol_sdi", "nbr_obs_int_id"])

    inf_interactome.to_hdf(unicode(cdd_bsa_path+"_bsa.h5"), "inferred", format='table', append=True, complevel=9, complib="bzip2")
    print(unicode(cdd_bsa_path+"_bsa.h5"))
    int_store.close()

def inferred_bsa_dask(cdd_obs_bsa, cdd_bsa_path):
    """Same method as above, but splits the dataframe across all available cores
    so a dataframe the cannot fir into memory can still be process. This method
    is slower so not the default"""

    inf_interactome = dd.read_hdf([unicode(cdd_bsa_path+".inferred_interactome")], "table")
    inf_interactome = inf_interactome.repartition(npartitions=NUM_WORKERS)
    inf_interactome = inf_interactome.merge(cdd_obs_bsa, how="left", on="nbr_obs_int_id")

    #Remove redundant interfaces
    inf_interactome = inf_interactome.groupby(["mol_sdi", "nbr_obs_int_id", "mol_sdi_from", "mol_sdi_to"])\
        .first().reset_index()

    bsa = inf_interactome.groupby(["mol_sdi", "nbr_obs_int_id"]).apply(
        get_asa,
        meta = pd.DataFrame({
            "mol_sdi":[np.float64],
            "nbr_obs_int_id":[np.float64],
            "c1_asa":[np.float64],
            "face1_asa": [np.float64],
            "bsa":[np.float64],
            "complex_asa":[np.float64],
            "pred_ratio":[np.float64],
            "ppi_type":[str]})
        )
    bsa = bsa.astype({
        "mol_sdi":np.float64,
        "nbr_obs_int_id":np.float64,
        "c1_asa":np.float64,
        "face1_asa": np.float64,
        "bsa":np.float64,
        "complex_asa":np.float64,
        "pred_ratio":np.float64,
        "ppi_type":str})

    inf_interactome = inf_interactome.merge(bsa, how="left", on=["mol_sdi", "nbr_obs_int_id"])

    inf_interactome.to_hdf(unicode(cdd_bsa_path+"_bsa.h5"), "inferred", table=True, format='table', append=True, complevel=9, complib="bzip2")
    print(unicode(cdd_bsa_path+"_bsa.h5"))

def start_toil(job, dataset_name, iteration=0, name="bsa"):
    if iteration > 1:
        return

    memory = "240000M" if iteration == 1 else None

    path = os.path.join(get_interfaces_path(dataset_name), "by_superfamily")
    for sfam_id in iter_unique_superfams(): #cdd, sfam_id in iter_cdd(use_id=True, group_superfam=True):
        sfam_path = os.path.join(path, str(int(sfam_id)), str(int(sfam_id)))
        if not os.path.isfile(sfam_path+".observed_interactome"):
            #Make sure there are defined interfaces. If observed exsits, then inferred may exist, but not always
            continue

        if os.path.isfile(sfam_path+"_bsa.h5"):
            try:
                store = pd.HDFStore(unicode(sfam_path+"_bsa.h5"))
                keys = store.keys()
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                job.log("FAILED READING {}".format(sfam_id))
                os.remove(sfam_path+"_bsa.h5")
                keys = []

            if "/observed" in keys or store.get("/observed")["ppi_type"].isnull().sum() > 0:
                store.close()
                os.remove(sfam_path+"_bsa.h5")
                keys = []
            elif "/inferred" in keys or store.get("/inferred")["ppi_type"].isnull().sum() > 0:
                df = store.get("/observed")
                df.to_hdf(unicode(sfam_path+"_bsa.h52"), "/observed")
                store.close()
                os.remove(sfam_path+"_bsa.h5")
                shutil.move(sfam_path+"_bsa.h52", sfam_path+"_bsa.h5")
                keys = ["/observed"]
            else:
                store.close()
                keys = []
        else:
            keys = []

        if "/observed" not in keys:
            cjob = job.addChildJobFn(observed_bsa, dataset_name, sfam_id, memory=memory)
            cjob.addFollowOnJobFn(inferred_bsa, dataset_name, sfam_id, memory=memory)
        elif "/inferred" not in keys:
            cjob = job.addChildJobFn(inferred_bsa, dataset_name, sfam_id, memory=memory)

    job.addFollowOnJobFn(start_toil, dataset_name, iteration=1)

        # cjob = job.addChildJobFn(observed_bsa, dataset_name, sfam_id)
        # if not os.path.isfile(sfam_path+".inferred_interactome"):
        #     continue
        # cjob.addFollowOnJobFn(inferred_bsa, dataset_name, sfam_id)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    dataset_name = options.jobStore.split(":")[-1]

    job = Job.wrapJobFn(start_toil, dataset_name)
    with Toil(options) as toil:
        toil.start(job)
