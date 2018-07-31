import os, sys

import subprocess
from itertools import izip
import gzip
import tempfile
import re
import shlex
import json

import pandas as pd
import numpy as np
import dask.dataframe as dd
from joblib import Parallel, delayed

import freesasa

from Bio.PDB import PDBList

from util import data_path_prefix, get_interfaces_path, iter_cdd

NUM_WORKERS = 4
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
    tmp_pdb_path = tmp_pdb.name
    pdb_path = os.path.join(RAW_PDB_PATH, pdb[1:3].lower(), "pdb{}.ent.gz".format(pdb.lower()))
    if os.path.isfile(pdb_path):
        with tmp_pdb as tmp, gzip.open(pdb_path, 'rb') as pdbf:
            for line in pdbf:
                tmp.write(line)
    else:
        PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=os.getcwd(), file_format="pdb")
        tmp_pdb_path = "pdb{}.ent".format(pdb.lower())

        if not os.path.isfile(tmp_pdb_path) or os.path.getsize(tmp_pdb_path) == 0:
            raise IOError("Cannot load "+pdb)

    tmp_pdb.close()

    prefix = "freesasa_{}_{}_{}_d1d2.pdb".format(pdb, chain1, chain2)
    tmp_pdb2 = tempfile.NamedTemporaryFile(prefix=prefix, delete=False)

    with open(tmp_pdb2.name, "w") as output:
        splitchain1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain1), tmp_pdb_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        splitdomain1 = subprocess.Popen(
            [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+list(sdi1),
            stdin=splitchain1.stdout,
            stdout=output,
            stderr=subprocess.PIPE)
        splitdomain1.communicate()
        if chain2 is not None and sdi2 is not None:
            splitchain2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_selchain.py"), "-{}".format(chain2), tmp_pdb_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            splitdomain2 = subprocess.Popen(
                [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]+list(sdi2),
                stdin=splitchain2.stdout,
                stdout=output,
                stderr=subprocess.PIPE)
            splitdomain2.communicate()

    os.remove(tmp_pdb.name)

    return tmp_pdb2.name

def run_freesasa(command):
    FNULL = open(os.devnull, 'w')
    freesasa = subprocess.check_output(command, stderr=FNULL)
    freesasa = "{"+freesasa.split("\n",1)[1]
    return json.loads(freesasa)

def calculate_buried_surface_area(pdb_file, pdb, chain1, chain2, sdi_sel1=None, sdi_sel2=None, face1=None, face2=None):
    chains = [chain1+chain2, chain1, chain2]

    command = ["freesasa", "--format=json", "--chain-groups={}".format("+".join(chains))]

    if sdi_sel1 is not None and sdi_sel1 is not None:
        command.append("--select='domain1, chain {} and resi {}'".format(chain1, "+".join(sdi_sel1)))
        command.append("--select='domain2, chain {} and resi {}'".format(chain2, "+".join(sdi_sel2)))

    if face1 is not None:
        residues = face1.replace(",", "+")
        selection = "binding-site1, chain {} and resi {}".format(chain1, residues)
        command.append("--select='{}'".format(selection))

    if face2 is not None:
        residues = face2.replace(",", "+")
        selection = "binding-site2, chain {} and resi {}".format(chain2, residues)
        command.append("--select='{}'".format(selection))

    command.append(pdb_file)

    freesasa = run_freesasa(command)

    if sdi_sel1 is not None and sdi_sel1 is not None:
        try:
            c1_asa = freesasa["results"][2]["structure"][0]["selections"][0]["area"]
        except (IndexError, KeyError):
            #Use full chain area, it is assumed the interacting domains are split out
            c1_asa = freesasa["results"][2]["structure"][0]["area"]["total"]

        try:
            c2_asa = freesasa["results"][3]["structure"][0]["selections"][0]["area"]
        except (IndexError, KeyError):
            #Use full chain area, it is assumed the interacting domains are split out
            c2_asa = freesasa["results"][3]["structure"][0]["area"]["total"]

        complex_asa = sum([sel["area"] for sel in freesasa["results"][1]["structure"][0]["selections"] \
            if "domain" in sel["name"]])
        print complex_asa, freesasa["results"][1]["structure"][0]["area"]["total"]

        if complex_asa == 0.0:
            complex_asa = freesasa["results"][1]["structure"][0]["area"]["total"]
            print "Error, used total", complex_asa, freesasa["results"][1]["structure"][0]["selections"], command
    else:
        c1_asa = freesasa["results"][2]["structure"][0]["area"]["total"]
        c2_asa = freesasa["results"][3]["structure"][0]["area"]["total"]
        complex_asa = freesasa["results"][1]["structure"][0]["area"]["total"]

    bsa = c1_asa+c2_asa-complex_asa

    for ppi_type, (low_cut, high_cut) in cutoffs.iteritems():
        if low_cut <= bsa < high_cut:
            return bsa, ppi_type, c1_asa, c2_asa, complex_asa
    else:
        return bsa, "unknown", c1_asa, c2_asa, complex_asa

def calculate_surface_area_chain(pdb_file, pdb, chain, sdi=None, face=None):
    sdi_from, sdi_to = sdi

    command = ["freesasa", "--chain-groups={}".format(chain)]

    if sdi_sel is not None:
        command.append("--select=domain, resi {}".format(sdi))

    if face is not None:
        residues = face.replace(",", "+")
        selection = "binding-site, chain {} and resi {}".format(chain, residues)
        command.append("--select={}".format(selection))

    command.append(pdb_file)

    print " ".join(command)

    freesasa = run_freesasa(command)

    if sdi_sel is not None:
        try:
            return freesasa["results"][1]["structure"][0]["selections"][0]["area"]
        except (IndexError, KeyError):
            pass
    return freesasa["results"][1]["structure"][0]["area"]["total"]

def get_bsa(df):
    r = df.iloc[0]

    if any(r[["mol_chain", "int_chain"]].isna()):
        #Cannot process chain, so bsa is unknown
        return pd.Series({"bsa":np.nan, "ppi_type":"unknown", "c1_asa":np.nan, "c2_asa":np.nan, "complex_asa":np.nan})

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
        return pd.Series({"bsa":np.nan, "ppi_type":"unknown", "c1_asa":np.nan, "c2_asa":np.nan, "complex_asa":np.nan})

    try:
        pdb_file = get_pdb(r["mol_pdb"], r["mol_chain"], r["int_chain"], sdi1, sdi2)
    except IOError:
        return pd.Series({"bsa":np.nan, "ppi_type":"unknown", "c1_asa":np.nan, "c2_asa":np.nan, "complex_asa":np.nan})

    bsa, ppi_type, c1_asa, c2_asa, complex_asa = calculate_buried_surface_area(
        pdb_file, r["mol_pdb"], r["mol_chain"], r["int_chain"],
        face1=r["mol_res"], face2=r["int_res"])

    try:
        os.remove(pdb_file)
    except OSError:
        pass

    return pd.Series({"bsa":bsa, "ppi_type":ppi_type, "c1_asa":c1_asa, "c2_asa":c2_asa, "complex_asa":complex_asa})

def observed_bsa(job, dataset_name, cdd, cores=NUM_WORKERS):
    cdd = cdd.replace("/", "")
    prefix = os.path.join(get_interfaces_path(dataset_name), cdd, cdd)
    cdd_interactome_path = prefix+".observed_interactome"

    cdd_interactome = pd.read_hdf(unicode(cdd_interactome_path), "table")
    cdd_interactome = cdd_interactome[cdd_interactome["mol_chain"] != cdd_interactome["int_chain"]]

    #Remove redundant interfaces
    cdd_interactome = cdd_interactome.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
        as_index=False).apply(lambda df: df.iloc[0]).reset_index(drop=True).copy()

    bsa = Parallel(n_jobs=NUM_WORKERS)(delayed(get_bsa)(group) for _, group in \
        cdd_interactome.groupby("mol_sdi_id", as_index=False))
    bsa = pd.DataFrame(bsa)

    cdd_interactome = pd.concat([cdd_interactome, bsa], axis=1)

    # df = dd.from_pandas(cdd_interactome, npartitions=NUM_WORKERS)
    #
    # meta = pd.DataFrame({"bsa":[1.], "ppi_type":["0"], "c1_asa":[1.], "c2_asa":[1.], "complex_asa":[1.]})
    # bsa = df.map_partitions(lambda _df: _df.apply(get_bsa, axis=1), meta=meta)\
    #     .compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)

    # cdd_interactome[bsa.columns] = bsa

    cdd_interactome.to_hdf(unicode(prefix+"_bsa.h5"), "observed", complevel=9, complib="bzip2")

def get_asa(df):
    r = df.iloc[0]
    sdi1 = ["{}-{}".format(int(row.mol_sdi_from), int(row.mol_sdi_to)) for row in df.itertuples()]
    try:
        pdb_file = get_pdb(r["mol_pdb"], r["mol_chain"], None, sdi1, None)
    except IOError:
        return pd.Series({"c2_asa_pred":np.nan, "pred_ratio":np.nan})

    asa = calculate_surface_area_chain(
        pdb_file, r["mol_pdb"], r["mol_chain"], face=row["mol_resi"])

    try:
        os.remove(pdb_file)
    except OSError:
        pass

    predicted_bsa = asa+r["c2_asa"]-r["complex_asa"]
    ratio = predicted_bsa/r["bsa"]

    for ppi_type, (low_cut, high_cut) in cutoffs.iteritems():
        if low_cut <= predicted_bsa < high_cut:
            break
    else:
        ppi_type = "unknown"

    return pd.Series({"c2_asa_pred":predicted_bsa, "pred_ratio":ratio, "ppi_type_pred":ppi_type})

def inferred_bsa(job, dataset_name, cdd, cores=NUM_WORKERS):
    cdd = cdd.replace("/", "")
    cdd_bsa_path = os.path.join(get_interfaces_path(dataset_name), cdd, cdd)
    cdd_obs_bsa = pd.read_hdf(unicode(cdd_bsa_path+"_bsa.h5"), "observed")
    cdd_obs_bsa = cdd_obs_bsa[["obs_int_id", "bsa", "c1_asa", "c2_asa", "complex_asa", "ppi_type"]]

    inf_interactome = pd.read_hdf(unicode(cdd_bsa_path+".inferred_interactome"), "table")

    inf_interactome = pd.merge(inf_interactome, cdd_obs_bsa, how="left", left_on="nbr_obs_int_id", right_on="obs_int_id")
    del inf_interactome["obs_int_id"]

    #Remove redundant interfaces
    inf_interactome = cdd_interactome.groupby(["obs_int_id", "mol_sdi_from", "mol_sdi_to"],
        as_index=False).apply(lambda df: df.iloc[0]).reset_index(drop=True).copy()

    bsa = Parallel(n_jobs=NUM_WORKERS)(delayed(get_asa)(group) for _, group in \
        inf_interactome.groupby("mol_sdi_id", as_index=False))
    bsa = pd.DataFrame(bsa)

    cdd_interactome = pd.concat([cdd_interactome, bsa], axis=1)

    df = dd.from_pandas(inf_bsa, npartitions=NUM_WORKERS)

    meta = pd.DataFrame({"c2_asa_pred":[1.], "pred_ratio":[1.], "ppi_type_pred":["1"]})
    bsa = df.map_partitions(lambda _df: _df.apply(get_asa, axis=1), meta=meta)\
        .compute(scheduler="multiprocessing", num_workers=NUM_WORKERS)

    inf_bsa[bsa.columns] = bsa

    inf_bsa.to_hdf(unicode(cdd_bsa_path+"_bsa.h5"), "inferred", complevel=9, complib="bzip2")

def start_toil(job, dataset_name, name="bsa"):
    for cdd in iter_cdd():
        cjob = job.addChildJobFn(observed_bsa, dataset_name, cdd)
        #cjob.addFollowOnJobFn(inferred_bsa, dataset_name, cdd)

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

    #x/0.1 = 28
#
#     if len(sys.argv) == 2:
#         submit_ibis_cdd(sys.argv[1])
#     elif len(sys.argv) == 4 and sys.argv[2] == "obs":
#         observed_bsa(sys.argv[1], shlex.split(sys.argv[3])[0])
#     elif len(sys.argv) == 4 and sys.argv[2] == "inf":
#         inferred_bsa(sys.argv[1], shlex.split(sys.argv[3])[0])
#     else:
#         print len(sys.argv), sys.argv
#
# def submit_ibis_cdd(dataset_name, job_name="calc_bsa", dependency=None):
#     obsjob = SlurmJob(job_name+"_obs", cpus=14, walltime="8:00:00")
#     infjob = SlurmJob(job_name+"_inf", cpus=14, walltime="8:00:00")
#     for cdd in iter_cdd():
#         obsjob += "{} {} {} obs \"{}\"\n".format(sys.executable, __file__, dataset_name, cdd)
#         infjob += "{} {} {} inf \"{}\"\n".format(sys.executable, __file__, dataset_name, cdd)
#
#
#     obs_jid = obsjob.run(dependency=dependency)
#     #print obs_jid
#
#     inf_jid = obsjob.run(dependency="afterok:"+obs_jid)
#     #print inf_jid
#     return inf_jid
