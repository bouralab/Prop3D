import os, sys
sys.path.append("/data/draizene/pdb-tools")
sys.path.append("/data/draizene/molmimic")

import subprocess
from itertools import izip
import gzip
import tempfile
import re

import pandas as pd
import numpy as np

from Bio.PDB import PDBList

from molmimic.calculate_features import SwarmJob
from molmimic.calculate_ibis_dataset import parse_ibis_from_pdb
from molmimic.util import get_interfaces_path, iter_cdd

cutoffs = {
    "weak transient": (0, 1500),
    "transient": (1500, 2500),
    "permanent": (2500, float("inf"))
}

def get_pdb(pdb, chain1, chain2):
    tmp_pdb_path = "freesasa_{}_{}_{}.pdb".format(pdb, chain1, chain2)
    pdb_path = '/pdb/pdb/{}/pdb{}.ent.gz'.format(pdb[1:3].lower(), pdb.lower())
    if os.path.isfile(pdb_path):
        with open(tmp_pdb_path, 'w') as tmp, gzip.open(pdb_path, 'rb') as pdbf:
            for line in pdbf:
                tmp.write(line)
    else:
        PDBList().retrieve_pdb_file(pdb, obsolete=True, pdir=os.getcwd(), file_format="pdb")
        tmp_pdb_path = "pdb{}.ent".format(pdb.lower())
        if os.path.getsize(tmp_pdb_path) == 0:
            raise IOError("Cannot load "+pdb)

    return tmp_pdb_path

def calculate_buried_surface_area(pdb_file, pdb, chain1, chain2, sdi1, sdi2, face1=None, face2=None):
    chains = [chain1+chain2, chain1, chain2]
    sdi1_from, sdi1_to = sdi1
    sdi2_from, sdi2_to = sdi2

    command = ["freesasa", "--chain-groups={}".format("+".join(chains))]
    command.append("--select='domain1, chain {} and resi {}-{}'".format(chain1, sdi1_from, sdi1_to))
    command.append("--select='domain2, chain {} and resi {}-{}'".format(chain2, sdi2_from, sdi2_to))

    if face1 is not None:
        residues = face1.replace(",", "+")
        selection = "binding-site1, chain {} and resi {}".format(chain1, residues)
        command.append("--select='{}'".format(selection))

    if face2 is not None:
        residues = face2.replace(",", "+")
        selection = "binding-site2, chain {} and resi {}".format(chain2, residues)
        command.append("--select='{}'".format(selection))

    command.append(pdb_file)

    FNULL = open(os.devnull, 'w')
    freesasa = subprocess.check_output(command, stderr=FNULL)

    freesasa = iter(freesasa.splitlines())
    chain_results = {}
    for line in freesasa:
        try:
            chain_result = process_chain_section(freesasa)
        except StopIteration:
            break
        chain_results[chain_result["chain_names"]] = chain_result

    try:
        c1_asa = chain_results[frozenset([chain1])]["selections"]["domain1"]
    except KeyError:
        c1_asa = chain_results[frozenset([chain1])]["total"]

    try:
        c2_asa = chain_results[frozenset([chain2])]["selections"]["domain2"]
    except KeyError:
        c2_asa = chain_results[frozenset([chain2])]["total"]

    try:
        complex_asa = chain_results[frozenset(chain1+chain2)]["selections"]
        complex_asa = complex_asa["selections"]["domain1"]+complex_asa["selections"]["domain2"]
    except KeyError:
        complex_asa = chain_results[frozenset(chain1+chain2)]["total"]

    bsa = c1_asa+c2_asa-complex_asa

    for ppi_type, (low_cut, high_cut) in cutoffs.iteritems():
        if low_cut <= bsa < high_cut:
            return bsa, ppi_type, c1_asa, c2_asa, complex_asa
    else:
        return bsa, "unknown", c1_asa, c2_asa, complex_asa

def calculate_surface_area_chain(pdb_file, chain, sdi, face=None):
    sdi_from, sdi_to = sdi

    command = ["freesasa", "--separate-chains"]
    command.append("--select=domain, resi {}-{}".format(sdi_from, sdi_to))

    if face is not None:
        residues = face.replace(",", "+")
        selection = "binding-site, chain {} and resi {}".format(chain, residues)
        command.append("--select={}".format(selection))

    command.append(pdb_file)

    command.append(pdb_file)

    FNULL = open(os.devnull, 'w')
    freesasa = subprocess.check_output(command, stderr=FNULL)

    freesasa = iter(freesasa.splitlines())
    chain_results = {}
    for line in freesasa:
        try:
            chain_result = process_chain_section(freesasa)
        except StopIteration:
            print chain_results

        chain_results[chain_result["chain_names"]] = chain_result

    try:
        return chain_results[frozenset([chain])]["selections"]["domain"]
    except KeyError:
        return chain_results[frozenset([chain])]["total"]


def process_chain_section(freesasa_chain):
    line = next(freesasa_chain)
    while not line.startswith("INPUT"):
        line = next(freesasa_chain)

    chain = {}
    chain["source"] = next(freesasa_chain).split(":")[1].strip()
    chain["chain_names"] = frozenset(list(next(freesasa_chain).split(":")[1].strip()))
    chain["model"] = int(next(freesasa_chain).split(":")[1].strip())
    chain["atoms"] = int(next(freesasa_chain).split(":")[1].strip())

    line = next(freesasa_chain)
    while not line.startswith("RESULTS"):
        line = next(freesasa_chain)

    chain["total"] = float(next(freesasa_chain).split(":")[1].strip())
    chain["Apolar"] = float(next(freesasa_chain).split(":")[1].strip())
    chain["polar"] = float(next(freesasa_chain).split(":")[1].strip())
    chain["chains"] = {}
    for c in chain["chain_names"]:
        line = next(freesasa_chain)
        chain["chains"][c] = float(line.split(":")[1].strip())

    process_selections = False
    while True:
        try:
            line = next(freesasa_chain)
        except StopIteration:
            break

        if line.startswith("SELECTIONS"):
            process_selections = True
            chain["selections"] = {}
            continue
        elif line.startswith("#"):
            break

        if process_selections:
            try:
                name, value = line.replace(" ", "").split(":")
            except ValueError:
                continue
            chain["selections"][name] = float(value)


    return chain

def get_bsa(row):
    pdb_file = get_pdb(row["mol_pdb"], row["mol_chain"], row["int_chain"])
    bsa, ppi_type, c1_asa, c2_asa, complex_asa = calculate_buried_surface_area(
            pdb_file, row["mol_pdb"], row["mol_chain"], row["int_chain"],
            (row["mol_sdi_from"], row["mol_sdi_to"]),
            (row["int_sdi_from"], row["int_sdi_to"]),
            row["mol_res"], row["int_res"])
    os.remove(pdb_file)
    return pd.Series({"bsa":bsa, "ppi_type":ppi_type, "c1_asa":c1_asa, "c2_asa":c2_asa, "complex_asa":complex_asa})

def observed_bsa(dataset_name, cdd):
    import dask.dataframe as dd
    from dask.multiprocessing import get

    prefix = os.path.join(get_interfaces_path(dataset_name), cdd.replace("/", ""))
    cdd_interactome_path = prefix+".observed_interactome"

    cdd_interactome = pd.read_hdf(cdd_interactome_path, "table") #read_table(cdd_interactome_path, header=None, names=names)
    cdd_interactome = cdd_interactome[cdd_interactome["mol_chain"] != cdd_interactome["int_chain"]]

    df = dd.from_pandas(cdd_interactome, npartitions=14)

    meta = pd.DataFrame({"bsa":[0], "ppi_type":[0], "c1_asa":[0], "c2_asa":[0], "complex_asa":[0]})
    bsa = df.map_partitions(lambda _df: _df.apply(get_bsa, axis=1), meta=meta).compute(get=get)

    cdd_interactome[bsa.columns] = bsa

    cdd_interactome.to_hdf(prefix+".observed_bsa", "table", complevel=9, complib="bzip2")

def inferred_bsa(dataset_name, cdd, chunksize=100):
    cdd_bsa_path = os.path.join(get_interfaces_path(dataset_name), "{}")
    cdd_obs_bsa = pd.read_table("{}.obs_bsa".format(cdd_bsa_path, bsa_type))

    st_domain_file = os.path.join(mmdb_path, "StDomain.csv")
    st_domain = pd.read_csv(st_domain_file, usecols=["sdi", "pdbId", "chnLett"], dtype="str")

    st_domain_intvl_file = os.path.join(mmdb_path, "StDomainIntvl.csv")
    st_domain_intvl = pd.read_csv(st_domain_intvl_file, usecols=["sdi", "pdbId", "chnLett"], dtype="str")
    st_domain_intvl.columns = ["sdi", "sdi_from", "sdi_to"]
    st_domain = pd.merge(st_domain, st_domain_intvl, how="left", on="sdi")

    inf_interactome = pd.read_csv(cdd_interactome_path+".inferred_interactome", chunksize=chunksize)

    for inf_int in inf_interactome:
        inf_int = pd.merge(inf_int, st_domain, how="left", left_on="mol_sdi_id", right_on="sdi")
        inf_int = pd.merge(inf_int, cdd_obs_interactome, how="left", left_on="nbr_obs_int_id", right_on="obs_int_id", suffixes=['_inf', '_obs'])
        
        for _, row in inf_int.iterrows():
            pdb_file = get_pdb(row["pdbId_inf"], row["chnLett_inf"], None)
            asa = calculate_surface_area_chain(
                pdb_file, row["pdbId_inf"], row["chnLett_inf"], 
                (row["sdi_from_inf"], row["sdi_to_inf"]),
                row["seqloc_inf"])

            predicted_bsa = asa+inf_int["c2_asa"]-inf_interactome["complex_asa"]
            ratio = predicted_bsa/inf_int["bsa"]

            if np.isclose(predicted_bsa, 1.0, atol=0.01):
                pass

def submit_ibis_cdd(dataset_name, obs, job_name="calc_bsa", dependency=None):
    job = SwarmJob(job_name, cpus=14, walltime="8:00:00")
    for cdd in iter_cdd():
        job += "/data/draizene/3dcnn-torch-py2 python {} {} {} \"{}\"\n".format(__file__, dataset_name, "obs" if obs else "inf", cdd)

    jid = job.run(dependency=dependency)
    print jid

if __name__ == "__main__":
    if len(sys.argv) == 2:
        submit_ibis_cdd(sys.argv[1])
    elif len(sys.argv) == 3 and sys.argv[1] == "obs":
        submit_ibis_cdd(sys.argv[2], obs=True)
    elif len(sys.argv) == 4 and sys.argv[1] == "obs":
        observed_bsa(sys.argv[2], sys.argv[3])
    else:
        print len(sys.argv), sys.argv
