import os, sys
import shutil
import glob
import traceback
from itertools import groupby
import numpy as np
import pandas as pd
from calculate_features import SwarmJob
from Bio import PDB
from molmimic.biopdbtools import Structure, InvalidPDB
from molmimic.util import initialize_data, get_interfaces_path

def split_ibis(f):
    for line in f:
        yield line.rstrip().split(":")

def parse_ibis(pdb_ibis):
    header = next(pdb_ibis).split(":")
    for entries in split_ibis(pdb_ibis):
        yield {header[i]:value for i, value in enumerate(entries)}

def get_ibis_path(pdb):
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "interfaces", "ibisdown", pdb[1:3].upper(), pdb.upper()+".txt")

def parse_ibis_from_pdb(pdb):
    path = get_ibis_path(pdb)
    with open(path) as f:
        for entry in parse_ibis(f):
            yield entry

def get_ibis(dataset_name, pdb_ibis_file, use_cdd_domain=None, multimers=False, check_binding_sites=False):
    out_dir = get_interfaces_path(dataset_name)
    seen_cdd = set()
    with open(pdb_ibis_file) as pdb_ibis:
        for pdb_chain, entries in groupby(parse_ibis(pdb_ibis), key=lambda l: l["Query"]):
            pdb, chain = pdb_chain[:-1], pdb_chain[-1]

            if check_binding_sites:
                try:
                    structure = Structure.from_pdb(pdb, chain)
                except (KeyboardInterrupt, SystemExit) as e:
                    raise
                except InvalidPDB:
                    continue
                except:
                    trace = traceback.format_exc()
                    print "Error:", trace
                    continue

            for entry in entries:
                if entry["Interaction_type"] not in ["PPI"]: continue #, "LIG"

                cdd = entry["Query_Domain"]
                partner = entry["Interaction_Partner"]
                residues = entry["PDB_Residue_No"].lstrip().replace(" ", ",")
                residue_str = entry["Binding_Site_Residues"]
                observed = entry["Is_Observed"] == "1"
                pdb_evidence = entry["PDB_Evidence"]
                is_multimer = "1" if cdd == partner else "0"


                if multimers and not is_multimer:
                    #Query domain should not be target domain
                    continue

                if use_cdd_domain is not None and cdd != use_cdd_domain:
                    continue

                if check_binding_sites:
                    #Check if positions match structure
                    pdb_seq = ""
                    for i, r in enumerate(residues.split(",")):
                        if r == "X": continue

                        res = structure.get_residue_from_resseq(r)

                        if res is None:
                            residue_str = residue_str[:i]+residue_str[i+1:]
                        else:
                            pdb_seq += PDB.Polypeptide.three_to_one(res.get_resname())

                    if pdb_seq != residue_str:
                        print "{} {} does not match {} =/= {}".format(entry["Query"], residues, pdb_seq, residue_str)
                        continue

                seen_cdd.add(cdd)

                with open(os.path.join(out_dir, cdd.replace("/", ""), "{}.tsv".format(pdb)), "a+") as f:
                    print >> f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(pdb, chain, residues, residue_str, is_multimer, cdd, partner, pdb_evidence, int(observed))

    # CDD = pd.read_csv("MMDB/StructDomSfam.csv", usecols=["label"]).drop_duplicates().dropna()
    # CDD = set(CDD["label"].tolist())
    # cdds_with_interfaces = CDD-seen_cdd

    # for cdd_no_ibis in cdds_with_interfaces:
    #     with open(os.path.join(out_dir, "{}.raw".format(cdd_no_ibis.replace("/", ""))), "w"):
    #         pass

def merge_ibis(dataset_name, cdd):
    out_dir = get_interfaces_path(dataset_name)
    with open(os.path.join(out_dir, "{}.raw".format(cdd.replace("/", ""))), "w") as raw_cdd:
        print >> raw_cdd, "pdb\tchain\tresi\tresn\tis_multimer\tcdd\tpartner\tpdb_evidence\tis_observed"
        for pdb in glob.glob(os.path.join(out_dir, cdd.replace("/", ""), "*.tsv")):
            with open(pdb) as f:
                for line in f:
                    raw_cdd.write(line)
    shutil.rmtree(os.path.join(out_dir, cdd.replace("/", "")))

def create_ibis(dataset_name, ibis_dir, cdd=None):
    initialize_data(dataset_name)
    if os.path.basename(ibis_dir) == "ibisdown":
        ibis_dir = os.path.join(ibis_dir, "*")
    files = glob.glob(os.path.join(ibis_dir, "*.txt"))
    for i, ibis in enumerate(files):
        get_ibis(dataset_name, ibis, use_cdd_domain=cdd)

def submit_ibis(dataset_name, ibis_data, job_name="build_ibis", dependency=None):
    out_path = os.path.join(os.path.dirname(__file__), "..", "data", "interfaces", dataset_name)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    CDD = pd.read_csv("MMDB/StructDomSfam.csv", usecols=["label"]).drop_duplicates().dropna()
    CDD = set(CDD["label"].tolist())
    for cdd in CDD:
        if not os.path.exists(os.path.join(out_path, cdd.replace("/", ""))):
            os.makedirs(os.path.join(out_path, cdd.replace("/", "")))

    job = SwarmJob(job_name, walltime="1-00:00:00", mem="5")

    pdb_groups = next(os.walk(ibis_data))[1]
    for pdb_group in pdb_groups:
        job += "python {} create {} {}\n".format(__file__, dataset_name, os.path.join(ibis_data, pdb_group))
    jid1 = job.run(dependency=dependency)

    merge_job = SwarmJob(job_name, walltime="6:00:00", mem="5")
    for cdd in CDD:
        merge_job += "python {} merge {} \"{}\"\n".format(__file__, dataset_name, cdd.replace("/", ""))
    jid2 = merge_job.run(dependency="afterany:"+jid1)

    return jid2

if __name__ == "__main__":
    if len(sys.argv) == 2:
        submit_ibis(sys.argv[1], "/data/draizene/molmimic/molmimic/data/ibisdown/")
    elif len(sys.argv) == 3:
        submit_ibis(sys.argv[1], sys.argv[2])
    elif len(sys.argv) in [4, 5] and "create" in sys.argv:
        create_ibis(*sys.argv[2:])
    elif len(sys.argv) == 4 and "merge" in sys.argv:
        merge_ibis(*sys.argv[2:])
    elif len(sys.argv) == 4 and "get" in sys.argv:
        get_ibis(sys.argv[-2], sys.argv[-1])
