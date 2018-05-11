import os, sys
import glob
import traceback
from itertools import groupby
import numpy as np
import pandas as pd
from calculate_features import SwarmJob
from Bio import PDB
from molmimic.biopdbtools import Structure, InvalidPDB
from molmimic.util import initialize_data

def split_ibis(f):
    f.next()
    for line in f:
        yield line.rstrip().split(":")

def get_ibis(dataset_name, pdb_ibis_file, use_cdd_domain=None, multimers=False, check_binding_sites=False):
    out_dir = os.path.join(os.path.dirname(__file__), "..", "data", "interfaces", dataset_name)
    with open(pdb_ibis_file) as pdb_ibis:
        for pdb_chain, entries in groupby(split_ibis(pdb_ibis), key=lambda l: l[0]):
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
                if entry[1] not in ["PPI"]: continue #, "LIG"

                is_multimer = "1" if entry[-6] == entry[-2] else "0"
                cdd = entry[-2]
                residues = entry[3].lstrip().replace(" ", ",")
                residue_str = entry[4]

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
                        print "{} {} does not match {} =/= {}".format(entry[0], residues, pdb_seq, residue_str)
                        continue

                with open(os.path.join(out_dir, "{}.raw.tsv".format(cdd.replace("/", ""))), "a+") as f:
                    print >> f, "{}\t{}\t{}\t{}\t{}".format(pdb, chain, residues, is_multimer, cdd)

                #info.append("{}\t{}\t{}\t{}\t{}".format(pdb, chain, residues, is_multimer, entry[-2]))

                #has_interactions = True

            # if not has_interactions:
            #     continue
            #
            # if use_structure_filters:
            #     try:
            #         structure = Structure(pdb, chain, force_feature_calculation=True, unclustered=True)
            #     except (KeyboardInterrupt, SystemExit) as e:
            #         raise
            #     except:
            #         trace = traceback.format_exc()
            #         print "Error:", trace
            #         continue
            #
            #     coords = structure.get_coords()
            #     min_coord = np.min(coords, axis=0)
            #     max_coord = np.max(coords, axis=0)
            #     dist = int(np.ceil(np.linalg.norm(max_coord-min_coord)))
            #
            #     if dist > max_dist:
            #         continue
            #
            #     if calculate_features:
            #         for atom in structure.structure.get_atoms():
            #             structure.get_features_for_atom(atom, preload=False)
            #
            #         structure.precalc_features.flush()
            #
            # if save_binding_site_sequences:
            #     with open("/data/draizene/molmimic/molmimic/data/ibisdown_data/{}/{}_{}.tsv".format(pdb[1:3].lower(), pdb, chain), "w") as f:
            #         for line in info:
            #             print >> f, line

def post_process(force=True):
    import re
    def get_digits(s):
        try:
            return int(re.findall('\d+', s)[0])
        except (IndexError, ValueError):
            print s
            raise

    full_data_dir = "/data/draizene/molmimic/molmimic/data"
    full_data = full_data_dir+"/ibis_full.tsv"
    if not os.path.exists(full_data) or force:
        with open(full_data, "w") as out:
            print >> out, "\t".join(["pdb", "chain", "residues", "ppi"])

            for f in glob.glob(full_data_dir+"/ibisdown_data/*"):
                with open(f) as info:
                    for line in info:
                        print >> out, line.rstrip()

    ibis_full = pd.read_table(full_data)
    ibis_full_ppi = ibis_full[ibis_full["ppi"] == 1]
    #ibis_full_pli = ibis_full[ibis_full["ppi"] == 0]

    for name, df in [("ppi", ibis_full_ppi)]: #, ("pli", ibis_full_pli)]:
        groups = df.groupby(["pdb", "chain"])

        with open(full_data_dir+"/ibis_{}.tsv".format(name), "w") as f:
            print >> f, "pdb\tchain\t{}_residues".format(name)
            for (pdb, chain), group in groups:
                print >> f, "{}\t{}\t{}".format(pdb, chain, ",".join(map(str, set(sorted(map(get_digits, ",".join(group["residues"]).split(",")))))))

    #ppi = pd.read_table(full_data_dir+"/ibis_ppi.tsv")
    #pli = pd.read_table(full_data_dir+"/ibis_pli.tsv")

    #ppi_pli = pd.merge(ppi, pli, on=['pdb', 'chain'])
    #ppi_pli.to_csv(full_data_dir+"/ibis_ppi_pli.tsv", sep="\t", index=False)

def post2():
    pass


def create_ibis(dataset_name, ibis_dir, cdd=None):
    initialize_data(dataset_name)
    files = glob.glob("{}/*/*.txt".format(ibis_dir))
    pdbs = []
    for i, ibis in enumerate(files):
        get_ibis(dataset_name, ibis, use_cdd_domain=cdd)


def submit_ibis(dataset_name, ibis_data, job_name="build_ibis", dependency=None):
    out_path = os.path.join(os.path.dirname(__file__), "..", "data", "interfaces", dataset_name)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    job = SwarmJob(job_name, individual=True)
    job += "python {} create {} {}\n".format(__file__, dataset_name, ibis_data)
    job_id = job.submit_individual(dependency=dependency)
    return job_id

if __name__ == "__main__":
    if len(sys.argv) == 2:
        submit_ibis(sys.argv[1], "/data/draizene/molmimic/molmimic/data/ibisdown/")
    elif len(sys.argv) == 3:
        submit_ibis(sys.argv[1], sys.argv[2])
    elif len(sys.argv) in [4, 5] and "create" in sys.argv:
        create_ibis(*sys.argv[2:])
    elif len(sys.argv) == 4 and "get" in sys.argv:
        get_ibis(sys.argv[-2], sys.argv[-1])
