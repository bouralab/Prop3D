import os, sys
import glob
import traceback
from itertools import groupby
import numpy as np
import pandas as pd
from calculate_features import SwarmJob
#from molmimic.biopdbtools import Structure

def split_ibis(f):
    f.next()
    for line in f:
        yield line.rstrip().split(":")

def get_ibis(pdb_ibis_file, multimers=False, max_dist=256.0, use_structure_filters=False, calculate_features=False, save_binding_site_sequences=True):
    with open(pdb_ibis_file) as pdb_ibis:
        for pdb_chain, entries in groupby(split_ibis(pdb_ibis), key=lambda l: l[0]):
            pdb, chain = pdb_chain[:-1], pdb_chain[-1]

            info = []
            has_interactions = False
            for entry in entries:
                if entry[1] not in ["PPI"]: continue #, "LIG"

                # if not multimers and entry[-6] == entry[-2]:
                #     #Query domain should not be target domain
                #     continue

                residues = entry[3].lstrip().replace(" ", ",")

                #is_ppi = "1" if entry[1] == "PPI" else "0"
                is_multimer = "1" if entry[-6] == entry[-2] else "0"
                cdd = entry[-2]

                with open("/data/draizene/molmimic/molmimic/data/ibis_full_by_cdd/{}.tsv".format(cdd.replace("/", "_")), "a+") as f:
                    print >> f, "{}\t{}\t{}\t{}\t{}".format(pdb, chain, residues, is_multimer, entry[-2])

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


def create_ibis(ibis_dir, csv=True, pdb=False):
    files = glob.glob("{}/*/*.txt".format(ibis_dir))
    pdbs = []
    for i, ibis in enumerate(files):
        pdb = os.path.splitext(os.path.basename(ibis))[0]
        if csv:
            out_path = "/data/draizene/molmimic/molmimic/data/ibisdown_data/{}".format(pdb[1:3].lower())
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            get_ibis(ibis)
        if pdb:
            pdbs.append(pdb)


def submit_ibis(ibis_data):
    job = SwarmJob("make_ibis", individual=True)
    job.add_individual_parameters()
    job += "python {} create {}\n".format(__file__, ibis_data)
    print job.submit_individual()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        submit_ibis("/data/draizene/molmimic/molmimic/data/ibisdown/")
    elif len(sys.argv) == 2:
        submit_ibis(sys.argv[1])
    elif len(sys.argv) == 3 and "create" in sys.argv:
        create_ibis(sys.argv[2])
