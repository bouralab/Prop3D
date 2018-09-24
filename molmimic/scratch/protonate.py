import sys
sys.path.append("/data/draizene/molmimic")

import os
import argparse
import time
from datetime import datetime
from itertools import groupby
import glob

def check_structures():
    pdb_path = "/data/draizene/molmimic/pdb2pqr_structures/PKc_like/"

    protonated_minimized = set([tuple(os.path.basename(f[:-8]).split("_")) for f in \
        glob.glob(os.path.join(pdb_path, "*", "*_*.min.pdb"))])
    print "Num protonated and minimized", len(protonated_minimized)

    protonated = set([tuple(os.path.basename(f[:-8]).split("_")) for f in \
        glob.glob(os.path.join(pdb_path, "*", "*_*.pqr.pdb"))])
    print "Num protonated", len(protonated)

    raw_chains = set([tuple(os.path.basename(f[:-4]).split("_")) for f in \
        glob.glob(os.path.join(pdb_path, "*", "*_*.pdb")) \
        if not f.endswith(".min.pdb") and not f.endswith(".pqr.pdb")])
    print "Num raw chains", len(raw_chains)

    combined = set([os.path.basename(f[:-4]) for f in \
        glob.glob(os.path.join(pdb_path, "*", "*.pdb")) \
        if "_" not in os.path.basename(f) and not f.endswith(".min.pdb") and
        not f.endswith(".pqr.pdb") and not f.endswith(".pqr_0001.pdb") and not f.endswith(".raw.pdb")])
    print "Num combined chains", len(combined)

    skip = (protonated-protonated_minimized).union(raw_chains-protonated-protonated_minimized)
    with open("/data/draizene/molmimic/pdb2pqr_structures/skip_pdbs.txt", "w") as f:
        for pdb in skip:
            print >> f, pdb

    return protonated_minimized, protonated, raw_chains, combined, skip

def to_int(s):
    return int("".join([d for d in s if d.isdigit()]))

def filter_training_data():
    protonated_minimized, = check_structures()
    all_cdd = pd.read_table("/data/draizene/molmimic/molmimic/data/ibis_full_multimers/Ig_176677", sep=",")
    group = all_cdd.groupby(["pdb", "chain"])
    residues = group["residues"].apply(lambda x:",".join(x))
    r_frame = residues.to_frame()
    best = []
    for i, r in r_frame.iterrows():
        if (r.name[0].lower(), r.name[1]) in protonated_minimized:
            best.append(pd.Series({"pdb":r.name[0], "chain":r.name[1], "residues":r.residues}))
    best = pd.DataFrame(best)
    best["residues"] = best["residues"].apply(lambda x: ",".join(map(str, sorted(map(to_int, set(x.split(",")))))))

def protonate_minimize_old():
    #unzip, tidy, and remove hetatms
    job = "gunzip -c /pdb/pdb/{0}/pdb{1}.ent.gz | python /data/draizene/pdb-tools/pdb_tidy.py | python /data/draizene/pdb-tools/pdb_striphet.py > {2}/{1}.pdb ; ".format(pdb_divided.lower(), pdb.lower(), pdb_path)

    #Split PDBs
    job += "python /data/draizene/pdb-tools/pdb_splitchain.py {}/{}.pdb ; ".format(pdb_path, pdb.lower())

    #Loop over all chains
    job += "for f in `ls -1 {}/{}_*.pdb`; do echo $f; ".format(pdb_path, pdb.lower())

    job += "base=${f%.*}; "

    #Remove altLocs
    job += "python /data/draizene/pdb-tools/pdb_delocc.py $f > $f.delocc; "
    job += "rm $f; mv $f.delocc $f; "

    #Add hydogens for chain $f
    job += "/data/draizene/3dcnn-torch-py2 shell pdb2pqr --ff=parse --ph-calc-method=propka --chain --drop-water $f $base.pqr.pdb; "

    job += "minimize -s $base.pqr.pdb "
    job += "-run:min_type lbfgs_armijo_nonmonotone -run:min_tolerance 0.001 -ignore_zero_occupancy false "
    job += "-out:file:scorefile $base.sc "
    job += "-out:path:pdb {0} -out:path:score {0} ; ".format(pdb_path)


    job += "python /data/draizene/pdb-tools/pdb_stripheader.py $base.pqr_0001.pdb > $base.min.pdb; "
    job += "mv $base.pqr.pdb.propka $base.propka; "
    job += "done;"
    return job

def load_ibis(ibis_data, minimize=True, cellular_organisms=False):
    from molmimic.calculate_features import SwarmJob
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data, cellular_organisms=cellular_organisms)
    data = dataset.data
    job = SwarmJob("pdb2pqr", cpus=1)

    key = lambda p: p[1:3]
    pdbs = sorted(set(data["pdb"].values.tolist()), key=key)

    i = 0
    for pdb_divided, pdb_group in groupby(pdbs, key=key):
        pdb_path = "/data/draizene/molmimic/pdb2pqr_structures/H3/{}".format(pdb_divided.lower())
        feature_path = "/data/draizene/molmimic/features_Ig/atom/{}".format(pdb_divided.lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
        if not os.path.exists(feature_path):
            os.makedirs(feature_path)
        for pdb in pdb_group:
            i+=1
            #if len(glob.glob(os.path.join(pdb_path, "{}_*.min.pdb".format(pdb.lower())))) > 0:
            #    continue
            job += "python /data/draizene/molmimic/molmimic/prepare_protein.py /pdb/pdb/{}/pdb{}.ent.gz \n".format(pdb_divided.lower(), pdb.lower())
            print "Running protein_prepare {}/{}: {}".format(i, len(pdbs), pdb)

    job.run()

if __name__ == "__main__":
    cellular_organisms = len(sys.argv) >= 3
    print cellular_organisms
    load_ibis(sys.argv[-1], cellular_organisms=cellular_organisms)
