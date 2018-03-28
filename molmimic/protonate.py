import sys
sys.path.append("/data/draizene/molmimic")

import os
import argparse
import time
from datetime import datetime
from itertools import groupby

def load_ibis(ibis_data, minimize=True):
    from molmimic.calculate_features import SwarmJob
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data)
    data = dataset.data
    job = SwarmJob("pdb2pqr", cpus=1)

    key = lambda p: p[1:3]
    pdbs = sorted(set(data["pdb"].values.tolist()), key=key)

    i = 0
    for pdb_divided, pdb_group in groupby(pdbs, key=key):
        pdb_path = "/data/draizene/molmimic/pdb2pqr_structures/pdbs/{}".format(pdb_divided.lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
        for pdb in pdb_group:
            print "Running {}/{}: {}".format(i+1, len(pdbs), pdb)

            #unzip, tidy, and remove hetatms
            job += "gunzip -c /pdb/pdb/{0}/pdb{1}.ent.gz | python /data/draizene/pdb-tools/pdb_tidy.py | python /data/draizene/pdb-tools/pdb_striphet.py > {2}/{1}.pdb ; ".format(pdb_divided.lower(), pdb.lower(), pdb_path)

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

            if minimize:
                job += "minimize -s $base.pqr.pdb "
                job += "-run:min_type lbfgs_armijo_nonmonotone -run:min_tolerance 0.001 -ignore_zero_occupancy false "
                job += "-out:file:scorefile $base.sc "
                job += "-out:path:pdb {0} -out:path:score {0} ; ".format(pdb_path)


                job += "python /data/draizene/pdb-tools/pdb_stripheader.py $base.pqr_0001.pdb > $base.min.pdb; "
                job += "mv $base.pqr.pdb.propka $base.propka; "
                #job += "mv $f $f.unminimized.pdb; ".format(pdb_path, pdb.lower())

                #
                #job += "mv $f.pqr_0001.pdb $base.pdb; "

            job += "done; \n"


            i += 1
    job.run()

if __name__ == "__main__":
    load_ibis(sys.argv[1])
