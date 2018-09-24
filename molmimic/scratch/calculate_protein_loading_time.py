import sys
sys.path.append("/data/draizene/molmimic")
import os

def calculate_time(pdb, chain, id):
    import time
    from molmimic.biopdbtools import Structure
    import numpy as np
    start = time.time()
    indices, _ = Structure.features_from_string(pdb, chain, id=id)
    total_time = time.time()-start
    
    num_atoms = indices.shape[0]
    
    min_coord = np.min(indices, axis=0)
    max_coord = np.max(indices, axis=0)
    dist = int(np.ceil(np.linalg.norm(max_coord-min_coord)))
    
    with open("pdbs/time_{}_{}.txt".format(pdb, chain), "w") as f:
        print >> f, "{}\t{}\t{}\t{}\t{}\t{}".format(pdb, chain, id, total_time, num_atoms, dist)

def load_ibis(ibis_data, course_grained=False):
    from molmimic.calculate_features import SwarmJob
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data, input_shape=(512,512,512))
    data = dataset.data
    job = SwarmJob("ibis_time")
    for i, row in data.iterrows():
        id = dataset.full_data.loc[(dataset.full_data["pdb"]==row["pdb"])&(dataset.full_data["chain"]==row["chain"])].iloc[0]["gi"]
        print "Running {}: {}.{}".format(id, row["pdb"], row["chain"])
        job += "/data/draizene/3dcnn-torch-py2 python {} {} {} {}\n".format(os.path.realpath(__file__), row["pdb"], row["chain"], id)
    job.run()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        load_ibis(sys.argv[1])
    elif len(sys.argv) == 4:
        calculate_time(*sys.argv[1:])
    else:
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain, id)")