from __future__ import print_function

import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/usr/share/pdb2pqr")

import os
import argparse
import time
from datetime import datetime
import subprocess

import torch
import sparseconvnet as scn

import numpy as np

def calculate_features(pdb, chain, id):
    from molmimic.biopdbtools import Structure
    inputSpatialSize = torch.LongTensor((264,264,264))
    struct = Structure(pdb, chain, id=id)
    for rotation in struct.rotate(1000):
        indices, data = struct.map_atoms_to_voxel_space()
        inputs = scn.InputBatch(3, inputSpatialSize)
        inputs.addSample()
        try:
            inputs.setLocations(torch.from_numpy(indices).long(), torch.from_numpy(data).float(), 0)
        except AssertionError:
            theta, phi, z = rotation[1:]
            min_coord = np.min(indices, axis=0)
            max_coord = np.max(indices, axis=0)
            dist = int(np.ceil(np.linalg.norm(max_coord-min_coord)))

            with open("{}_{}_{}.txt".format(pdb, chain, id), "a") as f:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    pdb, chain, id, theta, phi, z, dist, 
                    min_coord[0], min_coord[1], min_coord[2], 
                    max_coord[0], max_coord[1], max_coord[2]), file=f)
        del inputs
        del indices
        del data

def load_ibis(ibis_data):
    from molmimic.calculate_features import SwarmJob
    from molmimic.torch_model.torch_loader import IBISDataset
    dataset = IBISDataset(ibis_data, input_shape=(264,264,264))
    data = dataset.data
    job = SwarmJob("size_rots")
    for i, row in data.iterrows():
        id = dataset.full_data.loc[(dataset.full_data["pdb"]==row["pdb"])&(dataset.full_data["chain"]==row["chain"])].iloc[0]["gi"]
        print("Running {}: {}.{}".format(id, row["pdb"], row["chain"]))
        job += "/data/draizene/3dcnn-torch-py2 python {} {} {} {}\n".format(os.path.realpath(__file__), row["pdb"], row["chain"], id)
    job.run()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        load_ibis(sys.argv[1])
    elif len(sys.argv) == 4:
        calculate_features(*sys.argv[1:])
    else:
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain, resi)")
