import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/data/draizene/3DUnetCNN")
sys.path.append("/data/draizene/seaborn")

import traceback

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from molmimic.keras_model.pdb_generator import IBISGenerator
from molmimic.biopdbtools import Structure

def calc_size(ibis_data):
    max_dist = 0
    min_dist = 100

    data = IBISGenerator(ibis_data, input_shape=(96,96,96,59))

    with open("protein_distances.txt", "w") as data_file:
        print >> data_file, "pdb\tchain\tid\tmin\tmax\tdist"
        for i, row in data.data.iterrows():
            print "Running {}:{}.{} ({} of {})".format(row["unique_obs_int"], row["pdb"], row["chain"], i+1, data.data.shape[0])

            try:
                structure = Structure(row["pdb"], row["chain"], id=row["unique_obs_int"])
            except (KeyboardInterrupt, SystemExit) as e:
                raise
            except:
                trace = traceback.format_exc()
                print "Error:", trace
                continue

            coords = structure.get_coords()
            min_coord = np.min(coords, axis=0)
            max_coord = np.max(coords, axis=0)
            dist = int(np.ceil(np.linalg.norm(max_coord-min_coord)))

            print >> data_file, "{}\t{}\t{}\t{}\t{}\t{}".format(
                row["pdb"],
                row["chain"],
                row["unique_obs_int"],
                min_coord.tolist(),
                max_coord.tolist(),
                dist)

            max_dist = max(max_dist, dist)
            min_dist = min(min_dist, dist)

        print >> data_file, "#Max distance: {}".format(max_dist)
        print >> data_file, "#Min distance: {}".format(min_dist)
        print "Max distance: {}".format(max_dist)
        print "Min distance: {}".format(min_dist)

def plot_sizes(size_data):
    max_size = 0
    sizes = []
    pdbs = []
    with open(size_data) as f:
        f.next()
        for line in f:
            if line.startswith("#"):
                #Max size
                max_size = int(line.rstrip().split(" ")[-1])
                break
    
            pdb, chain, id, min_coord, max_coord, dist = line.rstrip().split("\t")
        
            if int(dist) >= 249:
                print pdb, chain, id, dist

            if (pdb, chain) in pdbs:
                continue

            sizes.append(int(dist))
            pdbs.append((pdb, chain))

    sns.set(style="white", palette="muted", color_codes=True)
    f, axes = plt.subplots(1, 1, figsize=(12, 12))
    
    sns.distplot(sizes, hist=False, color="g", kde_kws={"shade": True}, ax=axes)

    plt.setp(axes, yticks=[])
    #plt.tight_layout()
    plt.xlabel('Distance from max to min Cartesian coordinate', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        calc_size(sys.argv[1])
