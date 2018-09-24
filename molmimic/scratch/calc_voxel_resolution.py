import os, sys, errno
sys.path.append("/data/draizene/molmimic")

def atom_distances(structure):
    import numpy as np
    from molmimic.biopdbtools import Structure
    from sklearn.metrics.pairwise import euclidean_distances

    if isinstance(structure, tuple):
        from molmimic.biopdbtools import Structure
        structure = Structure(structure[0], structure[1])

    coords = s.get_coords(exclude_atoms=["N", "C", "O"])
    distances = euclidean_distances(coords, squared=True)
    mask = np.ones(distances.shape, dtype=bool)
    np.fill_diagonal(mask, 0)
    min_dist = np.sqrt(distances[mask].min())

    directory = "atom_distances"

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    with open(os.path.join(directory, "{}_{}_min_dist.txt".format(structure.pdb, structure.chain)), "w") as f:
        print >> f, min_dist

def atom_density(structure):
    import numpy as np
    if isinstance(structure, tuple):
        from molmimic.biopdbtools import Structure
        structure = Structure(structure[0], structure[1])
    coords = structure.get_coords()
    min_coord = np.min(coords, axis=0)
    max_coord = np.max(coords, axis=0)
    volume = np.prod(max_coord-min_coord)
    density = len(coords)/float(volume)

    directory = "atom_densities"

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    with open(os.path.join("{}_{}_density.txt".format(structure.pdb, structure.chain)), "w") as f:
        print >> f, density

def atom_density_by_sliding_box(structure, volume=5, stride=1, pad=True):
    import numpy as np
    from Bio.PDB.NeighborSearch import NeighborSearch

    if isinstance(structure, tuple):
        from molmimic.biopdbtools import Structure
        structure = Structure(structure[0], structure[1])

    kdtree = NeighborSearch(self.get_atoms(exclude_atoms=["N", "C", "O"]))

    coords = structure.get_coords()
    min_coord = np.min(coords, axis=0)
    max_coord = np.max(coords, axis=0)
    volume = np.prod(max_coord-min_coord)

    pad_min = stride/2
    pad_max = int(np.ciel(stride/2.))

    total_volume = volume*volume*volume

    densities = []
    for i, x in enumerate(xrange(int(min_coord[0])-pad_min, int(max_coord[0])+pad_max, stride)):
        for j, y in enumerate(xrange(int(min_coord[1])-pad_min, int(max_coord[1])+pad_max, stride)):
            for k, z in enumerate(xrange(int(min_coord[2])-pad_min, int(max_coord[2])+pad_max, stride)):
                neighbors = kdtree.search(((2*x+1)/2., (2*y+1)/2., (2*z+1)/2.), radius=volume/2.)
                densities.append(len(neighbors)/float(total_volume))

    density = np.mean(densities)
    std = np.std(densities)


    directory = "atom_densities2"

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    with open(os.path.join(directory, "{}_{}_density_{}.txt".format(structure.pdb, structure.chain, volume)), "w") as f:
        print >> f, "{}\t{}".format(density, std)

def atom_collisions(structure):
    from collections import Counter
    import numpy as np

    if isinstance(structure, tuple):
        from molmimic.biopdbtools import Structure
        structure = Structure(structure[0], structure[1])

    directory = "atom_collisions"

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for grid_size in xrange(25, 255, 5):
        grid_size /= 100.
        structure.set_voxel_size(grid_size)
        occupancy = Counter()
        for atom in structure.get_atoms(exclude_atoms=["N", "C", "O"]):
            grid = tuple(structure.get_grid_coord(atom).tolist())
            occupancy[grid] += 1
        with open(os.path.join(direcotry, "{}_{}_collisions_{}.tsv".format(structure.pdb, structure.chain, grid_size)), "w") as f:
            print f.name
            for grid, collisions in occupancy.iteritems():
                print >> f, "{}\t{}\t{}\t{}".format(grid[0], grid[1], grid[2], collisions)

def run(pdb, chain, density=True, collisions=False, atom_distances=False):
    from molmimic.biopdbtools import Structure
    structure = Structure(pdb, chain)

    if density:
        atom_density_by_sliding_box(structure)

    if collisions:
        atom_collisions(structure)

def plot(density=True, collisions=True, atom_distances=True):
    import glob
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    directory = "plots"

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    if density:
        densities = []
        for f in glob.glob("/data/draizene/molmimic/atom_densities/*_density.txt"):
            with open(f) as f_:
                densities.append(float(f_.read().rstrip()))
        sns.distplot(densities, hist=False, axlabel="Densities")
        plt.savefig(os.path.join(directory, "densities.pdf"))

    if collisions:
        collision_data = pd.DataFrame()
        for grid_size in xrange(25, 255, 5):
            grid_size /= 100.
            with open("/data/draizene/molmimic/atom_collisions/{}.txt".format(grid_size), "w") as out:
                for f in glob.glob("/data/draizene/molmimic/atom_collisions/*_collisions_{}.tsv".format(grid_size)):
                    with open(f) as f_:
                        for line in f_:
                            print >> out, line.rstrip().split("\t")[-1]
            print grid_size

        collision_data = pd.DataFrame()
        for grid_size in xrange(25, 255, 5):
            grid_size /= 100.
            data = pd.read_table("/data/draizene/molmimic/atom_collisions/{}.txt".format(grid_size), names="Occupancy")
            sns.distplot(data, hist=False, axlabel="Occupancy")
            plt.savefig(os.path.join(directory, "collisions_{}.pdf".format(grid_size)))
            data["Grid Size"] = grid_size
            collision_data = pd.concat((collision_data, data))

        sns.swarmplot(x="Grid Size", y="Occupancy", data=collision_data)
        plt.savefig("collisions.pdf".format(grid_size))

    if atom_distances:
        distances = []
        for f in glob.glob("/data/draizene/molmimic/atom_distances_dir/*_min_dist.txt"):
            with open(f) as f_:
                distances.append(float(f_.read().rstrip()))
        sns.distplot(distances, hist=False, axlabel="Distances")
        plt.savefig(os.path.join(directory, "distances.pdf"))

def voxel_resolution2(pdb, chain):
    import numpy as np
    from molmimic.biopdbtools import Structure
    from sklearn.metrics.pairwise import euclidean_distances
    s = Structure(pdb, chain)
    coords = s.get_coords(exclude_atoms=["N", "C", "O"])
    distances = euclidean_distances(atoms)
    mask = np.ones(distances.shape, dtype=bool)
    np.fill_diagonal(mask, 0)
    min_dist = distances[mask].min()

    with open("{}_{}_distances.txt".format(pdb, chain), "w") as f:
        print >> f, min_dist

def min_voxel_resolution():
    import glob
    min_dist = 100
    for fname in glob.glob("*_min_dist.txt"):
        with open(fname) as f:
            dist = float(f.read().rstrip())
            min_dist = min(min_dist, dist)
    print "Min distance is", min_dist
    with open("global_min_dist.txt", "w") as f:
        print >> f, min_dist

def load_ibis(ibis_data):
    import os
    from molmimic.calculate_features import SwarmJob
    from molmimic.torch_model.torch_loader import IBISDataset
    data = IBISDataset(ibis_data, input_shape=(512,512,512))
    data = data.data[['pdb', 'chain']].drop_duplicates()
    jobs = []
    for i, row in data.iterrows():
        print "Running {}.{}".format(row["pdb"], row["chain"])
        job = SwarmJob("dist_{}.{}".format(row["pdb"], row["chain"]))
        job += "/data/draizene/3dcnn-torch python {} {} {}\n".format(os.path.realpath(__file__), row["pdb"], row["chain"])
        jid = job.submit()
        jobs.append(jid)

    print "Running Reduce"
    job = SwarmJob("reduce")
    job += "/data/draizene/3dcnn python-torch {} reduce\n".format(os.path.realpath(__file__), row["pdb"], row["chain"])
    job.submit(hold_jid=jobs)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        if sys.argv[1] == "reduce":
            plot()
        else:
            load_ibis(sys.argv[1])
    elif len(sys.argv) == 3:
        run(*sys.argv[1:])
    else:
        raise RuntimeError("Must be path to ibis luca file or (pdb, chain)")
