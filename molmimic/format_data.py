import sys
from itertools import izip

import numpy as np

from Qsub import Qsub

def parse_cluster(cluster_num, num_clusters, members):
    print(cluster_num, "of", num_clusters, ":", members)

    protein_grid_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/results/molmimic/piface_v3/cluster_{}.npy".format(cluster_num)
    truth_grid_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/results/molmimic/piface_v3/cluster_{}_truth.npy".format(cluster_num)

    protein_grid = np.memmap(protein_grid_f, dtype='float32', mode='w+', shape=(2*len(members), 144,144,144,50))
    truth_grid = np.memmap(truth_grid_f, dtype='float32', mode='w+', shape=(2*len(members), 144,144,144, num_clusters*2+1))

    for protein_num, pdb_chains in enumerate(members):
        print(pdb_chains)
        pdb = pdb_chains[:4]
        chain1 = pdb_chains[-2]
        chain2 = pdb_chains[-1]

        try:
            structure = Structure.from_pdb(pdb, expand_to_p1=False)
        except IOError:
            continue

        structure = structure.extract_peptides()
        structure.orient_to_pai()

        protein1 = structure.extract_chain(chain1)
        protein2 = structure.extract_chain(chain2)

        binding_site1, binding_site2, distances = protein1.get_interface(protein2)

        for chain_num, protein, binding_site in ((0, protein1, binding_site1), (1, protein2, binding_site2)):
            for atom in protein.pdb_hierarchy.atoms():
                grid = protein1.get_grid_coord(atom)
                protein_grid[2*protein_num+chain_num, grid[0], grid[1], grid[2], :] = protein.get_features_for_atom(atom)
            for atom in binding_site.pdb_hierarchy.atoms():
                grid = protein1.get_grid_coord(atom)
                truth_grid[2*protein_num+chain_num, grid[0], grid[1], grid[2], cluster_num*2+chain_num-1] = 1

        # for binding_site in (binding_site1, binding_site2):
        #     min_coord = binding_site.pdb_hierarchy.atoms().extract_xyz().min()
        #     reset = binding_site.pdb_hierarchy.atoms().extract_xyz()-min_coord

        #     binding_site.pdb_hierarchy.atoms().set_xyz(reset)

        #     for atom in binding_site.pdb_hierarchy.atoms():
        #         grid_point = (np.floor(atom.xyz[0]), np.floor(atom.xyz[1]), np.floor(atom.xyz[2]))

        #         residue = binding_site.get_residue(atom, one_hot=True)
        #         charge = atom.charge_as_int()
        #         kd = binding_site.get_hydrophobicity(atom)
        #         biological = binding_site.get_hydrophobicity(atom, "biological")
        #         octanal = binding_site.get_hydrophobicity(atom, "octanal")

        #         rg_asa = binding_site.parent.get_accessible_surface_area(index=binding_site.parent_iseqs[atom.i_seq])

        #         try:
        #             grid[grid_point[0], grid_point[1], grid_point[2], :] = np.array(
        #                 residue+[charge, kd, biological, octanal, rg_asa])
        #         except IndexError:
        #             #Interface is too large!
        #             continue



def parse_full(pdb):
    pass


def parse_piface(num_clusters=10, filter=True, min_members=20, max_members=20):
    piface_file = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/PIface_01_24_2014.txt"

    if num_clusters is None:
        with open(piface_file) as f:
            num_clusters = sum((1 for line in f))

    with open(piface_file) as f:
        curr_cluster = 0
        for cluster_num, line in enumerate(f):
            if curr_cluster >= num_clusters: break

            if "3Q05" in line or "3ASN" in line: continue

            members = line.rstrip().rsplit("Members", 1)[-1].split()
            if filter and (len(members) < min_members or len(members) > max_members): continue

            qsub = Qsub("piface_{}".format(cluster_num), threads=1)
            qsub += "source /panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/phenix/phenix-dev-2890/phenix_env.sh\n"
            qsub += "iotbx.python {} {} {} {}\n".format(sys.argv[0], curr_cluster, num_clusters, " ".join(members))
            qsub.submit()
            curr_cluster += 1

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parse_piface()
    elif len(sys.argv) > 1:
        from pdbtools import Structure
        parse_cluster(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3:])
    else:
        print "?", sys.argv
