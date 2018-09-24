from biopdbtools import Structure
from visualize import plot_volume

pdbs_info = [
        #("4WJ5", "A", "7,8,9,11,22,24,26,31,34,47,95,96,97,112,114,115,116,119,121"), 
        ("3WKJ", "H", "80,84,87,96,97,100,104,105"), 
        #("3WOG", "B", "150,165,172,174,178,182,185,186,187,188,189,190,191,192,193")
        ]

pdbs = [Structure(id.lower(), chain) for id, chain, _ in pdbs_info]
pdbs = [s.extract_chain(s.chain) for s in pdbs]

binding_sites = [s.align_seq_to_struc(i[-1], return_residue=True) for s,i in zip(pdbs, pdbs_info)]


for s in pdbs:
    s.orient_to_pai()

volumes = [s.create_full_volume() for s in pdbs]

binding_sites = [s.align_seq_to_struc(i[-1], return_residue=True) for s,i in zip(pdbs, pdbs_info)]
binding_site_grid = [s.get_features(residue_list=b, return_data=False) for s,b in zip(pdbs, binding_sites)]

plot_volume(volumes[0])

plot_volume(volumes[0], binding_site=binding_site_grid[0])



