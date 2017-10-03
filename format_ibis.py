import sys

def parse_ibis_centroid(index, pdb, chain, resi):
    """
    """
    import h5py
    f = h5py.File('{}_{}_{}.h5'.format(index, pdb, chain), 'w')
    data_grid = f.create_dataset("data", (96, 96, 96, 50), fillvalue=0, compression='gzip', compression_opts=9)
    truth_grid = f.create_dataset("truth", (96, 96, 96, 1), fillvalue=0, compression='gzip', compression_opts=9)
    print index, pdb, chain, resi
    try:
        structure = Structure.from_pdb(pdb, expand_to_p1=False)
    except IOError:
        f.close()
        return

    structure = structure.extract_peptides()
    structure = structure.extract_chain(chain[0])
    structure.orient_to_pai()
    sele = "resseq {}".format(resi.replace(",", " or resseq "))
    binding_site = structure.extract_selection_from_string(sele)
    
    for atom in binding_site.pdb_hierarchy.atoms():
        grid = binding_site.get_grid_coord(atom)
        data_grid[grid[0], grid[1], grid[2], :] = structure.get_features_for_atom(binding_site.parent_iseqs[atom.i_seq], i_seq=True)
        truth_grid[grid[0], grid[1], grid[2], 0] = 1
    
    f.close()

def concatenate(h5_files_f):
    import h5py

    with open(h5_files_f) as f:
        num_binding_sites = sum(1 for _ in f)

    output_file = h5py.File('ibis_binding_sites.h5', 'w')
    data_grid = output_file.create_dataset("data", (num_binding_sites, 96, 96, 96, 50), fillvalue=0., compression='gzip', compression_opts=9)
    truth_grid = output_file.create_dataset("truth", (num_binding_sites, 96, 96, 96, 1), fillvalue=0., compression='gzip', compression_opts=9)

    with open(h5_files_f) as files:
        for i, file in enumerate(files):
            print "FILE:", file.rstrip()
            file = h5py.File(file.rstrip(), 'r')
            data_grid[i, ...] = file["data"]
            truth_grid[i, ...] = file["truth"]
            file.close()

    output_file.close()

def parse_ibis_centroids(tax_glob_group="A_eukaryota", num_represtatives=20):
    import pandas as pd
    from Qsub import Qsub

    observations_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/unique_obs_bs-8.tab"
    resfaces_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/unique_molresface_resi_resn_NR_all.tab"

    observations = pd.read_table(observations_f, header=0, usecols=(0,2))
    resfaces = pd.read_table(resfaces_f, header=0)
    data = pd.merge(resfaces, observations, left_on="unique_obs_int", right_on="uniq_obs_int_id")
    data = data.loc[(data["tax_glob_group"] == tax_glob_group) & (data["n"] >= num_represtatives)]
    
    h5_data_files = open("h5_data_files.txt", "w")
    job_ids = []
    for i, row in data.iterrows():
        qsub = Qsub("ibis_{}".format(i), threads=1)
        qsub += "source /panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/phenix/phenix-dev-2890/phenix_env.sh\n"
        qsub += "iotbx.python {} {} {} {} {}\n".format(sys.argv[0], i, row["pdb"], row["chain"], row["resi"])
        jid = qsub.submit()
        job_ids.append(jid)
        print >> h5_data_files, "{}_{}_{}.h5".format(i, row["pdb"], row["chain"])
    h5_data_files.close()

    qsub = Qsub("ibis_merge".format(i), threads=1)
    qsub += "python {} {}\n".format(sys.argv[0], "h5_data_files.txt")
    qsub.submit(hold_jid=job_ids)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parse_ibis_centroids()
    elif len(sys.argv) == 5:
        from pdbtools import Structure
        parse_ibis_centroid(*sys.argv[1:])
    elif len(sys.argv) == 2:
        concatenate(sys.argv[1])
    else:
        print "?", sys.argv


