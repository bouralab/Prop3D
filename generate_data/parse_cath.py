def get_seqres():
    seqres = pd.read_table("cath-domain-boundaries-seqreschopping.txt", header=None, names=["cath_id", "ranges"])
    seqres.loc[:, ["pdb", "chain", "domNo"]] = seqres["cath_id"].apply(
        lambda pcd: pd.Series({"pdb":pcd[:4], "chain":pcd[4], "domNo":int(pcd[5:])}))
    return seqres

def group_cath(sfam_id):
    names = pd.read_table("cath-names.txt", header=None, names=["node", "cath_id", "name"])

    all_sdoms = pd.read_hdf(os.path.join(data_path_prefix, "PDB.h5"), "merged")
    sdoms_groups = all_sdoms.groupby("sfam_id")

    try:
        sdoms = sdoms_groups.get_group(sfam_id).copy()
    except KeyError:
        job.log("No structural domains for {}".format(cdd))
        return

    del sdoms_groups
    sdoms = sdoms.rename(columns={"pdbId":"pdb", "chnLett":"chain"})

    sdoms = pd.merge(sdoms, get_seqres(), on=["pdb", "chain", "domNo"], how="left")

    cath_node = pd.merge(sdoms, names, on="cath_id")["node"].iloc[0]

    
