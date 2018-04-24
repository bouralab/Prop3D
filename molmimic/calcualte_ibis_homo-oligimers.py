import pandas as pd

cdds = ['ABC_CysA_sulfate_importer','secoisolariciresinol-DH_like_SDR_c','Cygb','PBP2_iGluR_Kainate_GluR6','Ngb','HGbI_like','IgC_MHC_I_alpha3','STKc_Aurora-A','CuRO_1_CuNIR','3beta-17beta-HSD_like_SDR_c','IgV_CTLA-4','IgV_TCR_delta','MDH-like_SDR_c','PBP2_iGluR_NMDA_Nr1','STKc_STK10','IgC_MHC_II_alpha','V-set','PTKc_Lck_Blk','ABC_MTABC3_MDL1_MDL2','ABC_MJ0796_LolCDE_FtsE','SDR_c5','Ig_CD8_alpha','STKc_p70S6K','SDR_c1','ABC_ATPase','Ig1_Nectin-1_like','ABC_cobalt_CbiO_domain1','BKR_3_SDR_c','STKc_PDK1','PBP2_iGluR_Kainate','TCP1_zeta','IgC_L','proteasome_beta_type_6','proteasome_beta_type_7','proteasome_beta_type_4','proteasome_beta_type_5','proteasome_beta_type_2','proteasome_beta_type_3','SDR_c11','proteasome_beta_type_1','SDR','STKc_PLK1','STKc_CaMKI_delta','PTKc_FAK','PTKc_InsR','STKc_p38alpha','protease_HslV','ABCF_EF-3','Ig2_VCAM-1','NR_LBD_RAR','Ig_CEACAM_D1_like','PTKc_EGFR','HBDH_SDR_c','STKc_LKB1','STKc_ASK','STKc_MAPK','NR_LBD','Ig_pIgR','NR_LBD_HNF4_like','CuRO_1_CuNIR_like','IgV_TCR_gamma','STKc_SnRK2-3','SDR_c','proteasome_alpha_type_4','STKc_IRE1','NR_LBD_PXR_like','NR_LBD_LXR','PR_SDR_c','STKc_MELK','TCP1_delta','PKc_MEK2','Periplasmic_Binding_Protein_Type_2','PKc_MEK1','NR_LBD_ER','STKc_OSR1_SPAK','PBP2_iGluR_Kainate_GluR5','IgC_TCR_delta','ABCC_Glucan_exporter_like','Ig4_SCFR','IGc1','Ga5DH-like_SDR_c','Ig1_LILRB1_like','PBP2_Arg_STM4351','STKc_PhKG2','HHT1','Ig1_Nectin-3_like','proteasome_alpha_type_2','IgC_CH1','PBP2_Dsm1740','STKc_Raf','Ig_SLAM-CD84_like_N','PBP2_GluR0','Ig_P0','STKc_PAK1','GroEL','AAA_2','haloalcohol_DH_SDR_c-like','XR_like_SDR_c','Globin_like','TCP1_eta','TER_DECR_SDR_a','3alpha_HSD_SDR_c','IgC','Hb','PTKc_Btk_Bmx','PTKc_RET','PTKc_Jak1_rpt2','Ig1_Necl-1','STKc_MAP3K-like','TCP1_epsilon','ABC_FeS_Assembly','Ig_MOG_like','BKR_like_SDR_like','ABC_ModC_like','ABC_SMC_barmotin','STKc_Twitchin_like','ChcA_like_SDR_c','meso-BDH-like_SDR_c','ENR_SDR','PBP2_Arg_3','IgV_CD8_beta','STKc_PknB_like','STKc_GAK','STKc_LIMK2','ABCC_MsbA','IGv','Ig2_FGFR','BphB-like_SDR_c','NR_LBD_MR','STKc_IRAK4','CuRO_1_FV_like','CR_SDR_c','H4','DHRS6_like_SDR_c','Ig2_Nectin-2_like','STKc_GSK3','ABC_SMC1_euk','TCP1_alpha','STKc_Chk1','Cupredoxin','Ig1_Nectin-4_like','A3DFK9-like_SDR_c','cpn60','proteasome_alpha_archeal','IgV_L_kappa','RDH_SDR_c','cyclohexanol_reductase_SDR_c','PTKc_TrkA','ABC_HisP_GlnQ','IgV_TCR_alpha_like','Ig_2','ABC_MalK_N','PBP2_iGluR_NMDA_Nr2','Ig1_PVR_like','Ig3_Contactin_like','HSD10-like_SDR_c','Hb-beta_like','STKc_CAMK','Ig_Myomesin_like_C','PKc_like','7_alpha_HSDH_SDR_c','STKc_TBK1','ABCC_Hemolysin','STKc_NIK','proteasome_protease_HslV','ABC_Carb_Solutes_like','PTKc_FGFR2','PTKc_FGFR1','CuRO_2_CuNIR','PTKc_Itk','NR_LBD_EcR','Ig_CEACAM_D1','STKc_AMPK_alpha','STKc_CK1_delta_epsilon','SDR_c12','STKc_SLK','HetN_like_SDR_c','THN_reductase-like_SDR_c','STKc_EIF2AK2_PKR','PTKc_Ror2','TR_SDR_c','STKc_CaMKII','Ig_CD3_epsilon','IgC_beta2m','CAD_SDR_c','ABCC_MRP_Like','ABC_YhbG','CuRO_1_2DMCO_NIR_like','IgV_TCR_beta','proteasome_alpha_type_5','CuRO_D2_2dMcoN_like','proteasome_alpha_type_7','proteasome_alpha_type_6','proteasome_alpha_type_1','STKc_CK2_alpha','proteasome_alpha_type_3','PK_STRAD_alpha','CBFD_NFYB_HMF','IgC_TCR_alpha','PBP2_PheC','PKc_Mps1','CuRO_1_2DMCO_NIR_like_2','STKc_PKB_beta','BKR_SDR_c','Ig2_FcgammaR_like','H2A','H2B','AAA','STKc_MARK','Mgc4172-like_SDR_c','PBP2_iGluR_AMPA_GluR4','chaperonin_like','Cu-oxidase_2','Ig_Titin_like','KDSR-like_SDR_c','IgC_TCR_gamma','IgC_CH3','IgC_CH2','TCP1_theta','proteasome_beta_archeal','NR_LBD_RXR_like','IgC_CH4','STKc_MAPK4_6','STKc_PhKG1','PKR_SDR_c','PBP2_iGluR_delta_2','class1_nsHb_like','STKc_CDK6','IgC_MHC_II_beta','STKc_IKK_beta','GlcDH_SDR_c','ABCC_CFTR1','ig','Mb_like','NR_LBD_TR2_like','Hb-alpha_like','Ig_FcalphaRI','STKc_EIF2AK4_GCN2_rpt2','Ig1_MRC-OX-2_like','STKc_RSK1_C','PBP2_iGluR_AMPA','TCP1_beta','IgV','STKc_SnRK3','Ntn_hydrolase','STKc_DAPK2','STKc_DAPK3','TCP1_gamma','STKc_DAPK1','IgC_CH2_IgE','NR_LBD_PPAR','CuRO_D1_2dMcoN_like','STKc_PRP4','CuRO_1_FVIII_like','STKc_Chk2','Ig1_CD4','SDR_c10','IgC_TCR_beta','STKc_ACVR1_ALK1','IgV_H','ABC_Rad50','DH-DHB-DH_SDR_c','NR_LBD_ERR','CuRO_4_FVIII_like','STKc_PKA','STKc_MAPK15-like','PTKc_Ack_like','IgV_L_lambda','STKc_ERK1_2_like','STKc_MAPKAPK2','ABC_Iron-Siderophores_B12_Hemin','SPR-like_SDR_c','Ig']

def get_homololigimers(cluster_size=500):
	tbl1 = pd.read_table("/data/draizene/molmimic/molmimic/data/unique_molresface_resi_resn_NR_all.tab.txt")
	tbl2 = pd.read_table("/data/draizene/molmimic/molmimic/data/unique_obs_bs-8.tab.txt")
	tbl2.columns = [u'unique_obs_int', u'sample_sdi', u'n', u'n_mmdbs', u'n_nr_mmdbs',
            u'cdds', u'sfs_ids', u'cdds_B', u'sfs_ids_B', u'single', u'self_lca',
            u'lca_level', u'lca_tax', u'lca_name', u'lca_depth']

	ibis_luca = pd.merge(tbl1, tbl2, on="unique_obs_int")
	large_families = ibis_luca[ibis_luca["n"]>cluster_size]

	homooligimers = [row for i, row in large_families.iterrows() if \
		not isinstance(row["cdds"], float) and not isinstance(row["cdds_B"], float) \
		and len(set(row["cdds"].split(",")).intersection(set(row["cdds_B"].split(",")))) > 0]

	homooligimers = pd.DataFrame(homooligimers)
	cdds = set([cdd for row in homooligimers["cdds"] for cdd in row.split(",")])
	return cdds

def process_ibis_down():
	full_data = "/data/draizene/molmimic/molmimic/data/ibis_full.tsv"
	ibis_full = pd.read_table(full_data)

	multimers = ibis_full[ibis_full["CDD"].isin(cdds)]

	pdb_chains = multimers.groupby(["pdb", "chain"])

	pdbs_with_one_family = pd.concat((grp for name, grp in pdb_chains if grp["CDD"].drop_duplicates().shape[0]==1))