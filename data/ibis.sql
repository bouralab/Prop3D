with t2 as 
	(
		select i2.*, count(i2.inf_id) OVER(partition by i2.clust_inf_id) as totcnt
    	from [Intrac1].[dbo].[Inferred1] i2
    	where i2.mol_sdi_id = id
    ),
    t1 as 
    (
    	select i2.*, o.int_sequence, o.int_cd_pssmid,
               o.int_superfam_acc, o.int_gi, o.mmdb_id, o.mol_mol_id, o.int_mol_id, o.mol_cd_pssmid,
               s.cid, s.sidname, s.cidname, s.exclude, s.mesh, s.pharmaction, s.active_assays, s.assays,
               dbo.fCheckAnnotType(o.mol_cd_cur_annot,i2.nbr_type) as cur_annot
        from Intrac1.dbo.Inferred1 i2 INNER JOIN Intrac.dbo.ObsInt o
        on i2.nbr_obs_int_id = o.obs_int_id\n"
        LEFT OUTER JOIN Intrac.dbo.SidCidInfo s on o.int_sid = s.sid
    )

select inf_id,mol_sdi_id,vast_id,nbr_sdi_id,nbr_obs_int_id,nbr_type,
       nbr_score,nbr_taxid,nbr_cd_from,nbr_cd_to,nbr_superfam_id,int_sid,
       int_superfam_id,clust_id,is_clust_rep,face_seq,
       int_superfam_acc,int_gi,cid,sidname,cidname,conf_idx,mmdb_id,exclude,pcnt_id,
       contact,conser_idx,int_sequence,mesh,pharmaction,cur_annot,
       comb_scr,mol_mol_id,int_mol_id,face_align,align_colors,consensus,clust_inf_id,
       curcnt,nrxtal,active_assays,assays,pssm_score,totcnt,1,
       int_cd_pssmid,vast_rmsd,mol_cd_pssmid
from 
	(
    	select i.*, c.comb_scr, c.conf_idx, c.pcnt_id, c.contact, c.conser_idx,
               c.align_colors, c.consensus, c.pssm_score
        from t1 i INNER JOIN Intrac.dbo.ObsInt o
        on i.nbr_obs_int_id = o.obs_int_id
        INNER JOIN Intrac1.dbo.ClustRep1 c on i.clust_inf_id = c.inf_id
    ) dt
order by nbr_type, pcnt_id desc, clust_id,"
         is_clust_rep desc, case vast_id when 0 then 0 else 1 end, nbr_score desc

