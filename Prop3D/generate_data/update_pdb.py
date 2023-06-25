import time
import urllib.request
from datetime import datetime

import h5pyd
import requests
import pandas as pd
from biotite.database.rcsb import FieldQuery, CompositeQuery, search, count, IdentityGrouping

from toil.realtimeLogger import RealtimeLogger


from Prop3D.generate_data.create_data_splits import split_dataset_at_level

def update_pdbs(job, full_pdb_h5):
    with h5pyd.open(full_pdb_h5, use_cache=False) as f:
        try:
            last_updated = f.attrs['last_updated']
        except KeyError:
            #No need to update, start from scratch
            return None
    
    last_updated_query = FieldQuery("rcsb_accession_info.deposit_date", greater_or_equal=last_updated)
    today = datetime.today().strftime('%Y-%m-%d')
    today_query = FieldQuery("rcsb_accession_info.deposit_date", less_or_equal=today)
    experimental_query = FieldQuery("rcsb_entry_info.structure_determination_methodology", exact_match="experimental")
    #entity_query = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
    query = CompositeQuery(prev_query, today_query, experimental_query)

    new_pdbs = search(query, return_type="polymer_entity")

    last_updated_revised_query = FieldQuery("rcsb_accession_info.revision_date", greater_or_equal=last_updated)
    today_revised_query = FieldQuery("rcsb_accession_info.revision_date", less_or_equal=today)
    revised_pdbs = search(last_updated_revised_query & today_revised_query, return_type="polymer_entity")

    with h5pyd.open(full_pdb_h5, use_cache=False) as f:
        for rpdb in revised_pdbs:
            del f[f"domain/{rpdb}"]

    return set(new_pdbs)&set(revised_pdbs)

def get_all_pdbs(job, full_pdb_h5, update=False):
    pdbs = None
    if update:
        pdbs = update_pdbs(job, full_pdb_h5)

    if pdbs is None:
        experimental_query = FieldQuery("rcsb_entry_info.structure_determination_methodology", exact_match="experimental")
        #entity_query = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
        #query = CompositeQuery(experimental_query, entity_query)
        pdbs = search(experimental_query, return_type="polymer_entity")
    
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group("domains")
        for pdb in pdbs:
            store.require_group(f"domains/{pdb}")
    
    job.addFollowOnJobFn(create_splits_for_levels, full_pdb_h5)
    
def get_custom_pdbs(job, pdbs, full_pdb_h5, update=False):
    if not isinstance(pdbs, (list, tuple)):
        pdbs = [pdbs]
    
    #Try polymer_entity first:
    entity_query = FieldQuery("rcsb_entry_info.polymer_entity", is_in=pdbs)
    entities = search(entity_query, return_type="polymer_entity")

    if len(entities) != len(pdbs):
        #Then try old chains names (polymer_instance or asym id)
        chain_query = FieldQuery("rcsb_entry_info.polymer_instance", is_in=pdbs)
        chains = search(entity_query, return_type="polymer_entity")

        if len(chains) != len(pdbs):
            #Not sure why, but might've used just the 4-letter PDB codes?
            entry_query = FieldQuery("rcsb_entry_container_identifiers.entry_id", is_in=pdbs)
            entries = search(entity_query, return_type="polymer_entity")

            if len(entries) < len(pdbs):
                #Maybe combined the three types?
                entities = set(entities).intersection(chains).intersection(entries)
                RealtimeLogger.info(f"Unable to match all chains (polymer_instance or asym id) to polymer_entity. Using {len(entities)/len(pdbs)}: {entities}")
            else:
                entities = entries
        else:
            entities = chains
            
    assert len(entities) > 0

    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group("domains")
        for pdb in entities:
            store.require_group(f"domains/{pdb}")
        
    job.addFollowOnJobFn(create_splits_for_levels, full_pdb_h5, custom_pdbs=entities)

def create_representatives(job, clusters, full_pdb_h5):
    key = "representatives"
    representatives = clusters["representatives"].drop_duplicates()
    missing_domains = []

    with h5py.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        group = store.require_group(key)

        for domain in representatives:
            try:
                group[domain] = store[f"domains/{domain}'"]
            except KeyError:
                missing_domains.append(domain)

        if len(missing_domains) > 0:
            store[key].attrs["missing_domains"] = missing_domains
        store[key].attrs["total_domains"] = len(representatives)

def update_clusters(job, full_pdb_h5, pct_id=30, update=False, custom_pdbs=None):
    assert pct_id in [30, 40, 50, 70, 90, 95, 100]
    
    if False:
        url = f"https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{pct_id}.txt"
        conn = urllib.request.urlopen(url, timeout=30)

        if update:
            time_struct = time.strptime(conn.headers['last-modified'], '%a, %d %b %Y %H:%M:%S %Z')
            with h5pyd.open(full_pdb_h5, use_cache=False) as f:
                try:
                    last_updated = time.strptime(f.attrs['last_updated'], '%Y-%m-%d')
                except KeyError:
                    #No need to update, start from scratch
                    return None
            
            if time_struct<=last_updated:
                return
            
            clusters = None
            with requests.get(url, stream=True) as r:
                r.raw.decode_content = True
                for i, line in enumerate(r.content):
                    cluster = line.rstrip().split()
                    cluster = pd.DataFrame({"cath_domain":cluster, "cluster_name":i})
                    if clusters is not None:
                        clusters = pd.concat((clusters, cluster))
                    else:
                        clusters = cluster
            clusters.reset_index()
    
    query = rcsb.FieldQuery("rcsb_entry_info.structure_determination_methodology", exact_match="experimental")
    # entity_query = rcsb.FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
    # query = experimental_query & entity_query
    if custom_pdbs is not None:
        entity_query = FieldQuery("rcsb_entry_info.polymer_entity", is_in=custom_pdbs)
        query = query & entity_query
    pdbs = search(query, group_by=IdentityGrouping(30), return_groups=True, return_type="polymer_entity")
    clusters = pd.DataFrame([
        (pdb,group_name,group[0]) for group_name, group in pdbs.items() for pdb in group], 
        columns=["cath_domain", "cluster_name", "representative"])

    job.addFollowOnJob(split_dataset_at_level, full_pdb_h5, "", clusters, pct_id, str(pct_id),
        split_size={"train":0.8, "validation":0.1, "test":0.1})
    
    if pct_id == 30:
        job.addFollowOnJob(create_representatives, clusters, full_pdb_h5)

def create_splits_for_levels(job, full_pdb_h5, custom_pdbs=None):
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group(f"data_splits")

    for level in [30, 40, 50, 70, 90, 95, 100]:
        job.addChildJobFn(update_clusters, full_pdb_h5, level, custom_pdbs=custom_pdbs)
