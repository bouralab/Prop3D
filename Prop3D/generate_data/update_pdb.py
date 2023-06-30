import time
import urllib.request
from datetime import datetime

import h5pyd
import requests
import pandas as pd
from biotite.database.rcsb import FieldQuery, CompositeQuery, search, count, IdentityGrouping

from Prop3D.util.toil import map_job
from toil.realtimeLogger import RealtimeLogger


from Prop3D.generate_data.create_data_splits import split_dataset_at_level

def chain2entity_old(pdbs):
    RealtimeLogger.info("Converting all chain ids to entity ids")
    values = []
    for pdb in pdbs:
        chain_query = FieldQuery("rcsb_polymer_entity_instance_container_identifiers.rcsb_id", exact_match=pdb)
        chains_result = [(pdb, entity) for entity in search(chain_query, return_type="polymer_instance")]
        values += chains_result
    values = pd.DataFrame(values, columns=["cath_domain", "entity_id"])
    RealtimeLogger.info("Converted all chain ids to entity ids")
    return values

def _chain2entity(pdbs, attempts=1):
    query = f"""{{
    polymer_entity_instances(instance_ids: ["{'", "'.join(pdbs)}"]) {{
    rcsb_id
    rcsb_polymer_entity_instance_container_identifiers {{
        entry_id
        entity_id
    }}
    }}
}}"""
    query = query.replace("\n", "")
    with requests.get(f"https://data.rcsb.org/graphql?query={query}") as r:
        try:
            results = r.json()["data"]["polymer_entity_instances"]
        except (requests.exceptions.JSONDecodeError, KeyError):
            if attempts > 0:
                return _chain2entity(pdbs, attempts=attempts-1)
            else:
                return []

    return [(entity["rcsb_id"], f'{entity["rcsb_id"][:4]}_{entity["rcsb_polymer_entity_instance_container_identifiers"]["entity_id"]}') for entity in results]
    

def chain2entity(pdbs, chunk_size=8000):
    values = []
    for pdb_chunk_start in range(0, len(pdbs), chunk_size):
        RealtimeLogger.info(f"Converting all chain ids to entity ids: {pdb_chunk_start}/{len(pdbs)}")
        pdb_chunk = pdbs[pdb_chunk_start:pdb_chunk_start+chunk_size]
        values += _chain2entity(pdb_chunk)
    values = pd.DataFrame(values, columns=["cath_domain", "entity_id"])   
    RealtimeLogger.info("Converted all chain ids to entity ids")
    return values

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
    RealtimeLogger.info("Downloading all PDBs names")
    pdbs = None
    if update:
        pdbs = update_pdbs(job, full_pdb_h5)

    if pdbs is None:
        experimental_query = FieldQuery("rcsb_entry_info.structure_determination_methodology", exact_match="experimental")
        #entity_query = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
        #query = CompositeQuery(experimental_query, entity_query)
        pdbs = search(experimental_query, return_type="polymer_instance")

    assert len(pdbs)>0
    
    RealtimeLogger.info("Adding all PDB names")
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group("domains")
        run_pdbs = set(pdbs)-set(store['domains'].keys())

    batch_size = 500
    batches = [pdbs[start:start+batch_size] for start in range(0,len(run_pdbs),batch_size)]
    map_job(job, create_groups, batches, full_pdb_h5)
        
        # for i, pdb in enumerate(run_pdbs):
        #     if i%500==0:
        #         RealtimeLogger.info(f"Adding {i}/{len(run_pdbs)}")
        #     store.require_group(f"domains/{pdb}")
    job.addFollowOnJobFn(finish_section, full_pdb_h5, "completed_names")

    #Map chains to entity ids
    pdb_chains = chain2entity(pdbs)
    
    job.addFollowOnJobFn(create_splits_for_levels, full_pdb_h5, pdb_chains)

def create_groups(job, pdbs, full_pdb_h5):
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        for i, pdb in enumerate(pdbs):
            if i%50==0:
                RealtimeLogger.info(f"Adding {i}/{len(pdbs)}")
            store.require_group(f"domains/{pdb}")
         
def get_custom_pdbs(job, pdbs, full_pdb_h5, update=False):
    if not isinstance(pdbs, (list, tuple)):
        pdbs = [pdbs]

    pdb_chains = set()

    #Then try old chains names (polymer_instance or asym id) (e.g. 1XYZ.A)
    possible_chains = [p for p in pdbs if "." in p]
    if len(possible_chains) > 0:    
        chain_query = FieldQuery("rcsb_polymer_entity_instance_container_identifiers.rcsb_id", is_in=pdbs)
        chains_result = set(search(chain_query, return_type="polymer_instance"))
        pdb_chains |= chains_result
    
    #Try polymer_entity first (e.g. 1XYZ_1)    
    possible_entities = [p.split("_") for p in pdbs if "_" in p]
    if len(possible_entities) > 0:
        # pdb_codes, pdb_entities = zip(*possible_entities)
        # entity_query = FieldQuery("rcsb_entry_container_identifiers.entry_id", is_in=list(set(pdb_codes))) & FieldQuery("rcsb_polymer_entity_instance_container_identifiers.entity_id", is_in=list(set(pdb_entities)))

        entity_query = None
        for pdb_code, pdb_entity in possible_entities:
            q1 = FieldQuery("rcsb_entry_container_identifiers.entry_id", exact_match=pdb_code)
            q2 = FieldQuery("rcsb_polymer_entity_instance_container_identifiers.entity_id", exact_match=pdb_entity)
            q3 = q1 & q2
            if entity_query is None:
                entity_query = q3
            else:
                entity_query |= q3
        entity_result = set(search(entity_query, return_type="polymer_entity"))
        pdb_chains |= entity_result

    #Might've used just the 4-letter PDB codes? (e.g. 1XYZ -> 1XYZ_1,1XYZ_2)
    possible_ids = [p for p in pdbs if len(p) == 4]
    if len(possible_ids) > 0:
        entry_query = FieldQuery("rcsb_entry_container_identifiers.entry_id", is_in=pdbs)
        entry_results = set(search(entry_query, return_type="polymer_instance"))
        pdb_chains |= entry_results

    if len(pdb_chains) < len(pdbs):
        RealtimeLogger.info(f"Unable to match all inputs to PDB.CHAIN. Using {len(pdb_chains)/len(pdbs)}: {pdb_chains}")

    pdb_chains = list(pdb_chains)

    assert len(pdb_chains)>0

    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group("domains")
        for pdb in pdb_chains:
            store.require_group(f"domains/{pdb}")

    pdb_chains = chain2entity(pdb_chains)
        
    job.addFollowOnJobFn(create_splits_for_levels, full_pdb_h5, pdb_chains, custom=False)

def create_representatives(job, clusters, full_pdb_h5):
    #Get representatives as entity id
    representatives = clusters["representative"].drop_duplicates()
    
    #Get representatives as chain id
    representatives = clusters[clusters.entity_id.isin(representatives)].cath_domain.drop_duplicates()

    
    missing_domains = []
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        group = store.require_group("representatives")

        for domain in representatives:
            try:
                group[domain] = store[f"domains/{domain}'"]
            except KeyError:
                missing_domains.append(domain)

        if len(missing_domains) > 0:
            store["representatives"].attrs["missing_domains"] = missing_domains
        store["representatives"].attrs["total_domains"] = len(representatives)
    
    job.addFollowOnJobFn(finish_section, full_pdb_h5, "completed_representatives")

def update_clusters(job, full_pdb_h5, pdbs, pct_id=30, custom=False, update=False):
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
    
    query = FieldQuery("rcsb_entry_info.structure_determination_methodology", exact_match="experimental")
    # entity_query = rcsb.FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
    # query = experimental_query & entity_query
    # if custom_pdbs is None:
    #     custom_pdbs = search(query, return_type="polymer_instance")
    if custom:
        query &= FieldQuery("rcsb_polymer_entity_instance_container_identifiers.rcsb_id", is_in=pdbs)


    search_results = search(query, group_by=IdentityGrouping(pct_id), return_groups=True, return_type="polymer_entity")
    clusters = pd.DataFrame([
        (pdb,group_name,group[0]) for group_name, group in search_results.items() for pdb in group], 
        columns=["entity_id", "cluster_name", "representative"])
    
    clusters = pd.merge(pdbs, clusters, how="left", on="entity_id")

    job.addFollowOnJobFn(split_dataset_at_level, full_pdb_h5, "", clusters, "cluster_name", str(pct_id),
        split_size={"train":0.8, "validation":0.1, "test":0.1})
    
    if pct_id == 30:
        job.addFollowOnJobFn(create_representatives, clusters, full_pdb_h5)

def create_splits_for_levels(job, full_pdb_h5, pdbs, custom=False):
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group(f"data_splits")

    for level in [100, 95, 90, 70, 50, 30]:#[30, 40, 50, 70, 90, 95, 100]:
        job.addChildJobFn(update_clusters, full_pdb_h5, pdbs, level, custom=custom)
    
    job.addFollowOnJobFn(finish_section, full_pdb_h5, "completed_domain_splits")

def finish_section(job, cath_full_h5, attribute):
    with h5pyd.File(cath_full_h5, mode="a", use_cache=False, retries=100) as store:
        store.attrs[attribute] = True
