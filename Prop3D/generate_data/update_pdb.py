import time
import urllib.request
from typing import Union, Any
from datetime import datetime

import h5pyd
import requests
import pandas as pd
from toil.job import Job
from biotite.database.rcsb import FieldQuery, CompositeQuery, search, count, IdentityGrouping

from Prop3D.util.toil import map_job
from toil.realtimeLogger import RealtimeLogger

from Prop3D.generate_data.create_data_splits import split_dataset_at_level


def _chain2entity(pdbs: list[str], attempts: int = 1) -> list[tuple[str,str]]:
    """Convert chain name to entity id to bridge PDB naming schemes using the PDB graph query language

    Parameters
    ----------
    pdbs : list of str
        PDB.CHAIN ids you want to convert
    attempts : int
        Number of extra attempts to make if the conversion fails. Defualt is 1

    Returns
    -------
    A list of tuple the original id and the converted ids in the same order as the input
    """
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
    

def chain2entity(pdbs: list[str], chunk_size: int = 8000) -> list[tuple[str,str]]:
    """Convert chain name to entity id to bridge PDB naming schemes using the PDB graph query language. 
    If the number of PDB chian ids is greater than 8000 (or chink_size) they will be split into 
    multiple chunks to avoid errors.

    Parameters
    ----------
    pdbs : list of str
        PDB.CHAIN ids you want to convert
    chunk_size : int
        Number of PDB ids in a single chunk. Defualt is 8000 since that is what seems to work best without failing.

    Returns
    -------
    A list of tuple the original id and the converted ids in the same order as the input
    """
    values = []
    for pdb_chunk_start in range(0, len(pdbs), chunk_size):
        RealtimeLogger.info(f"Converting all chain ids to entity ids: {pdb_chunk_start}/{len(pdbs)}")
        pdb_chunk = pdbs[pdb_chunk_start:pdb_chunk_start+chunk_size]
        values += _chain2entity(pdb_chunk)
    values = pd.DataFrame(values, columns=["cath_domain", "entity_id"])   
    RealtimeLogger.info("Converted all chain ids to entity ids")
    return values

def update_pdbs(job: Job, full_pdb_h5: str, return_query: bool = False) -> Union[set[str], FieldQuery]:
    """Get new PDBs since last update

    Parameters
    ----------
    job : toil Job
        Current job
    full_pdb_h5 : str
        Path to h5 file on hsds endpoint
    return_query : bool
        Only return the query used to obtain new PDBs. Defualt false
    
    Returns
    -------
    Either a set of new PDB ids or a biotite FieldQuery
    """
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
    is_prot = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
    query = CompositeQuery(prev_query, today_query, experimental_query, is_prot)

    new_pdbs = search(query, return_type="polymer_entity")

    last_updated_revised_query = FieldQuery("rcsb_accession_info.revision_date", greater_or_equal=last_updated)
    today_revised_query = FieldQuery("rcsb_accession_info.revision_date", less_or_equal=today)
    full_query = last_updated_revised_query & today_revised_query

    if return_query:
        return full_query
    
    revised_pdbs = search(full_query, return_type="polymer_entity")

    with h5pyd.open(full_pdb_h5, use_cache=False) as f:
        for rpdb in revised_pdbs:
            del f[f"domain/{rpdb}"]

    return set(new_pdbs)&set(revised_pdbs)

def get_all_pdbs(job: Job, full_pdb_h5: str, update: bool = False, create_groups: bool = True) -> Union[None, pd.DataFrame]:
    """Add the list of All PDB_CHAIN from the enitre PDB to the h5 file. 

    Parameters
    ----------
    job : toil.Job
        Currently running job
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    update : bool
        Add in new entries since last update
    create_groups : bool
        Create groups for all PDB_CHAIN or wait to be added later

    Returns
    -------
    Either nothing or just the list of PDB IDs to run
    """
    RealtimeLogger.info("Downloading all PDBs names")
    pdbs = None
    if update:
        pdbs = update_pdbs(job, full_pdb_h5)

    if pdbs is None:
        experimental_query = FieldQuery("rcsb_entry_info.structure_determination_methodology", exact_match="experimental")
        is_prot = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
        #entity_query = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)
        #query = CompositeQuery(experimental_query, entity_query)
        pdbs = search(experimental_query&is_prot, return_type="polymer_instance")

    assert len(pdbs)>0
    
    RealtimeLogger.info("Adding all PDB names")
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group("domains")
        run_pdbs = set(pdbs)-set(store['domains'].keys())

    #Map chains to entity ids
    pdb_chains = chain2entity(pdbs)

    if not create_groups:
        return pdb_chains

    batch_size = 500
    batches = [pdbs[start:start+batch_size] for start in range(0,len(run_pdbs),batch_size)]
    map_job(job, create_groups, batches, full_pdb_h5)
        
    job.addFollowOnJobFn(finish_section, full_pdb_h5, "completed_names")
    
    job.addFollowOnJobFn(create_splits_for_levels, full_pdb_h5, pdb_chains)

def create_groups(job: Job, pdbs: list[str], full_pdb_h5: str) -> None:
    """Create subjob for adding a subset of PDB_CHAINS to the h5 hsds store

    Parameters
    ----------
    job : toil.Job
        Currently running job
    pdbs: list of str
        The list of PDB ids. Must be PDB_CHAIN.
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    """
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        for i, pdb in enumerate(pdbs):
            if i%50==0:
                RealtimeLogger.info(f"Adding {i}/{len(pdbs)}")
            store.require_group(f"domains/{pdb}")
         
def get_custom_pdbs(job: Job, pdbs: Union[str,list[str]], full_pdb_h5: str, update: bool = False, create_groups: bool = True) -> None:
    """Add the list of custom PDBs to the h5 file. 

    Parameters
    ----------
    job : toil.Job
        Currently running job
    pdbs: list of str
        The list of PDB ids. Can be PDB_CHAIN, or PDB_ENTITY which maps to PDB_CHAIN, or 
        just 4- letter PDB code, which maps to all PDB_CHAINS present in protein.
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    update : bool
        Add in new entries since last update
    create_groups : bool
        Create groups for all PDB_CHAIN or wait to be added later

    Returns
    -------
    Either nothing or just the list of PDB IDs to run
    """
    if not isinstance(pdbs, (list, tuple)):
        pdbs = [pdbs]

    pdb_chains = set()

    if update:
        update_query = update_pdbs(job, full_pdb_h5, return_query=True)

    def search_(query, **kwds):
        if update:
            query &= update_query
        
        return search(query, **kwds)
    
    is_prot = FieldQuery("rcsb_entry_info.polymer_entity_count_protein", greater=0)

    #Then try old chains names (polymer_instance or asym id) (e.g. 1XYZ.A)
    possible_chains = [p for p in pdbs if "." in p]
    if len(possible_chains) > 0:    
        chain_query = FieldQuery("rcsb_polymer_entity_instance_container_identifiers.rcsb_id", is_in=pdbs)
        chains_result = set(search_(chain_query&is_prot, return_type="polymer_instance"))
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
        entity_result = set(search_(entity_query&is_prot, return_type="polymer_entity"))
        pdb_chains |= entity_result

    #Might've used just the 4-letter PDB codes? (e.g. 1XYZ -> 1XYZ_1,1XYZ_2)
    possible_ids = [p for p in pdbs if len(p) == 4]
    if len(possible_ids) > 0:
        entry_query = FieldQuery("rcsb_entry_container_identifiers.entry_id", is_in=pdbs)
        entry_results = set(search_(entry_query&is_prot, return_type="polymer_instance"))
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

    if not create_groups:
        return pdb_chains
    
    job.addFollowOnJobFn(create_splits_for_levels, full_pdb_h5, pdb_chains, custom=False)

    return pdb_chains

def create_representatives(job: Job, clusters: pd.DataFrame, full_pdb_h5: str) -> None:
    """Add hard links from domain representatives to its equivalent in the "domains" branch hierarchy.

    Parameters
    ----------
    job : toil.Job
        Currently running job
    clusters : pd.DrataFrame
        A dataframe mapping pdb_chain to clusters with cols ["cath_id, "entity_id", "cluster_name", "representative"]
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    """
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

def update_clusters(job: Job, full_pdb_h5: str, pdbs: pd.DataFrame, pct_id: int = 30, custom: bool = False) -> None:
    """Create jobs to create data splits based on mmseqs %seq ids

    Parameters
    ----------
    job : toil.Job
        Currently running job
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    pdbs : pd.DataFrame with cols ["cath_domain", "entity_id"] 
        A data frame mapping where cath ID (Just PDB_CHAIN) to its entity_id
    pct_id : int
        MMseqs cluster level provided by the PDB. Must be 30, 40, 50, 70, 90, 95, 100. Default is 30.
    custom : bool
        Use a custom set of PDBs. Defualt false, for all PDBs
    """
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

def create_splits_for_levels(job: Job, full_pdb_h5: str, pdbs: list[str], custom: bool = False) -> None:
    """Create jobs to create data splits based on mmseqs %seq ids

    Parameters
    ----------
    job : toil.Job
        Currently running job
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    pdbs : list of str
        list of PDB.ENTITY
    custom : bool
        Use a custom set of PDBs. Defualt false, for all PDBs
    """
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.require_group(f"data_splits")

    for level in [100, 95, 90, 70, 50, 30]:#[30, 40, 50, 70, 90, 95, 100]:
        job.addChildJobFn(update_clusters, full_pdb_h5, pdbs, level, custom=custom)
    
    job.addFollowOnJobFn(finish_section, full_pdb_h5, "completed_domain_splits")

def finish_section(job: Job, full_pdb_h5: str, attribute: str) -> None:
    """Add attribute to entire dataset (root level) specifying which sections have been completed

    Parameters
    ----------
    job : toil.Job
        Currently running job
    full_pdb_h5 : str
        Path to h5 file on HSDS endpoint
    atrribute : str
        Name of complete section
    """
    with h5pyd.File(full_pdb_h5, mode="a", use_cache=False, retries=100) as store:
        store.attrs[attribute] = True
