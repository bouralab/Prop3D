import os, sys
import glob
from multiprocessing.pool import ThreadPool
import shutil
import itertools as it

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(*args, **kwds):
        return args

from Bio.PDB.Polypeptide import three_to_one

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import iter_unique_superfams, get_file, filter_hdf, filter_hdf_chunks
from molmimic.generate_data.job_utils import map_job, map_job_rv, map_job_rv_list
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS
from molmimic.parsers.eppic import get_interfaces, get_interface_residues

from toil.realtimeLogger import RealtimeLogger

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def collapse_binding_site(df):
    columns = ["pdb", "interfaceId", "firstChain", "secondChain"]

    #Add columns for Domain:Domain, Domain:Loop, and Loop:Domain Interactions
    columns += [c for c in df.columns if "firstCath" in c or "secondCath" in c]
    # if "secondCathDomainClosest1" in df.columns:
    #     #Domain:Loop Interactions
    #     columns += ["firstCathDomain", "firstCathCode", "secondCathDomainClosest1"]
    # elif "firstCathDomainClosest1" in df.columns:
    #     #Loop:Domain Interactions
    #     columns = ["firstCathDomainClosest1", "secondCathDomain"]
    # elif "firstCathDomain" in df.columns:
    #     #Domain:Domain Interactions
    #     columns += ["firstCathDomain", "firstCathCode", "secondCathDomain", "secondCathCode"]
    #

    binding_site = df[columns].iloc[0]
    #domain_residues["pdbResidueNumber"] = domain_residues["pdbResidueNumber"].apply(
    #    lambda x: "{}{}".format(int(x[1]), x[2].strip()))
    mol_resi = df["firstPdbResNumber"].drop_duplicates().apply(
        lambda x: "{}{}".format(int(x[1]), x[2].strip()))
    int_resi = df["secondPdbResNumber"].drop_duplicates().apply(
        lambda x: "{}{}".format(int(x[1]), x[2].strip()))
    mol_resn = df["firstResType"].loc[mol_resi.index].apply(three_to_one)
    int_resn = df["secondResType"].loc[int_resi.index].apply(three_to_one)

    binding_site = pd.concat((binding_site, pd.Series({
        "firstResi": mol_resi.str.cat(sep=','),
        "firstResn": mol_resn.str.cat(sep=','),
        "secondResi": int_resi.str.cat(sep=','),
        "secondResn": int_resn.str.cat(sep=',')
    })), axis=0)

    return binding_site

def get_binding_sites(contacts):
    """Collapse interactions to one binding site per line"""

    if "secondCathDomainClosest1" in contacts.columns:
        #Domain:Loop Interactions
        columns = ["firstCathDomain", "secondCathDomainClosest1"]
        order = [
            "pdb", "interfaceId", "reverse", "firstCathCode", "firstCathDomain",
            "firstChain", "firstResi", "firstResn", "secondCathCodeClosest1",
            "secondCathDomainClosest1", "secondCathDomainDistance1","secondChain",
            "secondResi", "secondResn"
        ]
    elif "firstCathDomainClosest1" in contacts.columns:
        #Loop:Domain Interactions
        columns = ["firstCathDomainClosest1", "secondCathDomain"]
        order = [
            "pdb", "interfaceId", "reverse", "firstCathCodeClosest1", "firstCathDomainClosest1",
            "firstCathDomainDistanceClosest1","firstChain", "firstResi", "firstResn",
            "secondCathCode", "secondCathDomain", "secondChain", "secondResi", "secondResn"
        ]
    else:
        #Domain:Domain Interactions
        columns = ["firstCathDomain", "secondCathDomain"]
        order = [
            "pdb", "interfaceId", "reverse", "firstCathCode", "firstCathDomain",
            "firstChain", "firstResi", "firstResn", "secondCathCode",
            "secondCathDomain", "secondChain", "secondResi", "secondResn"
        ]

    binding_sites = contacts.groupby(columns, as_index=False)
    binding_sites = binding_sites.apply(collapse_binding_site).reset_index(drop=True)
    binding_sites = binding_sites.assign(reverse=False)

    #Duplicate binding sites by swaping first and second chain
    binding_sites2 = binding_sites.copy()
    binding_sites2 = binding_sites2.rename(columns=
        {c:c.replace("first", "second") if c[0]=="f" else c.replace("second", "first") \
         for c in binding_sites2.columns if "first" in c or "second" in c})
    binding_sites2["reverse"] = True
    binding_sites = binding_sites.append(binding_sites2).reset_index(drop=True)

    #Reorder collapsed DDI columns
    binding_sites = binding_sites[order]

    return binding_sites

def update_contact_columns(df):
    contacts = df.drop(columns={
        'asa', 'bsa', 'bsaPercentage', 'residueNumber', 'residueType',
        'side', 'pdb'})
    contacts = contacts.rename(columns={
        "pdbResidueNumber":resSide[side]+"PdbResNumber",
        "entropyScore":resSide[side]+"EntropyScore",
        "region":resSide[side]+"Region",
        "chain":resSide[side]+"Chain",
        "cath_domain":resSide[side]+"CathDomain",
        "cathcode":resSide[side]+"CathCode",
        "uid_x":"contactUid",
        "uid_y":"firstUid",
        "uid":"secondUid"
        "pdbCode":"pdb",
        "interfaceId_x":"interfaceId"})
    #Reorder columns
    contacts = contacts[[
        'contactUid', 'pdb', 'interfaceId',
        'firstResNumber', 'firstPdbResNumber', 'firstResType', 'firstChain', 'firstCathDomain', 'firstCathCode',
        'secondResNumber', 'secondPdbResNumber', 'secondResType', 'secondChain', 'secondCathDomain', 'secondCathCode',
        'firstBurial', 'firstEntropyScore', 'firstRegion', 'firstUid',
        'secondBurial', 'secondEntropyScore', 'secondRegion', 'secondUid',
        'isClash', 'isDisulfide', 'minDistance', 'numAtoms', 'numHBonds',
    ]]

    return contacts

def get_top_cath_domains_to_loop_contacts(df_):
    df_ = df_.sort_values("distance")

    #Get info on first closest cath domain
    best1 = df_[["cath_domain", "cathcode", "distance"]].iloc[0].rename(columns={
        "cath_domain":"secondCathDomainClosest1",
        "cathcode":"secondCathCodeClosest1",
        "distance":"secondCathDomainDistance1"
    })

    #Get info on second cloests cath domain if they are equidistant
    if len(dist)>1 and df_["distance"].iloc[:2].nunique()==1:
        best2 = df_[["cath_domain", "cathcode", "distance"]].iloc[1].rename(columns={
            "cath_domain":"secondCathDomainClosest2",
            "cathcode":"secondCathCodeClosest2",
            "distance":"secondCathDomainDistance2"
        })
    else:
        best2 = pd.Series({
            "secondCathDomainClosest2":np.nan,
            "secondCathCodeClosest2":np.nan,
            "secondCathDomainDistance2":np.nan
        })

    best = pd.concat((best1, best2))

    #Select relevetn info of contact from first row
    contact_info = df_[[
        'contactUid', 'pdb', 'interfaceId',
        'firstResNumber', 'firstPdbResNumber', 'firstResType', 'firstChain', 'firstCathDomain', 'firstCathCode',
        'secondResNumber', 'secondPdbResNumber', 'secondResType', 'secondChain',
        'firstBurial', 'firstEntropyScore', 'firstRegion', 'firstUid',
        'secondBurial', 'secondEntropyScore', 'secondRegion', 'secondUid'
        'isClash', 'isDisulfide', 'minDistance', 'numAtoms', 'numHBonds'
    ]].iloc[0]

    contact_info = pd.concat((contact_info, best))

    #Reorder columns
    contact_info = contact_info[[
        'contactUid', 'pdb', 'interfaceId',
        'firstResNumber', 'firstPdbResNumber', 'firstResType', 'firstChain', 'firstCathDomain', 'firstCathCode',
        'secondResNumber', 'secondPdbResNumber', 'secondResType', 'secondChain',
        'secondCathDomainClosest1', 'secondCathCodeClosest1', 'secondCathDomainDistance1',
        'secondCathDomainClosest2', 'secondCathCodeClosest2', 'secondCathDomainDistance2',
        'firstBurial', 'firstEntropyScore', 'firstRegion', 'firstUid',
        'secondBurial', 'secondEntropyScore', 'secondRegion', 'secondUid'
        'isClash', 'isDisulfide', 'minDistance', 'numAtoms', 'numHBonds'
    ]].iloc[0]

    return contact_info

def process_pdb(job, pdbId, cathFileStoreID, work_dir=None):
    work_dir = job.fileStore.getLocalTempDir()

    if isinstance(cathFileStoreID, pd.DataFrame):
        cath = cathFileStoreID[cathFileStoreID["pdb"]==pdbId]
        should_remove = False
    else:
        cath_file, ftype = get_file(job, "cath.h5", cathFileStoreID, return_type=True)
        cath = filter_hdf(cath_file, "table", pdb=pdbId)
        should_remove = ftype == "fileStoreID"

    #Get all chain:chain interfaces for PDB from EPPIC
    interfaces = get_interfaces(pdbId.lower(), bio=True)

    resSide = ["first", "second"]

    all_domain_residues = pd.DataFrame()
    all_loop_residues = pd.DataFrame()
    all_binding_sites = pd.DataFrame()
    for interface in interfaces.itertuples():
        #Get all residues chain:chain interface from EPPIC
        interfaceResidues = get_interface_residues(pdbId.lower(), interface.interfaceId))

        #Only keep residues that are Totally buried, Surface, Rim, Core geometry.
        #or Core evolutionary
        interfaceResidues = interfaceResidues[interfaceResidues["region"]>=0]

        #Make PDB residues sortable
        interfaceResidues["pdbResidueNumber"] = interfaceResidues["pdbResidueNumber"].apply(
            natural_keys)

        interfaceResidues = interfaceResidues.assign(interfaceId=interfaceId)
        interfaceResidues["side"] = interfaceResidues["side"].astype(int)
        chainLetts = interfaceResidues[["chain1", "chain2"]].tolist()

        interfaceContacts = get_contacts(pdbId.lower(), interface.interfaceId)
        domain_domain_contacts = interfaceContacts

        all_interface_domain_residues = pd.DataFrame()
        all_interface_loop_residues = pd.DataFrame()
        for side, chain in interfaceResidues.groupby("side"):
            chainLett = chainLetts[side]
            chain_residues = chain.assign(pdb=pdbId.lower(), chain=chainLett)
            cath_chain = cath[cath["chain"]==chainLett]
            cath_chain = pd.merge(
                cath_chain[["cath_domain", "cathcode", "pdb", "chain", "srange_start", "srange_stop"]],
                chain_residues,
                how="left", on=["pdb", "chain"])

            domain_residues = cath_chain[
                (cath_chain["pdbResidueNumber"] >= cath_chain["srange_start"]) & \
                (cath_chain["pdbResidueNumber"] <= cath_chain["srange_stop"])]
            domain_residues = domain_residues.drop(columns=["srange_start", "srange_stop"])
            #domain_residues["pdbResidueNumber"] = domain_residues["pdbResidueNumber"].apply(
            #    lambda x: "{}{}".format(int(x[1]), x[2].strip()))
            all_interface_domain_residues = all_interface_domain_residues.append(domain_residues)

            loop_residues = cath_chain[
                (cath_chain["pdbResidueNumber"] >= cath_chain["srange_start"]) & \
                (cath_chain["pdbResidueNumber"] <= cath_chain["srange_stop"])]
            loop_residues = domain_residues.drop(columns=["srange_start", "srange_stop"])
            #loop_residues["pdbResidueNumber"] = loop_residues["pdbResidueNumber"].apply(
            #    lambda x: "{}{}".format(int(x[1]), x[2].strip()))
            all_loop_residues = all_loop_residues.append(domain_residues)

            #Create DDI interactome by merging with contacts
            domain_domain_contacts = pd.merge(
                domain_domain_contacts, domain_residues, how="inner",
                left_on=resSide[side]+"ResNumber", right_on="residueNumber")
            domain_domain_contacts = update_contact_columns(domain_domain_contacts)

            #Create Loop:Loop interactome by merging with contacts
            loop_loop_contacts = pd.merge(
                loop_loop_contacts, loop_residues, how="inner",
                left_on=resSide[side]+"ResNumber", right_on="residueNumber")
            loop_loop_contacts = update_contact_columns(loop_loop_contacts)

            #Complete Domain:Loop and Loop:Domain interactions
            if side == 0:
                #First side is either DDI or LLI w/o the second residue info
                loop_domain_contacts = loop_loop_contacts
                domain_loop_contacts = domain_domains_contacts
            else:
                #Second side is merged with first sie of DDI or LLI w/ residue info
                loop_domain_contacts = pd.merge(
                    loop_domain_contacts, domain_residues, how="inner",
                    left_on=resSide[side]+"ResNumber", right_on="residueNumber")
                loop_domain_contacts = update_contact_columns(loop_domain_contacts)
                domain_loop_contacts = pd.merge(
                    domain_domain_contacts, loop_residues, how="inner",
                    left_on=resSide[side]+"ResNumber", right_on="residueNumber")
                domain_loop_contacts = update_contact_columns(domain_loop_contacts)




        # #Collapse domain domain interaction to one binding site per line
        # ddi_groups = domain_domain_contacts.groupby(["firstCathDomain", "secondCathDomain"], as_index=False)
        # binding_sites = ddi_groups.apply(collapse_binding_site).reset_index(drop=True)
        # binding_sites = binding_sites.assign(reverse=False)
        #
        # #Duplicate binding sites by swaping first and second chain
        # binding_sites2 = binding_sites.copy()
        # binding_sites2 = binding_sites2.rename(columns=
        #     {c:c.replace("first", "second") if c[0]=="f" else c.replace("second", "first") \
        #      for c in binding_sites2.columns if "first" in c or "second" in c})
        # binding_sites2["reverse"] = True
        # binding_sites = binding_sites.append(binding_sites2).reset_index(drop=True)

        #Reorder collapsed DDI columns
        # binding_sites = binding_sites[[
        #     "pdb", "interfaceId", "reverse", "firstCathCode", "firstCathDomain",
        #     "firstChain", "firstResi", "firstResn", "secondCathCode",
        #     "secondCathDomain", "secondChain", "secondResi", "secondResn"
        # ]]

        #Get Domain:Loop interactions closest domains
        domain_loop_contacts_to_cath = pd.merge(
            domain_loop_contacts,
            cath_chain[["pdb", "chain", "cath_domain", "cathcode", "srange_start", "srange_stop"]],
            left_on=["pdb", "secondChain"],
            right_on=["pdb", "chain"])
        domain_loop_contacts_to_cath = domain_loop_contacts_to_cath[
            (domain_loop_contacts_to_cath["secondPdbResNumber"]>domain_loop_contacts_to_cath["srange_stop"]) &\
            (domain_loop_contacts_to_cath["secondPdbResNumber"]<domain_loop_contacts_to_cath["srange_start"])]

        dist_to_start = domain_loop_contacts_to_cath["srange_start"]-domain_loop_contacts_to_cath["secondPdbResNumber"]
        dist_to_end = domain_loop_contacts_to_cath["secondPdbResNumber"]-domain_loop_contacts_to_cath["srange_stop"]
        distance = pd.concat((dist_to_start, dist_to_end), axis=1).abs()
        domain_loop_contacts_to_cath = domain_loop_contacts_to_cath.assign(
            distance=distance.min(), position=distance.idxmin())

        top_cath_domain_loops_contacts = domain_loop_contacts_to_cath.groupby(
            ["firstCathDomain", "firstResNumber", "secondResNumber"]).apply(
                get_top_cath_domains_to_loop_contacts)

        domain_domain_binding_sites = get_binding_sites(domain_domain_contacts)
        loop_loop_binding_sites = get_binding_sites(loop_loop_contacts)
        domain_loop_binding_sites = get_binding_sites(domain_loop_contacts)
        loop_domain_binding_sites = get_binding_sites(domain_loop_contacts)

        domain_domain_contacts
        loop_loop_contacts
        loop_domain_contacts
        domain_loop_contacts

    if should_remove:
        try:
            os.remove(cath_file)
        except OSError:
            pass

def process_pdb_group(job, pdb_group, cathFileStoreID, further_parallize=False):
    work_dir = job.fileStore.getLocalTempDir()

    if further_parallelize:
        map_job(job, process_pdb, pdb_group, cathFileStoreID)
    else:
        cath_file = get_file(job, "cath.h5", cathFileStoreID)
        cath = filter_hdf(cath_file, "table", query="pdb in pdbId")
        for pdbId in pdb_group:
            try:
                process_pdb(job, pdbId, cath, work_dir=work_dir)
            except (SystemExit, KeyboardInterrupt):
                raise
            except Exception as e:
                RealtimeLogger.info("Failed getting interactome for {}".format(pdbId))

        try:
            os.remove(cath_file)
        except OSError:
            pass


def start_toil(job, further_parallelize=False):
    work_dir = job.fileStore.getLocalTempDir()

    in_store = IOStore.get("aws:us-east-1:molmimic-cath")
    sfam_file = os.path.join(work_dir, "cath.h5")
    in_store.read_input_file("cath-domain-description-file-small.h5", sfam_file)

    cathFileStoreID = job.fileStore.writeGlobalFile(sfam_file)

    pdb = filter_hdf_chunks(sfam_file, "table", columns=["pdb"],
       drop_duplicates=True).sort_values("pdb")
    pdb = pdb.assign(group=pdb["pdb"].str[1:3])
    pdb_groups = pdb.groupby("group")["pdb"].apply(list)
    map_job(job, run_cath_hierarchy, classes, cathFileStoreID, update_features)

    try:
        os.remove(sfam_file)
    except (FileNotFoundError, OSError):
        pass

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    job = Job.wrapJobFn(start_toil)
    with Toil(options) as workflow:
        workflow.start(job)

    #     domain_residues = cath_chain[
    #         (cath_chain["pdbResidueNumber"] >= cath_chain["srange_start"]) &\
    #         (cath_chain["pdbResidueNumber"] <= cath_chain["srange_stop"])]
    #     domain_residues = domain_residues.drop(columns=["srange_start", "srange_stop"])
    #     domain_residues["pdbResidueNumber"] = domain_residues1["pdbResidueNumber"].apply(
    #         lambda x: "{}{}".format(int(x[1]), x[2].strip()))
    #     all_domain_residues = pd.concat((all_domain_residues, domain_residues), axis=0)
    #
    #     loop_residues = cath_chain[
    #         (cath_chain["pdbResidueNumber"] < cath_chain["srange_start"]) &\
    #         (cath_chain["pdbResidueNumber"] > cath_chain["srange_stop"])]
    #     loop_residues = domain_residues.drop(columns=["srange_start", "srange_stop"])
    #     loop_residues["pdbResidueNumber"] = domain_residues1["pdbResidueNumber"].apply(
    #         lambda x: "{}{}".format(int(x[1]), x[2].strip()))
    #     all_loop_residues = pd.concat((all_loop_residues, loop_residues), axis=0)
    #
    #
    # contacts = pd.DataFrame(get_contacts(pdbId.lower(), interface.interfaceId))
    # del contacts["interfaceId"]
    #
    # all_domain_residues = pd.DataFrame()
    # all_loop_residues = pd.DataFrame()
    #
    # chain_interactome = []
    # domain_interactome = []
    #
    # for interface in interfaces:
    #     interfaceId = interface["interfaceId"]
    #     chain1 = interface["chain1"]
    #     chain2 = interface["chain1"]
    #     chains = [chain1, chain2]
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #     for score in interface["interfaceScores"]["interfaceScore"]:
    #         if score["name"] == "eppic":
    #             interface_type = score["callName"]
    #             break
    #     else:
    #         continue
    #         interface_type = "xtal"
    #
    #     if interface_type == "xtal":
    #         continue
    #
    #     interfaceResidues = get_interface_residues(pdbId.lower(), interfaceId)
    #
    #     residues = [{"region":[], "resi":[], "resn":[]}]*2
    #
    #     for interfaceResidue in interfaceResidues:
    #         if interfaceResidue["region"] == -1:
    #             continue
    #         side = int(interfaceResidue["side"])
    #         residues[side]["region"].append(str(interfaceResidue["region"]))
    #         residues[side]["resi"].append(interfaceResidue["pdbResidueNumber"])
    #         residues[side]["resn"].append(three_to_one[interfaceResidue["residueType"]])
    #
    #     chain_dfs = []
    #     for chain_num in range(2):
    #         chainLett = chains[chain_num]
    #         chain_residues = residues[chain_num]
    #         chain = pd.Series({
    #             "pdb":pdbId,
    #             "chain":chainLett,
    #             "resi":chain_residues["resi"],
    #             "resn":chain_residues["resn"],
    #             "resr":chain_residues["region"]
    #         })
    #         chain_dfs.append(chain)
    #
    #         #Map Binding sites to CATH
    #         cath_chain = cath[cath["chain"]==chain]
    #         chain = pd.merge(cath_chain, chain[["pdb", "chain", "resi"]],
    #             how="inner", on=["pdb", "chain"])
    #         chain = pd.concat((chain, pd.DataFrame(chain.resi.values.tolist())),
    #             axis=1)
    #
    #         #Make each residue one row
    #         chain_residues = pd.melt(sfam_cath_mol, id_vars=[
    #             "mol_sdi_id", "pdb", "chain", "srange_start", "srange_start"],
    #             val_name="residue")
    #
    #         #Get resiudes that fall inside CATH domains
    #         domain_residues = chain_residues[
    #             (chain_residues["srange_start"] >= chain_residues["residue"]) &\
    #             (chain_residues["residue"] <= chain_residues["srange_stop"])]
    #
    #         #Get resiudes that fall outputse CATH domains
    #         loop_residues = chain_residues[~chain_residues.index.isin(
    #             domain_residues.index)]
    #
    #         loop_residues = pd.merge(loop_residues, first_last_domains, on=["pdb", "chain"]))
    #         loop_residues = loop_residues.assign(
    #             before_first=loop_residues["residue"]<loop_residues["min_start"],
    #             after_last=loop_residues["residue"]>loop_residues["max_stop"])
    #
    #         loop_residues.pivot_table(
    #             index=["mol_sdi_id", "pdb", "chain", "srange_start", "srange_start"],
    #             values=["residue", "before_first", "after_last"]).reset_index()
