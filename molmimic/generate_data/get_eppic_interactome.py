import os, sys
import glob
import logging
import shutil
import json
import itertools as it
from multiprocessing.pool import ThreadPool
from itertools import groupby

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(*args, **kwds):
        return args

from Bio.PDB.Polypeptide import three_to_one

from molmimic.generate_data.iostore import IOStore
from molmimic.generate_data.util import iter_unique_superfams, get_file, \
    filter_hdf, filter_hdf_chunks, natural_keys
from molmimic.generate_data.job_utils import map_job, map_job_rv, map_job_rv_list
from molmimic.generate_data.map_residues import decode_residues, InvalidSIFTS
#from molmimic.generate_data.parse_cath import run_cath_hierarchy
from molmimic.parsers.eppic import EPPICApi

from toil.realtimeLogger import RealtimeLogger

logging.getLogger('boto3').setLevel(logging.WARNING)
logging.getLogger('botocore').setLevel(logging.WARNING)
logging.getLogger('s3transfer').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)
#
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)

# Helper fubtions for pandas apply

def collapse_binding_site(df):
    columns = ["pdb", "interfaceId", "firstChain", "secondChain"]

    #Add columns for Domain:Domain, Domain:Loop, and Loop:Domain Interactions
    columns += [c for c in df.columns if "firstCath" in c or "secondCath" in c]

    binding_site = df[columns].iloc[0]
    mol_resi = df["firstPdbResNumber"].drop_duplicates().astype(str)
    int_resi = df["secondPdbResNumber"].drop_duplicates().astype(str)

    def conv_resi_codes(r):
        try:
            return three_to_one(4)
        except:
            return r
    mol_resn = df["firstResType"].loc[mol_resi.index].apply(conv_resi_codes).astype(str)
    int_resn = df["secondResType"].loc[int_resi.index].apply(conv_resi_codes).astype(str)

    binding_site = pd.concat((binding_site, pd.Series({
        "firstResi": mol_resi.str.cat(sep=','),
        "firstResn": mol_resn.str.cat(sep=','),
        "secondResi": int_resi.str.cat(sep=','),
        "secondResn": int_resn.str.cat(sep=',')
    })), axis=0)

    return binding_site

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

class OutputWriter(object):
    def __init__(self, work_dir, store):
        self.work_dir = work_dir
        self.store = store

    def write(self, df, output_type):
        self.write_pdb(df, output_type)
        self.write_cath(df, output_type)

    def write_pdb(self, df, output_type):
        pdb, interfaceId = list(df.iloc[0][["pdb", "interfaceId"]])
        key = "pdb/{}/{}_{}.h5".format(pdb, interfaceId, output_type)
        file = os.path.join(self.work_dir, key.replace("/", "_"))
        df.to_hdf(file, "table", format="table", table=True, complevel=9,
            complib="bzip2", min_itemsize=1024)
        self.store.write_output_file(file, key)

        try:
            os.remove(file)
        except OSError:
            pass

    def write_cath(self, df, output_type):
        pdb, interfaceId = list(df.iloc[0][["pdb", "interfaceId"]])
        for (cathcode, cath_domain), cath_df in df.groupby(
          ["firstCathCode", "firstCathDomain"], as_index=False):
            key = "cath/{}/{}_{}_{}.h5".format(cathcode.replace(".", "/"), cath_domain,
                interfaceId, output_type)
            file = os.path.join(self.work_dir, key.replace("/", "_"))
            cath_df.to_hdf(file, "table", format="table", table=True, complevel=9,
                complib="bzip2", min_itemsize=1024)
            self.store.write_output_file(file, key)
            try:
                os.remove(file)
            except OSError:
                pass

class Status(dict):
    TEMPLATE = {
        "ddi": None,
        "ddi_residues": None,
        "lli": None,
        "lli_residues": None,
        "dli": None,
        "dli_residues": None
    }

    def __init__(self, pdbId, store, manual_status=False):
        self.pdbId = pdbId
        self.store = store
        self.exists = False
        self.work_dir = work_dir if work_dir is not None else os.getcwd()

        self.key = "pdb/{}/status.json".format(pdbId)
        self.path = os.path.join(self.work_dir, "{}_status.json".format(pdbId))
        self["numInterfaces"] = 0

        if store.exists(key):
            self.read()
        elif manual_status:
            self.manual_status()

    def __get__(self, key):
        try:
            return super(Status, self).__get__(key)
        except KeyError:
            self[key] = status.TEMPLATE.copy()
            return super(Status, self).__get__(key)

    def write(self):
        with open(self.path, "w") as fh:
            json.dump(self, fh)

        self.store.write_output_file(self.path, self.key)

        try:
            os.remove(self.path)
        except:
            pass

    def read(self):
        self.store.read_input_file(self.key, self.file)
        with open(self.file) as fh:
            status = json.load(fh)

        self.update(status)
        self.exists = True

    def manual_status(self):
        pdb_key = "pdb/{}".format(self.pdbId)
        for key in self.store.list_input_directory(pdb_key):
            if "status" in key: continue

            fname = os.path.splitext(key[len(pdb_key)+1:])[0]
            intId, output_type = fname.split("_", 1)

            fname = os.path.join(work_dir, fname)
            self.store.read_input_file(key, fname)
            pd.read_hdf(fname, "table")
            df = pd.HDFStore(fname)
            self[intId][output_type] = df.get_storer('table').nrows
            df.close()

            try:
                os.remove(fname)
            except OSError:
                pass

    def finished_interfaces(self):
        finished_intIds = []
        for intId, int_status in self.items():
            if intId == "numInterfaces":
                pass
            elif intId == "error":
                #Interface failed previous version due to not being in API
                skip_intIds.append(intId)
            elif isinstance(int_status, dict) and all(isinstance(v, int) for v in \
              int_status.values()):
                #Interface is done
                skip_intIds.append(intId)
            else:
                #Interface not calculated, must recalculate
                continue

        return finished_intIds

    def is_complete(self, intIds=None):
        skip_intIds = len(intIds) if intIds is None else len(self.finished_interfaces())
        if "numInterfaces" in self and self["numInterfaces"]==len(skip_intIds):
            if status["numInterfaces"] > 0 or (status["numInterfaces"] == 0 and \
              "error" in "status"):
                #PDB is Complete
                if not self.exists:
                    self.write()
            return True
        return False

class EPPICInteractome(object):
    resSide = ["first", "second"]

    def __init__(self, pdbId, cathFileStoreID, store, manual_status=True, work_dir=None):
        self.cathFileStoreID = cathFileStoreID
        self.manual_status = manual_status
        self.work_dir = work_dir if work_dir is not None else os.getcwd()

        self.status = Status(pdbId, store, work_dir=work_dir)
        self.writer = OutputWriter(work_dir, store)

        self.skip_intIds = self.status.finished_interfaces()

        #Start new connection EPPIC api
        eppic_store = IOStore.get("aws:us-east-1:molmimic-eppic-service")
        pdbe_store = IOStore.get("aws:us-east-1:molmimic-pdbe-service")
        self.eppic = EPPICApi(pdbId.lower(), eppic_store, pdbe_store,
            use_representative_chains=False, work_dir=work_dir)

    def run(self):
        if self.status.is_complete(self.skip_intIds):
            return

        #Get all chain:chain interfaces for PDB from EPPIC
        interfaces = self.eppic.get_interfaces(bio=True)

        if interfaces is None:
            #save empty status file
            self.status["error"] = "interfaces is None"
            self.status.write()
            return

        RealtimeLogger.info("IFACES {}".format(interfaces["interfaceId"]))
        if len(self.skip_intIds)>0:
            interfaces = interfaces[~interfaces.isin(self.skip_intIds)]

        if interfaces.empty:
            #save empty status file
            self.status["error"] = "interfaces is empty"
            self.status.write()
            return

        #Convert CATH start/strop ranges to numeric indices
        residue_info = self.eppic.get_residue_info()

        if residue_info is None:
            #save empty status file
            self.status["error"] = "residue_info is None"
            self.status.write()
            return

        residue_info = residue_info.drop(residue_info.columns.difference(
            ["chain", "residueNumber", "pdbResidueNumber"]), 1)
        residue_info["pdbResidueNumber"] = residue_info["pdbResidueNumber"].astype(str)

        #Get CATH ranges
        self.cath = self.get_cath_ranges(residue_info)
        del residue_info

        for interface in interfaces.itertuples():
            self.process_interface(interface)

        self.status.write()

    def process_interface(self, interface):
        intId = str(interface.interfaceId)
        RealtimeLogger.info("intID {}".format(intId))

        #Get all residues chain:chain interface from EPPIC
        interfaceResidues = self.eppic.get_interface_residues(interface.interfaceId)
        if interfaceResidues is None or interfaceResidues.empty:
            #save empty status file
            self.status[intID]["error"] = "interfaceResidues is {}".format(None if \
                interfaceResidues is None else "empty")
            return

        #Only keep residues that are "Totally buried", "Surface", "Rim",
        #"Core geometry", or "Core evolutionary"
        interfaceResidues = interfaceResidues[interfaceResidues["region"]>=0]

        interfaceResidues = interfaceResidues.assign(interfaceId=interface.interfaceId)
        interfaceResidues["side"] = interfaceResidues["side"].astype(int)
        chainLetts = (interface.chain1, interface.chain2)

        interfaceContacts = self.eppic.get_contacts(interface.interfaceId)
        if interfaceContacts is None or interfaceContacts.empty:
            #save empty status file
            self.status[intID]["error"] = "interfaceContacts is {}".format(None if \
                interfaceContacts is None else "empty")
            return

        self.domain_domain_contacts = self.loop_loop_contacts = interfaceContacts

        for side, chain in interfaceResidues.groupby("side"):
            self.process_interface_side(intId, side, chain, chainLetts[side])

        get_domain_loop_interactions(intId)
        get_binding_sites(intId, domain_domain_contacts, "ddi")
        get_binding_sites(intId, loop_loop_contacts, "lli")

        self.status["numInterfaces"] += 1

    def process_interface_side(self, intId, side, chain, chainLett):
        chain_residues = chain.assign(pdb=self.pdbId.lower(), chain=chainLett)
        self.cath_chain = self.cath[self.cath["chain"]==chainLett]
        self.cath_chain = pd.merge(
            self.cath_chain[["cath_domain", "cathcode", "pdb", "chain",
                "cathStartResidueNumber", "cathStopResidueNumber"]],
            chain_residues,
            how="left", on=["pdb", "chain"])

        domain_residues = self.cath_chain[
            (self.cath_chain["residueNumber"] >= self.cath_chain["cathStartResidueNumber"]) & \
            (self.cath_chain["residueNumber"] <= self.cath_chain["cathStopResidueNumber"])]
        domain_residues = domain_residues.drop(columns=
            ["cathStartResidueNumber", "cathStopResidueNumber"])

        loop_residues = self.cath_chain[
            (self.cath_chain["residueNumber"] < self.cath_chain["cathStartResidueNumber"]) & \
            (self.cath_chain["residueNumber"] > self.cath_chain["cathStopResidueNumber"])]
        loop_residues = loop_residues.drop(columns=
            ["cathStartResidueNumber", "cathStopResidueNumber"])

        #Create DDI interactome by merging with contacts
        self.domain_domain_contacts = pd.merge(
            self.domain_domain_contacts, domain_residues, how="inner",
            left_on=[self.resSide[side]+"ResNumber", "pdbCode" if "pdbCode" in \
                self.domain_domain_contacts.columns else "pdb", "interfaceId"],
            right_on=["residueNumber", "pdb", "interfaceId"])
        self.domain_domain_contacts = self._update_contact_columns(self.domain_domain_contacts, side)

        #Create Loop:Loop interactome by merging with contacts
        self.loop_loop_contacts = pd.merge(
            self.loop_loop_contacts, loop_residues, how="inner",
            left_on=[self.resSide[side]+"ResNumber", "pdbCode" if "pdbCode" in \
                loop_loop_contacts.columns else "pdb", "interfaceId"],
            right_on=["residueNumber", "pdb", "interfaceId"])
        self.loop_loop_contacts = self._update_contact_columns(self.loop_loop_contacts, side)

        #Complete Domain:Loop and Loop:Domain interactions
        if side == 0:
            #First side is either DDI or LLI w/o the second residue info
            pass
        else:
            #Second side is merged with first sie of DDI or LLI w/ residue info
            self.loop_domain_contacts = pd.merge(
                self.loop_loop_contacts, domain_residues, how="inner",
                left_on=[self.resSide[side]+"ResNumber", self.resSide[side]+"Chain", "pdb", "interfaceId"],
                right_on=["residueNumber", "chain", "pdb", "interfaceId"])
            self.loop_domain_contacts = self.loop_domain_contacts.drop(columns=[
                'secondChain', 'secondCathDomain', 'secondCathCode',
                'secondUid', 'secondRegion', 'secondPdbResNumber',
                'secondEntropyScore'])
            self.loop_domain_contacts = self._update_contact_columns(self.loop_domain_contacts, side)
            self.domain_loop_contacts = pd.merge(
                self.domain_domain_contacts, loop_residues, how="inner",
                left_on=[self.resSide[side]+"ResNumber", self.resSide[side]+"Chain", "pdb", "interfaceId"],
                right_on=["residueNumber", "chain", "pdb", "interfaceId"])
            self.domain_loop_contacts = self.domain_loop_contacts.drop(columns=[
                'secondChain', 'secondCathDomain', 'secondCathCode',
                'secondUid', 'secondRegion', 'secondPdbResNumber',
                'secondEntropyScore'])
            self.domain_loop_contacts = self._update_contact_columns(self.domain_loop_contacts, side)

    def get_domain_loop_interactions(self, intId):
        #Get Domain:Loop interactions closest domains
        domain_loop_contacts_to_cath = pd.merge(
            self.domain_loop_contacts,
            self.cath_chain[["pdb", "chain", "cath_domain", "cathcode", "cathStartResidueNumber", "cathStopResidueNumber"]],
            left_on=["pdb", "secondChain"],
            right_on=["pdb", "chain"])
        domain_loop_contacts_to_cath = domain_loop_contacts_to_cath[
            (domain_loop_contacts_to_cath["secondResNumber"]>domain_loop_contacts_to_cath["cathStopResidueNumber"]) &\
            (domain_loop_contacts_to_cath["secondResNumber"]<domain_loop_contacts_to_cath["cathStartResidueNumber"])]

        #Get distance from 2nd PDB Res to start and ends of closest CATH domains
        dist_to_start = domain_loop_contacts_to_cath["cathStartResidueNumber"]-domain_loop_contacts_to_cath["secondResNumber"]
        dist_to_end = domain_loop_contacts_to_cath["secondResNumber"]-domain_loop_contacts_to_cath["cathStopResidueNumber"]

        if len(dist_to_start) > 0 and len(dist_to_end) > 0:
            distance = pd.concat((dist_to_start, dist_to_end), axis=1).abs()
            domain_loop_contacts_to_cath = domain_loop_contacts_to_cath.assign(
                distance=distance.min(), position=distance.idxmin())

            top_cath_domain_loops_contacts = domain_loop_contacts_to_cath.groupby(
                ["firstCathDomain", "firstResNumber", "secondResNumber"]).apply(
                    get_top_cath_domains_to_loop_contacts)

            #Save DLI binding sites for PDB file
            self.domain_loop_binding_sites = get_binding_sites(self.domain_loop_contacts)
            self.writer.write(self.domain_loop_binding_sites, "dli")
            self.status[intId]["dli"] = len(self.domain_loop_binding_sites)

            #Save DLI contacts for PDB file with info on each residue
            self.writer.write(self.domain_loop_contacts, "dli_residues")
            self.status[intId]["dli_residues"] = len(self.domain_loop_contacts)

        else:
            self.status[intId]["dli"] = 0
            self.status[intId]["dli_residues"] = 0

    def get_binding_sites(self, intId, contacts, binding_site_type):
        if not contacts.empty:
            #Save binding sites for PDB file (collapsed residues for one site)
            binding_sites = self._get_binding_sites(contacts)
            self.writer.write(binding_sites, binding_site_type)
            self.status[intId][binding_site_type] = len(binding_sites)

            #Save DDI contacts for PDB file with info on each residue
            contact_type = "{}_residues".format(binding_site_type)
            self.writer.write(contacts, contact_type)
            self.status[intId][contact_type] = len(contacts)

        else:
            contact_type = "{}_residues".format(binding_site_type)
            self.status[intId][binding_site_type] = 0
            self.status[intId][contact_type] = 0

    def get_cath_ranges(self, residue_info):
        if isinstance(self.cathFileStoreID, pd.DataFrame):
            cath = self.cathFileStoreID[self.cathFileStoreID["pdb"]==self.pdbId]
        else:
            cath_file = fileStore.readGlobalFile(self.cathFileStoreID)
            cath = filter_hdf(cath_file, "table", pdb=pdbId)

        cath["srange_start"] = cath["srange_start"].astype(str)
        cath["srange_stop"] = cath["srange_stop"].astype(str)

        for i, cath_side in enumerate(("start", "stop")):
            cath_col = "srange_{}".format(cath_side)
            cath = pd.merge(
                cath,
                residue_info.rename(columns={"pdbResidueNumber":cath_col}),
                on=["chain", cath_col],
                how="inner")
            cath = cath.rename(columns={"residueNumber":"cath{}ResidueNumber".format(
                cath_side.title())})

        return cath

    def _update_contact_columns(self, df, side):
        if "pdb_x" in df.columns:
            df = df.rename(columns={"pdb_x":"pdb"}).drop(columns=["pdb_y"])

        contacts = df.drop(columns={
            'asa', 'bsa', 'bsaPercentage', 'residueNumber', 'residueType', 'side'})

        if "pdbCode" in contacts.columns and "pdb" in contacts.columns:
            contacts = contacts.drop(columns=["pdb"])

        contacts = contacts.rename(columns={
            "pdbResidueNumber":self.resSide[side]+"PdbResNumber",
            "entropyScore":self.resSide[side]+"EntropyScore",
            "region":self.resSide[side]+"Region",
            "chain":self.resSide[side]+"Chain",
            "cath_domain":self.resSide[side]+"CathDomain",
            "cathcode":self.resSide[side]+"CathCode",
            "uid_x":"contactUid",
            "uid_y":"firstUid",
            "uid":"secondUid",
            "pdbCode":"pdb",
            "interfaceId_x":"interfaceId"})

        #Reorder columns
        cols = [
            'contactUid', 'pdb', 'interfaceId',
            'firstResNumber', 'firstPdbResNumber', 'firstResType', 'firstChain',
            'firstCathDomain', 'firstCathCode']
        if "secondResNumber" in contacts.columns and "secondCathCode" in contacts.columns:
            cols += [
                'secondResNumber', 'secondPdbResNumber', 'secondResType', 'secondChain',
                'secondCathDomain', 'secondCathCode',
            ]
        else:
            cols += ['secondResNumber', 'secondResType', 'secondBurial']

        cols += ['firstBurial', 'firstEntropyScore', 'firstRegion', 'firstUid']

        if "secondResNumber" in contacts.columns and "secondCathCode" in contacts.columns:
            cols += ['secondEntropyScore', 'secondRegion', 'secondUid']

        cols += ['isClash', 'isDisulfide', 'minDistance', 'numAtoms', 'numHBonds']

        contacts = contacts[cols]

        return contacts

    def _get_binding_sites(contacts):
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

def process_pdb(job, pdbId, cathFileStoreID, manual_status=True, work_dir=None):
    work_dir = work_dir if work_dir is not None else job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-eppic-interfaces")
    interactome = EPPICInteractome(pdbId, cathFileStoreID, store,
        manual_status=manual_status, work_dir=work_dir)
    interactome.run()

def process_pdb_group(job, pdb_group, cathFileStoreID, further_parallelize=False):
    work_dir = job.fileStore.getLocalTempDir()

    if further_parallelize:
        map_job(job, process_pdb, pdb_group, cathFileStoreID)
    else:
        cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
        cath = filter_hdf(cath_file, "table", columns=["cath_domain", "cathcode",
            "pdb", "chain", "srange_start", "srange_stop"], drop_duplicates=True)
        cath = cath[cath["pdb"].isin(pdb_group)]
        for pdbId in pdb_group:
            try:
                process_pdb(job, pdbId, cath, work_dir=work_dir)
            except (SystemExit, KeyboardInterrupt):
                raise
            except Exception as e:
                import traceback
                RealtimeLogger.info("Failed getting interactome for {}: {} - {}".format(
                    pdbId, e.__class__.__name__, e))
                RealtimeLogger.info(traceback.format_exc())

def merge_cath(job, cathFileStoreID, further_parallelize=False):
    work_dir = job.fileStore.getLocalTempDir()
    store = IOStore.get("aws:us-east-1:molmimic-eppic-interfaces")
    store.download_input_directory("cath/{}".format(cathcode.replace(".", "/")))

    for file_type in ("ddi", "lli", "dli"):
        files = list(glob.glob("*_{}.h5".format(file_type)))
        binding_sites = dd.read_hdf(files, "table")
        binding_sites = binding_sites.repartition(nparts=20)
        binding_sites.to_hdf("{}.h5".format(file_type), "table", format="table",
            table=True, complevel=9, complib="bzip2", min_itemsize=1024)
        for f in files:
            try:
                os.remove(f)
            except:
                pass

def start_toil(job, cathFileStoreID, check=False, further_parallelize=True):
    work_dir = job.fileStore.getLocalTempDir()
    cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)

    RealtimeLogger.info("Filter CATH")
    pdb = pd.read_hdf(cath_file, "table", columns=["pdb", "cath_domain"]).drop_duplicates()

    if check:
        store = IOStore.get("aws:us-east-1:molmimic-eppic-interfaces")
        keys = list(store.list_input_directory())

        done_pdbs = []
        for pdbId, files in groupby(store.list_input_directory(), lambda k: k.split("/")[1]):
            files = [os.path.splitext("".join(k.split("/")[2:]))[0] for k in files]
            if len(files) > 1 and "status" in files:
                done_pdbs.append(pdbId)

        pdb = pdb[~pdb["pdb"].isin(done_pdbs)]

    RealtimeLogger.info("Filtered CATH {}".format(len(pdb)))
    pdb = pdb.assign(group=pdb["pdb"].str[:3])
    pdb_groups = pdb.groupby("group")["pdb"].apply(list)
    map_job(job, process_pdb_group, pdb_groups, cathFileStoreID, further_parallelize)

    #map_job(job, process_pdb, pdb["pdb"], cathFileStoreID)

    #job.addFollowOnJobFn(run_cath_hierarchy, merge_cath, None, cathFileStoreID)

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--check", default=False, action="store_true")
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    if not os.path.isfile("cath.h5"):
        in_store = IOStore.get("aws:us-east-1:molmimic-cath")
        in_store.read_input_file("cath-domain-description-file-small.h5", "cath.h5")

    with Toil(options) as workflow:
        cathFileURL = 'file://' + os.path.abspath("cath.h5")
        cathFileID = workflow.importFile(cathFileURL)
        workflow.start(Job.wrapJobFn(start_toil, cathFileID, options.check))

    # status_key = "pdb/{}/status.json".format(pdbId)
    # status_file = os.path.join(work_dir, "{}_status.json".format(pdbId))
    #
    # status = Status(status_file, status_key, store)
    # writer = OutputWriter(work_dir, store)
    #
    # skip_intIds = status.finished_interfaces()
    # if status.is_complete(skip_intIds):
    #     return
    #
    # #Start new connection EPPIC api
    # eppic_store = IOStore.get("aws:us-east-1:molmimic-eppic-service")
    # pdbe_store = IOStore.get("aws:us-east-1:molmimic-pdbe-service")
    # eppic = EPPICApi(pdbId.lower(), eppic_store, pdbe_store,
    #     use_representative_chains=False, work_dir=work_dir)
    #
    # #Get all chain:chain interfaces for PDB from EPPIC
    # interfaces = eppic.get_interfaces(bio=True)
    #
    # if interfaces is None:
    #     #save empty status file
    #     status["error"] = "interfaces is None"
    #     status.write()
    #     return
    #
    # RealtimeLogger.info("IFACES {}".format(interfaces["interfaceId"]))
    # if len(skip_intIds)>0:
    #     interfaces = interfaces[~interfaces.isin(skip_intIds)]
    #
    # if interfaces.empty:
    #     #save empty status file
    #     status["error"] = "interfaces is empty"
    #     status.write()
    #     return
    #
    # #Convert CATH start/strop ranges to numeric indices
    # residue_info = eppic.get_residue_info()
    #
    # if residue_info is None:
    #     #save empty status file
    #     status["error"] = "residue_info is None"
    #     status.write()
    #     return
    #
    # residue_info = residue_info.drop(residue_info.columns.difference(
    #     ["chain", "residueNumber", "pdbResidueNumber"]), 1)
    # residue_info["pdbResidueNumber"] = residue_info["pdbResidueNumber"].astype(str)
    #
    # #Get CATH ranges
    # cath = get_cath_ranges(pdbId, cathFileStoreID, residue_info)
    #
    # del residue_info
    #
    # resSide = ["first", "second"]
    #
    # for interface in interfaces.itertuples():
    #     intId = str(interface.interfaceId)
    #     RealtimeLogger.info("intID {}".format(intId))
    #     if intId not in status:
    #         status[intId] = status.TEMPLATE.copy()
    #
    #     #Get all residues chain:chain interface from EPPIC
    #     interfaceResidues = eppic.get_interface_residues(interface.interfaceId)
    #     if interfaceResidues is None or interfaceResidues.empty:
    #         #save empty status file
    #         status[intID]["error"] = "interfaceResidues is {}".format(None if \
    #             interfaceResidues is None else "empty")
    #         continue
    #
    #
    #     #Only keep residues that are "Totally buried", "Surface", "Rim",
    #     #"Core geometry", or "Core evolutionary"
    #     interfaceResidues = interfaceResidues[interfaceResidues["region"]>=0]
    #
    #     interfaceResidues = interfaceResidues.assign(interfaceId=interface.interfaceId)
    #     interfaceResidues["side"] = interfaceResidues["side"].astype(int)
    #     chainLetts = (interface.chain1, interface.chain2)
    #
    #     interfaceContacts = eppic.get_contacts(interface.interfaceId)
    #     if interfaceContacts is None or interfaceContacts.empty:
    #         #save empty status file
    #         status[intID]["error"] = "interfaceContacts is {}".format(None if \
    #             interfaceContacts is None else "empty")
    #         continue
    #
    #     domain_domain_contacts = loop_loop_contacts = interfaceContacts
    #     domain_domain_contacts.eppic.type = "ddi"
    #
    #     for side, chain in interfaceResidues.groupby("side"):
    #         chainLett = chainLetts[side]
    #         chain_residues = chain.assign(pdb=pdbId.lower(), chain=chainLett)
    #         cath_chain = cath[cath["chain"]==chainLett]
    #         cath_chain = pd.merge(
    #             cath_chain[["cath_domain", "cathcode", "pdb", "chain",
    #                 "cathStartResidueNumber", "cathStopResidueNumber"]],
    #             chain_residues,
    #             how="left", on=["pdb", "chain"])
    #
    #         domain_residues = cath_chain[
    #             (cath_chain["residueNumber"] >= cath_chain["cathStartResidueNumber"]) & \
    #             (cath_chain["residueNumber"] <= cath_chain["cathStopResidueNumber"])]
    #         domain_residues = domain_residues.drop(columns=
    #             ["cathStartResidueNumber", "cathStopResidueNumber"])
    #
    #         loop_residues = cath_chain[
    #             (cath_chain["residueNumber"] < cath_chain["cathStartResidueNumber"]) & \
    #             (cath_chain["residueNumber"] > cath_chain["cathStopResidueNumber"])]
    #         loop_residues = loop_residues.drop(columns=
    #             ["cathStartResidueNumber", "cathStopResidueNumber"])
    #
    #         #Create DDI interactome by merging with contacts
    #         domain_domain_contacts = pd.merge(
    #             domain_domain_contacts, domain_residues, how="inner",
    #             left_on=[resSide[side]+"ResNumber", "pdbCode" if "pdbCode" in \
    #                 domain_domain_contacts.columns else "pdb", "interfaceId"],
    #             right_on=["residueNumber", "pdb", "interfaceId"])
    #         domain_domain_contacts = update_contact_columns(domain_domain_contacts, side)
    #
    #         #Create Loop:Loop interactome by merging with contacts
    #         loop_loop_contacts = pd.merge(
    #             loop_loop_contacts, loop_residues, how="inner",
    #             left_on=[resSide[side]+"ResNumber", "pdbCode" if "pdbCode" in \
    #                 loop_loop_contacts.columns else "pdb", "interfaceId"],
    #             right_on=["residueNumber", "pdb", "interfaceId"])
    #         loop_loop_contacts = update_contact_columns(loop_loop_contacts, side)
    #
    #         #Complete Domain:Loop and Loop:Domain interactions
    #         if side == 0:
    #             #First side is either DDI or LLI w/o the second residue info
    #             pass
    #         else:
    #             #Second side is merged with first sie of DDI or LLI w/ residue info
    #             loop_domain_contacts = pd.merge(
    #                 loop_loop_contacts, domain_residues, how="inner",
    #                 left_on=[resSide[side]+"ResNumber", resSide[side]+"Chain", "pdb", "interfaceId"],
    #                 right_on=["residueNumber", "chain", "pdb", "interfaceId"])
    #             loop_domain_contacts = loop_domain_contacts.drop(columns=[
    #                 'secondChain', 'secondCathDomain', 'secondCathCode',
    #                 'secondUid', 'secondRegion', 'secondPdbResNumber',
    #                 'secondEntropyScore'])
    #             loop_domain_contacts = update_contact_columns(loop_domain_contacts, side)
    #             domain_loop_contacts = pd.merge(
    #                 domain_domain_contacts, loop_residues, how="inner",
    #                 left_on=[resSide[side]+"ResNumber", resSide[side]+"Chain", "pdb", "interfaceId"],
    #                 right_on=["residueNumber", "chain", "pdb", "interfaceId"])
    #             domain_loop_contacts = domain_loop_contacts.drop(columns=[
    #                 'secondChain', 'secondCathDomain', 'secondCathCode',
    #                 'secondUid', 'secondRegion', 'secondPdbResNumber',
    #                 'secondEntropyScore'])
    #             domain_loop_contacts = update_contact_columns(domain_loop_contacts, side)
    #
    #     #Get Domain:Loop interactions closest domains
    #     domain_loop_contacts_to_cath = pd.merge(
    #         domain_loop_contacts,
    #         cath_chain[["pdb", "chain", "cath_domain", "cathcode", "cathStartResidueNumber", "cathStopResidueNumber"]],
    #         left_on=["pdb", "secondChain"],
    #         right_on=["pdb", "chain"])
    #     domain_loop_contacts_to_cath = domain_loop_contacts_to_cath[
    #         (domain_loop_contacts_to_cath["secondResNumber"]>domain_loop_contacts_to_cath["cathStopResidueNumber"]) &\
    #         (domain_loop_contacts_to_cath["secondResNumber"]<domain_loop_contacts_to_cath["cathStartResidueNumber"])]
    #
    #     #Get distance from 2nd PDB Res to start and ends of closest CATH domains
    #     dist_to_start = domain_loop_contacts_to_cath["cathStartResidueNumber"]-domain_loop_contacts_to_cath["secondResNumber"]
    #     dist_to_end = domain_loop_contacts_to_cath["secondResNumber"]-domain_loop_contacts_to_cath["cathStopResidueNumber"]
    #
    #     if len(dist_to_start) > 0 and len(dist_to_end) > 0:
    #         distance = pd.concat((dist_to_start, dist_to_end), axis=1).abs()
    #         domain_loop_contacts_to_cath = domain_loop_contacts_to_cath.assign(
    #             distance=distance.min(), position=distance.idxmin())
    #
    #         top_cath_domain_loops_contacts = domain_loop_contacts_to_cath.groupby(
    #             ["firstCathDomain", "firstResNumber", "secondResNumber"]).apply(
    #                 get_top_cath_domains_to_loop_contacts)
    #
    #         #Save DLI binding sites for PDB file
    #         domain_loop_binding_sites = get_binding_sites(domain_loop_contacts)
    #         writer.write(domain_loop_binding_sites, "dli")
    #         status[intId]["dli"] = len(domain_loop_binding_sites)
    #
    #         #Save DLI contacts for PDB file with info on each residue
    #         writer.write(domain_loop_contacts, "dli_residues")
    #         status[intId]["dli_residues"] = len(domain_loop_contacts)
    #
    #         del domain_loop_contacts
    #         del domain_loop_binding_sites
    #     else:
    #         status[intId]["dli"] = 0
    #         status[intId]["dli_residues"] = 0
    #
    #     if not domain_domain_contacts.empty:
    #         #Save DDI binding sites for PDB file (collapsed residues for one site)
    #         domain_domain_binding_sites = get_binding_sites(domain_domain_contacts)
    #         writer.write(domain_domain_binding_sites, "ddi")
    #         status[intId]["ddi"] = len(domain_domain_binding_sites)
    #
    #         #Save DDI contacts for PDB file with info on each residue
    #         writer.write(domain_domain_contacts, "ddi_residues")
    #         status[intId]["ddi_residues"] = len(domain_domain_contacts)
    #
    #         del domain_domain_binding_sites
    #     else:
    #         status[intId]["ddi"] = 0
    #         status[intId]["ddi_residues"] = 0
    #
    #     if not loop_loop_contacts.empty:
    #         #Save LLI binding sites for PDB file (collapsed residues for one site)
    #         loop_loop_binding_sites = get_binding_sites(loop_loop_contacts)
    #         writer.write(loop_loop_binding_sites, "lli")
    #         status[intId]["lli"] = len(loop_loop_binding_sites)
    #
    #         #Save LLI contacts for PDB file with info on each residue
    #         writer.write(loop_loop_contacts, "lli_residues")
    #         status[intId]["lli_residues"] = len(loop_loop_contacts)
    #
    #         del loop_loop_binding_sites
    #     else:
    #         status[intId]["lli"] = 0
    #         status[intId]["lli_residues"] = 0
    #
    #     status["numInterfaces"] += 1
    #
    #     del domain_domain_contacts
    #     del loop_loop_contacts
    #
    # status.write()


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
    #             "mol_sdi_id", "pdb", "chain", "cathStartResidueNumber", "cathStartResidueNumber"],
    #             val_name="residue")
    #
    #         #Get resiudes that fall inside CATH domains
    #         domain_residues = chain_residues[
    #             (chain_residues["cathStartResidueNumber"] >= chain_residues["residue"]) &\
    #             (chain_residues["residue"] <= chain_residues["cathStopResidueNumber"])]
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
    #             index=["mol_sdi_id", "pdb", "chain", "cathStartResidueNumber", "cathStartResidueNumber"],
    #             values=["residue", "before_first", "after_last"]).reset_index()
