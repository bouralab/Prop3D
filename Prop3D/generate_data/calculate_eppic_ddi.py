import os, sys
import shutil
from itertools import chain, combinations

import pandas as pd
from Bio.PDB.Polypeptide import three_to_one
from joblib import Parallel, delayed

from Prop3D.common.ProteinTables import three_to_one
from Prop3D.util import SubprocessChain, safe_remove, getcwd
from Prop3D.util.iostore import IOStore, FileS3IOStore
from Prop3D.util.cath import run_cath_hierarchy
from Prop3D.util.hdf import get_file, filter_hdf_chunks
from Prop3D.util.toil import map_job, map_job_follow_ons, loop_job, finish_group
from Prop3D.util.pdb import PDB_TOOLS, rottrans_from_symop

from Prop3D.parsers.eppic import EPPICApi
from Prop3D.parsers.MaxCluster import MaxCluster
from Prop3D.parsers.USEARCH import ClusterFast

from Prop3D.generate_data.prepare_protein import process_domain
from Prop3D.generate_data.calculate_features import calculate_features
from Prop3D.generate_data import data_stores

from toil.realtimeLogger import RealtimeLogger

import logging
logging.getLogger('boto3').setLevel(logging.WARNING)
logging.getLogger('botocore').setLevel(logging.WARNING)
logging.getLogger('s3transfer').setLevel(logging.WARNING)
logging.getLogger('urllib3').setLevel(logging.WARNING)
logging.getLogger('matplotlib').setLevel(logging.WARNING)



def download_pdb(job, cath_code, cath_domain, raw=False, work_dir=None):
    if work_dir is None:
        work_dir = os.getcwd()

    store = IOStore.get("aws:us-east-1:prepared-structures")

    prefix = "{}/{}.pdb".format(cath_code.replace(".", "/"), cath_domain)
    file = os.path.join(work_dir, os.path.basename(prefix))
    if raw:
        prefix += ".raw"

    assert in_store.exists(prefix)

    in_store.read_input_file(prefix, file)

    return file

def from_pdb(job, pdb, assembly, advanced_dock=False, work_dir=None):
    store = IOStore.get("aws:us-east-1:eppic-interfaces")

    pdbe_store = IOStore.get("aws:us-east-1:Prop3D-pdbe-service")
    eppic_store = IOStore.get("aws:us-east-1:Prop3D-eppic-service")
    eppic_api = EPPICApi(pdb, eppic_store, pdbe_store, work_dir=self.work_dir)
    sym_op = eppic_api.get_interfaces()[["interfaceId", "operator"]].set_index(interfaceId).T.to_dict("records")[0]

    interfaces_key = "pdb/{}/{}_ddi.h5".format(pdb, assembly)
    interfaces_file = "{}_{}_ddi.h5".format(pdb, assembly)
    store.read_input_file(interfaces_key, interfaces_file)

    interfaces = pd.read_hdf(interfaces_file, "table")
    interfaces = interfaces[interfaces["reverse"]==False]

    cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
    cath_domains = pd.read_hdf(cath_file, "table", columns=["chain", "cath_domain"],
        where="pdb={}".format(pdb)).drop_duplicates()

    all_domains = {}

    #Each DDI indivudally
    for i, ddi_interface in enumerate(interfaces):
        if ddi_interface["firstCathDomain"] not in all_domains:
            all_domains[ddi_interface["firstCathDomain"]] = download_pdb(d1)

        if ddi_interface["secondCathDomain"] not in all_domains:
            secondDomain = download_pdb(d2)
            all_domains[ddi_interface["secondCathDomain"]] = d2

        id = ddi_interface["interfaceId"]

        d2 = rottrans_from_symop(
            all_domains[ddi_interface["secondCathDomain"]],
            sym_ops[id])

        orig_ddi_complex_file = next(prep(
            all_domains[ddi_interface["firstCathDomain"]],
            d2,
            merge=True, work_dir=work_dir))

        face1 = ddi_interface["firstResi"].split(",", expand=True).astype(int)
        face2 = ddi_interface["secondResi"].split(",", expand=True).astype(int)

        int_id = "{}/{}/{}/{}:{}".format(pdb, assembly, id,
            ddi_interface["firstCathDomain"],
            ddi_interface["secondCathDomain"])

        original_complex = Complex(
            orig_ddi_complex_file,
            face1=face1,
            face2=face2,
            method="original",
            work_dir=work_dir,
            job=job)

        job.addChildJobFn(refine_complex, d1_file, face1, d2_file, face2,
            orig_ddi_complex_file, int_id=int_id, advanced_dock=advanced_dock,
            cores=4 if advanced_dock else 1,
            work_dir=work_dir)

    #Each chain_domains:chain_domains individually
    for id, ddi in interfaces.groupby("interfaceId"):
        total_chains = {c:cath_domains[cath_domain["chain"]==c].drop_duplicates() \
            for c in set(ddi["firstChain"])+set(ddi["secondChain"])}

        # if all(len(c["chain"].drop_duplicates()) == 1 for c in total_chains):
        #     #DDI is same as chain DDI, do not repeat
        #     continue

        int_id = "{}/{}/ddi_{}:ddi_{}".format(pdb, assembly,
            ddi["firstChain"].str.cat(sep=''),
            ddi["secondCathDomain"].str.cat(sep=''))

        domain_files = [all_domains.get(d, download_pdb(d)) for domains in \
            total_chains for d in set(domains["firstCathDomain"])+set(domains["secondCathDomain"])]

        orig_complex_file = next(prep(*domain_files, merge=True, work_dir=work_dir))
        face1 = interfaces["firstResi"].split(",", expand=True).astype(int).values.flatten()
        face2 = interfaces["secondResi"].split(",", expand=True).astype(int).values.flatten()

        job.addChildJobFn(refine_complex, None, face1, None, face2,
            orig_complex_file, int_id=int_id, advanced_dock=False,
            cores=1, work_dir=work_dir)

        for domains in total_chains:
            pass

        #Combinations
        firstDomains = interfaces.groupby('firstCathDomain')
        for num_first_domains in range(firstDomains.ngroups):
            for firstDomain in combination(firstDomains, num_first_domains):
                for num_first_domains in range(firstDomains.ngroups):
                    for firstDomain in combination(firstDomains, num_first_domains):
                        pass

    #All domains at once


    #All at once with unk domains

def extract_binding_sites(job, cath_domain, superfamily, cathFileStoreID, work_dir=None, output_dir=None):
    #Get all binding sites for cathID
    #Get all pdb files for domains with binding sites
    #Extract binding sites
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    eppic_interfaces_store = FileS3IOStore("us-east-1", "eppic-interfaces", file_path_dir=output_dir)

    command = [sys.executable, os.path.join(PDB_TOOLS, "pdb_rslice.py")]

    pdb_key = "{}/{}.pdb.raw".format(superfamily, cath_domain)
    pdb_file = os.path.join(work_dir, "{}.pdb".format(cath_domain))
    try:
        data_stores.prepared_cath_structures.read_input_file(pdb_key, pdb_file)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        #Raise warning
        return

    ddi = [ddi_f for ddi_f in eppic_interfaces_store.list_input_directory( \
        "cath/{}/{}".format(superfamily, cath_domain)) if ddi_f.endswith("_ddi.h5")]
    RealtimeLogger.info("DDI {} {}".format(ddi, superfamily))

    fasta_fname = os.path.join(work_dir, "{}.fasta".format(cath_domain))

    with open(fasta_fname, "w") as fasta_file:
        for bs_key in ddi:
            bs_file = os.path.join(work_dir, os.path.basename(bs_key))
            try:
                eppic_interfaces_store.read_input_file(bs_key, bs_file)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                raise

            binding_sites = pd.read_hdf(bs_file, "table")
            for i, binding_site in enumerate(binding_sites.itertuples()):
                extracted_bs_file = "{}_{}_{}.pdb".format(
                    os.path.splitext(bs_file)[0], i,
                    binding_site.interfaceId)
                extracted_bs_key = "binding_sites_3d/{}/{}".format(
                    superfamily, os.path.basename(extracted_bs_file))

                resi = binding_site.firstResi.split(",")
                resn = [three_to_one(r) for r in binding_site.firstResn.split(",")]

                print(">{}\n{}".format(
                    os.path.splitext(os.path.basename(extracted_bs_file))[0],
                    "".join(resn)),
                    file=fasta_file)

                if eppic_interfaces_store.exists(extracted_bs_key):
                    continue

                with open(extracted_bs_file, "w") as fh:
                    cmd = [command+resi+[pdb_file]]
                    RealtimeLogger.info("Running {}".format(cmd))
                    SubprocessChain(cmd, fh)

                eppic_interfaces_store.write_output_file(
                    extracted_bs_file, extracted_bs_key)

    eppic_interfaces_store.write_output_file(
        fasta_fname, "binding_sites_3d/{}/{}".format(
            superfamily, os.path.basename(fasta_fname)))

def cluster_binding_sites(job, cathFileStoreID, cathcode=None, seq=False, no_struct=False, work_dir=None, output_dir=None):
    if seq:
        job.addChildJobFn(cluster_binding_site_sequences, cathcode=cathcode, work_dir=work_dir,
            output_dir=output_dir)

    if no_struct:
        return

    job.addChildJobFn(cluster_binding_site_structures, cathFileStoreID, cathcode=cathcode, work_dir=work_dir,
        output_dir=output_dir)

def cluster_binding_site_sequences(job, cathcode=None, work_dir=None, output_dir=None, force=False):
    for seqID in [1.0, 0.95, 0.60, 0.35]:
        cluster_file = "binding_sites_3d/sequence_clusters/all_binding_site_clusters_uc_{}.h5".format(seqID)
        if not force and data_stores.eppic_interfaces.exists(cluster_file):
            continue
        job.addChildJobFn(cluster_binding_site_sequences_with_seqID,
            seqID, cathcode=cathcode, work_dir=work_dir, output_dir=output_dir)

def cluster_binding_site_sequences_with_seqID(job, seqID, cathcode=None, work_dir=None, output_dir=None, cores=20, memory="256G"):
    #Run max cluster on all binding site structures
    #Save clusters
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    eppic_interfaces_store = FileS3IOStore("us-east-1", "eppic-interfaces", file_path_dir=output_dir)

    all_binding_sites_key = "binding_sites_3d/all_binding_sites.fasta"
    all_binding_sites_file = os.path.join(output_dir, "eppic-interfaces", all_binding_sites_key)

    if not eppic_interfaces_store.exists(all_binding_sites_key):
        all_binding_sites = [os.path.join(output_dir, "eppic-interfaces", bs) for bs in \
            eppic_interfaces_store.list_input_directory("binding_sites_3d/"+key) if bs.endswith(".fasta")]

        with open(all_binding_sites_file, "w") as all_binding_sites_fh:
            for domain in all_binding_sites:
                with open(domain, "r") as domain_fh:
                    shutil.copyfileobj(domain_fh, all_binding_sites_fh)

        eppic_interfaces_store.write_output_file(all_binding_sites_file, all_binding_sites_key)

    cluster = ClusterFast(job=job, work_dir=work_dir, detach=True)
    cluster_files = cluster.cluster(fasta_file=all_binding_sites_file, id=seqID, threads=cores)

    def save_cluster_msa_file(cluster_msa_file):
        data_stores.eppic_interfaces.write_output_file(
            cluster_msa_file,
            "binding_sites_3d/sequence_clusters/seqID_{}/{}.fasta".format(
                seqID, os.path.basename(cluster_msa_file).split(".")[-1]
            ))
        safe_remove(cluster_msa_file)

    cluster_h5_base = "all_binding_site_clusters_uc_{}.h5".format(seqID)
    cluster_h5_file = os.path.join(work_dir, cluster_h5_base)
    cluster_h5_key = "binding_sites_3d/sequence_clusters/{}".format(cluster_h5_base)
    cluster_files["uc"].to_hdf(str(cluster_h5_file), "table", complevel=9, complib="bzip2")
    eppic_interfaces_store.write_output_file(cluster_h5_file, cluster_h5_key)
    safe_remove(cluster_h5_file)

    Parallel(n_jobs=cores)(delayed(save_cluster_msa_file)(f) for f in cluster_files["msaout"])

# def cluster_binding_site_structures(job, distance, cathcode=None, work_dir=None, output_dir=None):
#     for distance in ["maxsub", "rmsd", "tm"]:
#         job.addChildJobFn(distance, cathcode=cathcode, work_dir=work_dir, output_dir=output_dir)
#
def cluster_binding_site_structures(job, cathFileStoreID, cathcode=None, work_dir=None, output_dir="", cores=1, memory="1G"):
    #Run max cluster on all binding site structures
    #Save clusters
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    import psutil
    RealtimeLogger.info("VIRTUAL MEM {}".format(psutil.virtual_memory()))

    eppic_interfaces_store = FileS3IOStore("us-east-1", "eppic-interfaces", file_path_dir=output_dir)

    # if isinstance(cathcode, (list, tuple)):
    #     key = "/".join(cathcode)
    # elif isinstance(cathcode, str):
    #     key = cathcode.replace(".", "/")
    # else:
    #     key = ""
    bs_path = lambda bs: os.path.join(output_dir, "eppic-interfaces", bs)
    key = ""
    all_binding_sites = [bs_path(bs) for bs in \
        eppic_interfaces_store.list_input_directory("binding_sites_3d/"+key) if \
        bs.endswith(".pdb") and os.path.isfile(bs_path(bs)) and os.path.getsize(bs_path(bs))>0]

    if len(all_binding_sites) == 0:
        RealtimeLogger.info("ENDING MaxCluster since empty")
        return

    RealtimeLogger.info("RUNNING MaxCluster one_vs_all")

    run_cath_hierarchy(job, cathcode, cluster_binding_site_structures_in_superfamily,
        cathFileStoreID, all_binding_sites, output_dir=output_dir)

def cluster_binding_site_structures_in_superfamily(job, cathcode, cathFileStoreID, file_list, work_dir=None, output_dir="", cores=20, memory="96G"):
    if work_dir is None:
        work_dir = job.fileStore.getLocalTempDir()

    eppic_interfaces_store = FileS3IOStore("us-east-1", "eppic-interfaces", file_path_dir=output_dir)
    mc = MaxCluster(job=job, work_dir=work_dir, detach=True, intermediate_file_store=eppic_interfaces_store)

    cath_files_prefix = str(os.path.join(output_dir, "eppic-interfaces",
        "binding_sites_3d", cathcode.replace(".", "/")))+"/"
    RealtimeLogger.info("PREFIX {}".format(cath_files_prefix))

    cath_domain_files = [(i, domain_file) for i, domain_file in enumerate(file_list) \
        if domain_file.startswith(cath_files_prefix)][:1]

    RealtimeLogger.info("RUNNING MaxCluster one_vs_all for all {} {}".format(cathcode, cath_domain_files))

    for domain_file in cath_domain_files:
        RealtimeLogger.info("RUNNING MaxCluster one_vs_all for {} {}".format(cathcode, domain_file))
        logFileID = mc.one_vs_all(domain_file, file_list, log=True, cores=cores)

    # log_file = os.path.join(work_dir, "binding_site_clusters.log")
    # dist_file = os.path.join(work_dir, "binding_site_clusters.dist")
    #
    #
    # logFileID = mc.all_vs_all(all_binding_sites, C=1, P=10, distance="maxsub",
    #     R=dist_file, log=log_file, distributed=True, sequence_independent=True,
    #     bb=True, cores=cores, memory="72G")
    # return

    # for f in logFileID:
    #     eppic_interfaces_store.write_output_file(f, "binding_sites_3d/mc_sib/{}".format(os.path.basename(f)))

    #logFileID = mc.all_vs_all(all_binding_sites, C=1, P=10, log=log_file, distributed=True)

    #job.addFollowOnJobFn(save_clustered_binding_sites, logFileID, output_dir=output_dir)

def save_clustered_binding_sites(job, logFileID, output_dir=""):
    work_dir = job.fileStore.getLocalTempDir()
    eppic_interfaces_store = FileS3IOStore("us-east-1", "eppic-interfaces", file_path_dir=output_dir)

    log_file = job.fileStore.readGlobalFile(logFileID, cache=True)

    #Save raw log file
    log_key = "binding_sites_3d/binding_site_clusters.log"
    eppic_interfaces_store.write_output_file(log_file, log_key)

    #Parse log file into clusters
    max_cluster = MaxCluster(job=job, work_dir=work_dir, intermediate_file_store=eppic_interfaces_store)
    binding_site_clusters = max_cluster.get_clusters(log_file)

    #Save dataframe
    bs_cluster_file = os.path.join(work_dir, "binding_site_clusters.h5")
    bs_cluster_key = "binding_sites_3d/binding_site_clusters.h5"
    binding_site_clusters.to_hdf(bs_cluster_file, "table")
    eppic_interfaces_store.write_output_file(bs_cluster_file, bs_cluster_key)

    job.fileStore.deleteGlobalFile(logFileID)

def process_superfamily(job, superfamily, cathFileStoreID, output_dir="",
  force=False, further_parallize=False):
    cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
    cathcode = superfamily.replace("/", ".")

    # cath_domains = filter_hdf_chunks(
    #     cath_file,
    #     "table",
    #     columns=["cath_domain"],
    #     drop_duplicates=True,
    #     cathcode=cathcode)["cath_domain"].tolist()

    avail_domains = [os.path.basename(domain).split("_")[0] for domain in \
        data_stores.eppic_interfaces.list_input_directory("cath/"+superfamily) \
        if domain.endswith("_ddi.h5")]

    #RealtimeLogger.info("Available Doms {}".format(avail_domains))

    if not force:
        done_files = list(data_stores.eppic_interfaces.list_input_directory("binding_sites_3d/"+superfamily))
        # done_pdbs = [os.path.basename(domain).split("_")[0] for domain in \
        #     done_files if domain.endswith(".pdb")] #Doesn't make sense
        done_fasta = [os.path.splitext(os.path.basename(domain))[0] for domain in \
            done_files if domain.endswith(".fasta")]
        done_domains = set(done_fasta) #set(done_pdbs).intersection(set(done_fasta))

        n_domains = len(avail_domains)
        cath_domains = list(set(avail_domains)-done_domains) #avail_domains #
            #list((set(cath_domains)-set(done_domains)).intersection(set(avail_domains)))
        RealtimeLogger.info("Running {}/{} domains from {}".format(len(cath_domains), n_domains, cathcode))
    else:
        RealtimeLogger.info("Running {} domains from {}".format(len(cath_domains), cathcode))

    if further_parallize:
        map_job(job, extract_binding_sites, cath_domains,
            superfamily, cathFileStoreID, output_dir=output_dir)
    else:
        loop_job(job, extract_binding_sites, cath_domains,
            superfamily, cathFileStoreID, output_dir=output_dir)

    job.addFollowOnJobFn(finish_group, "binding_sites_3d/{}".format(superfamily),
        "file-aws:us-east-1:eppic-interfaces:{}".format(output_dir))


def start_toil(job, cathFileStoreID, cathcode=None, update_features=None, seq=False, no_struct=False, output_dir="", force=False):
    RealtimeLogger.info("started") #slurm-13374299.out

    # cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)
    #
    # cath_domains = filter_hdf_chunks(
    #     cath_file,
    #     "table",
    #     columns=["cath_domain", "cathcode"],
    #     drop_duplicates=True)
    eppic_interfaces_store = data_stores.eppic_interfaces_sync(output_dir)

    if not eppic_interfaces_store.exists("binding_sites_3d/all_binding_sites.fasta"):
        #Start CATH hiearchy
        run_cath_hierarchy(job, cathcode, process_superfamily, cathFileStoreID,
            output_dir=output_dir, force=force)

    job.addFollowOnJobFn(cluster_binding_sites, cathFileStoreID, seq=seq, no_struct=no_struct, cathcode=cathcode, output_dir=output_dir)


if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--seq", default=False, action="store_true")
    parser.add_argument("--no-struct", default=False, action="store_true")
    parser.add_argument("--force", default=False, action="store_true")
    parser.add_argument("--output-dir", default=getcwd())
    options = parser.parse_args()
    options.logLevel = "DEBUG"
    options.clean = "always"
    options.targetTime = 1

    if not os.path.isfile("cath.h5"):
        store = IOStore.get("aws:us-east-1:Prop3D-cath")
        store.read_input_file("cath-domain-description-file-small.h5", "cath.h5")

    with Toil(options) as workflow:
        cathFileURL = 'file://' + os.path.abspath("cath.h5")
        cathFileID = workflow.importFile(cathFileURL)
        workflow.start(Job.wrapJobFn(start_toil, cathFileID, seq=options.seq,
            no_struct=options.no_struct, output_dir=options.output_dir, force=options.force))
