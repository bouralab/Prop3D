import os
import re
import glob
import logging
import shutil
from collections import defaultdict
from itertools import groupby
from multiprocessing.pool import Pool


import requests
import pandas as pd

from botocore.exceptions import ClientError

from Prop3D.util.iostore import IOStore
from Prop3D.generate_data.job_utils import map_job
from Prop3D.util.cath import run_cath_hierarchy, download_cath_domain
from Prop3D.util.hdf import get_file, filter_hdf, filter_hdf_chunks

from toil.realtimeLogger import RealtimeLogger

def download_superfamily(job, sfam_id, cathFileStoreID, download_all=False,
  check_if_downloaded=False, skip_ids=None, further_parallelize=True):
    work_dir = job.fileStore.getLocalTempDir()
    cath_file = job.fileStore.readGlobalFile(cathFileStoreID, cache=True)

    RealtimeLogger.info("Filter CATH")
    cath_domains = filter_hdf(cath_file, "table", drop_duplicates=True,
        columns=["cath_domain"], cathcode=sfam_id)["cath_domain"]
    RealtimeLogger.info("Filtered CATH")
    RealtimeLogger.info("Filtered CATH {}".format(len(cath_domains)))
    #RealtimeLogger.info("CATH DOMAIN: {}".format(cath_domains))

    cath_store = IOStore.get("aws:us-east-1:Prop3D-cath-structure")

    if not check_if_downloaded:
        done_domains = set(os.path.basename(k)[:-4] for k in \
            cath_store.list_input_directory(sfam_id.replace(".", "/")) \
            if k.endswith(".pdb"))
        cath_domains = cath_domains[~cath_domains.isin(done_domains)]

    if skip_ids is not None and isinstance(skip_ids, dict) and sfam_id in skip_ids:
        cath_domains = cath_domains[~cath_domains.isin(skip_ids[sfam_id])]

    RealtimeLogger.info("Start dl {}".format(len(cath_domains)))

    cath_domains = cath_domains.iloc[:1]

    if further_parallelize:
        map_job(job, download_cath_domain, cath_domains, sfam_id, check)
    else:
        for cath_domain in cath_domains:
            RealtimeLogger.info("Running {}".format(cath_domain))
            try:
                download_cath_domain(job, cath_domain, sfam_id, check=check, work_dir=work_dir)
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                import traceback
                tb = traceback.format_exc()
                fail_file = os.path.join(work_dir, "{}.fail".format(cath_domain))
                fail_key = "{}/{}.fail".format(sfam_id.replace(".", "/"), cath_domain)
                with open(fail_file, "w") as fh:
                    fh.write(tb)
                cath_store.write_output_file(fail_file, fail_file)
                try:
                    os.remove(fail_file)
                except OSError:
                    pass
                RealtimeLogger.info("Failed running CATH DOMAIN {}/{}:".format(
                    sfam_id, cath_domain))
                RealtimeLogger.info(tb)
                continue

cath_desc = {
    "cath_domain":str,
    "pdb":str,
    "chain": str,
    "domain": str,
    "version": str,
    "verdate": str,
    "name": str,
    "source": str,
    "cathcode": str,
    "class": int,
    "architecture": int,
    "topology": int,
    "homology": int,
    "class_name": str,
    "architechture_name": str,
    "topology_name": str,
    "homology_name": str,
    "dlength": int,
    "dseqh": str,
    "dseqs": str,
    "nsegments": int,
    "nseg": int,
    "segment": str,
    "srange_start": str,
    "srange_stop": str,
    "slength": int,
    "sseqh": str,
    "sseqs": str
}

def parse_cath_chunk(job, chunk_info, cath_file, chunk_size=100):
    work_dir = job.fileStore.getLocalTempDir()

    if len(chunk_info) == 2 and len(chunk_info[1]) == 2:
        RealtimeLogger.info("Chunk info {}".format(chunk_info))
        chunk, (line_start, line_end) = chunk_info
        is_end = False
    elif len(chunk_info) == 2 and len(chunk_info[1]) == 3:
        chunk, (line_start, line_end, last_count) = chunk_info
        is_end = True
    else:
        raise RuntimeError("Invalid args")

    #chunk, (line_start, line_end) = chunk_info
    num_domains = 0
    nseg = 0
    current_domain = None
    _current_domain = None #Temp variable for multiple segments
    parsing = False

    srange_re = re.compile("START=([a-zA-Z0-9\-]+)\s+STOP=([a-zA-Z0-9\-]+)")

    cath_file = get_file(job, "cath.txt", cath_file, work_dir=work_dir)

    all_domains = pd.DataFrame()

    with open(cath_file) as cath:
        for i, line in enumerate(cath):
            if i < line_start:
                continue
            if i > line_end:
                break
            #print(line.rstrip())
            line = line.strip()
            if line.startswith("FORMAT"):
                print("START")
                parsing = True
                current_domain = defaultdict(str)
                #if num_domains%50==0: print(num_domains, "of", num_entries)
                num_domains += 1
                continue
            if parsing and line.startswith("//"):
                nseg = 0
                download_cath_domain(job, current_domain, work_dir)
                continue
            if parsing:
                field, value = line[:10].rstrip(), line[10:]
                if field == "DOMAIN":
                    current_domain["cath_domain"] = value
                    current_domain["pdb"] = value[:4]
                    current_domain["chain"] = value[4]
                    current_domain["domain"] = int(value[5:])
                elif field == "CATHCODE":
                    current_domain["cathcode"] = value
                    c, a, t, h = value.split(".")
                    current_domain["class"] = int(c)
                    current_domain["architechture"] = int(a)
                    current_domain["topology"] = int(t)
                    current_domain["homology"] = int(h)

                elif field == "CLASS":
                    current_domain["class_name"] += value
                elif field == "ARCH":
                    current_domain["architechture_name"] += value
                elif field == "TOPOL":
                    current_domain["topology_name"] += value
                elif field == "HOMOL":
                    current_domain["homology_name"] += value
                elif field in ["DLENGTH", "NSEGMENTS", "SLENGTH"]:
                    current_domain[field.lower()] = int(value)
                #elif field == "NSEGMENTS":
                #    _current_domain = current_domain.copy()
                #    continue
                elif field == "ENDSEG":
                    current_domain["nseg"] = nseg+1
                    #append_domain(all_domains, current_domain)
                    #current_domain = _current_domain
                    all_domains = all_domains.append(pd.Series(current_domain), ignore_index=True)
                    nseg += 1
                    print(chunk, all_domains.shape[0])
                elif field == "SRANGE":
                    m = srange_re.match(value)
                    current_domain["srange_start"] = m.group(1) if m else np.NaN
                    current_domain["srange_stop"] = m.group(2) if m else np.NaN

                else:
                    current_domain[field.lower()] += value

    correct_num = last_count if is_end else chunk_size
    assert all_domains["cath_domain"].drop_duplicates().shape[0] == correct_num

    cath_chunk_file = os.path.join(work_dir, "cath-domain-description-file.h5.{}".format(chunk))
    all_domains.to_hdf(cath_chunk_file, "table", format="table",
        table=True, complevel=9, complib="bzip2", min_itemsize=1024,
        data_coumns=["cath_domain", "pdb", "chain", "domain",
                     "cathcode", "class", "arch", "topology", "homology"])

    in_store = IOStore.get("aws:us-east-1:Prop3D-cath")
    in_store.write_output_file(cath_chunk_file, os.path.basename(cath_chunk_file))

    for f in (cath_file, cath_chunk_file):
        try:
            os.remove(f)
        except OSError:
            pass

def get_cath_split_jobs(cath_file, chunk_size=100):
    """
    """
    indices = []
    total_seen = 0
    with open(cath_file) as cath:
        num_seen = 0
        start_line = 0
        for i, line in enumerate(cath):
            if line.startswith("//"):
                if num_seen == chunk_size-1:
                    indices.append((start_line, i+1))
                    start_line = i+1
                    total_seen += 1
                    num_seen = 0
                else:
                    num_seen += 1

    if num_seen > 0:
        indices.append((start_line, i+1, num_seen))

    return indices

def create_small_description_file(job):
    work_dir = job.fileStore.getLocalTempDir()
    in_store = IOStore.get("aws:us-east-1:Prop3D-cath")

    small_desc_file = os.path.join(work_dir, "cath-domain-description-file-small.h5")
    for cath_chunk_key in in_store.list_input_directory("cath-domain-description-file.h5."):
        cath_chunk_file = os.path.join(work_dir, cath_chunk_key)
        in_store.read_input_file(cath_chunk_key, cath_chunk_file)

        chunk_id = int(cath_chunk_file.split(".")[-1])
        cath_chunk = pd.read_hdf(cath_chunk_file, "table", mode="r")
        cath_chunk = cath_chunk[[
            "cath_domain", "pdb", "chain", "domain",
            "cathcode", "class", "architechture", "topology", "homology",
            "nsegments", "nseg", "srange_start", "srange_stop", "slength"]]
        cath_chunk = cath_chunk.assign(chunk=chunk_id)
        cath_chunk.to_hdf(small_desc_file, "table",
            format="table", table=True, mode="a", append=True,
            complevel=9, complib="bzip2", min_itemsize=1024,
            data_columns=["cath_domain", "cathcode"])
        del cath_chunk

        try:
            os.remove(cath_chunk_file)
        except OSError:
            pass

    in_store.write_output_file(small_desc_file, os.path.basename(small_desc_file))

def start_toil(job, cathFileID, download_all=False, check=False, skip_ids=None):
    work_dir = job.fileStore.getLocalTempDir()
    in_store = IOStore.get("aws:us-east-1:Prop3D-cath")

    if download_all:
        RealtimeLogger.info("Downlaod all, Check? {}".format(check))
        cath_file = job.fileStore.readGlobalFile(cathFileID, cache=True)

        cathcodes = filter_hdf_chunks(cath_file, "table", columns=["cathcode"],
            drop_duplicates=True)["cathcode"]

        #cathcodes = ["1.10.530.40"]

        map_job(job, download_sfam, cathcodes, cathFileID, download_all, check, skip_ids)

        return

    sdoms_file = os.path.join(work_dir, "cath-domain-description-file.txt")
    in_store.read_input_file("cath-domain-description-file.txt", sdoms_file)
    cathFileStoreID = job.fileStore.writeGlobalFile(sdoms_file)

    chunks = [(i, chunk) for i, chunk in enumerate(get_cath_split_jobs(sdoms_file)) \
        if not in_store.exists("cath-domain-description-file.h5.{}".format(i))]

    map_job(job, parse_cath_chunk, chunks, cathFileStoreID)

    job.addFollowOnJobFn(create_small_description_file)

    try:
        os.remove(sdoms_file)
    except OSError:
        pass

def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

def map_cath_cdd(job, sfam_id, pdbFileStoreID, cathFileStoreID):
    cdd_file = get_file(job, "PDB.h5", pdbFileStoreID)
    cdd_domains = filter_hdf(
        pdb_file, "merged",
        columns = ["sdi", "domNo", "gi", "pdbId", "chnLett", "from", "to"],
        sfam_id = sfam_id).drop_duplicates()

    cath_file = get_file(job, "cath.h5", pdbFileStoreID)
    cath_domains = pd.read_hdf(cath_file, "table", where="pdb in cdd_domain[pdbId] and chain in cdd_domain[chnLett]")

    for pdb_chain, group in cdd_domains.groupby(["sdi"]):
        merged = pd.merge(cdd_domains, cath_domains, left_on=["pdbId", "chnLett"], right_on=["pdb", "chain"])
        merged["overlap"] = merged.apply(lambda row: overlap(row["from"], row["to"],
            row["srange_start"], row["srange_stop"]), axis=0)
        merged.groupby(["from", "to"])["overlap"].idxmax()

if __name__ == "__main__":
    from toil.common import Toil
    from toil.job import Job

    parser = Job.Runner.getDefaultArgumentParser()
    parser.add_argument("--download-all", default=False, action="store_true")
    parser.add_argument("--check", default=False, action="store_true")
    parser.add_argument("--skip-ids", default=None)

    options = parser.parse_args()
    options.logLevel = "INFO"
    options.clean = "always"
    options.targetTime = 1

    if options.skip_ids is not None:
        with open(options.skip_ids) as fh:
            options.skip_ids = {sfam_id.replace("/", "."): list(pdb_keys) \
                for sfam_id, pdb_keys in groupby(sorted(fh), \
                key=lambda k: k.rsplit("/", 1)[0])}


    detector = logging.StreamHandler()
    logging.getLogger().addHandler(detector)

    if options.download_all and not os.path.isfile("cath.h5"):
        in_store = IOStore.get("aws:us-east-1:Prop3D-cath")
        in_store.read_input_file("cath-domain-description-file-small.h5", "cath.h5")

    with Toil(options) as workflow:
        if options.download_all:
            cathFileURL = 'file://' + os.path.abspath("cath.h5")
            cathFileID = workflow.importFile(cathFileURL)
        else:
            cathFileID = None
        workflow.start(Job.wrapJobFn(start_toil, cathFileID, options.download_all, False, options.skip_ids))


# if __name__ == "__main__":
#     indicies = get_cath_split_jobs("cath-domain-description-file.txt")
#     if len(indices) > 0:
#         Parallel(n_jobs=20)(delayed(parse_cath_chunk)((i, s, e), cath_file) for i, (s, e) in enumerate(indices))

    # from toil.common import Toil
    # from toil.job import Job
    #
    # parser = Job.Runner.getDefaultArgumentParser()
    # options = parser.parse_args()
    # options.clean = "always"
    # options.targetTime = 1
    #
    # job = Job.wrapJobFn(start_toil)
    # with Toil(options) as workflow:
    #     pdbs = [(os.path.basename(f), workflow.importFile('file://' + f)) for f in glob.glob("/root/ig/*/*.pdb")]
    #     workflow.start(Job.wrapJobFn(start_toil, pdbs))







#
# def get_seqres():
#     seqres = pd.read_table("cath-domain-boundaries-seqreschopping.txt", header=None, names=["cath_id", "ranges"])
#     seqres.loc[:, ["pdb", "chain", "domNo"]] = seqres["cath_id"].apply(
#         lambda pcd: pd.Series({"pdb":pcd[:4], "chain":pcd[4], "domNo":int(pcd[5:])}))
#     return seqres
#
# def group_cath(sfam_id):
#     names = pd.read_table("cath-names.txt", header=None, names=["node", "cath_id", "name"])
#
#     all_sdoms = pd.read_hdf(os.path.join(data_path_prefix, "PDB.h5"), "merged")
#     sdoms_groups = all_sdoms.groupby("sfam_id")
#
#     try:
#         sdoms = sdoms_groups.get_group(sfam_id).copy()
#     except KeyError:
#         job.log("No structural domains for {}".format(cdd))
#         return
#
#     del sdoms_groups
#     sdoms = sdoms.rename(columns={"pdbId":"pdb", "chnLett":"chain"})
#
#     sdoms = pd.merge(sdoms, get_seqres(), on=["pdb", "chain", "domNo"], how="left")
#
#     cath_node = pd.merge(sdoms, names, on="cath_id")["node"].iloc[0]
