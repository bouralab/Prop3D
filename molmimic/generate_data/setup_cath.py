from __future__ import print_function
import re
from collections import defaultdict
import pandas as pd
from joblib import Parallel, delayed

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
    "architechture": int,
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

def parse_cath_chunk(chunk, cath_file, line_start, line_end):
    num_domains = 0
    nseg = 0
    current_domain = None
    parsing = False

    srange_re = re.compile("START=([a-zA-Z0-9\-]+)\s+STOP=([a-zA-Z0-9\-]+)")

    all_domains = pd.DataFrame()
    def append_domain(domain):
        assert set(domain.keys())==set(cath_desc.keys())
        all_domains = all_domains.append(domain)
        if all_domains.shape[0] >= 10000:
            save_domains()
            all_domains = pd.DataFrame()

    def save_domains():
        all_domains.to_hdf(str("cath-domain-description-file.h5.{}".format(chunk)), "table", format="table",
            table=True, complevel=9, complib="bzip2", min_itemsize=1024,
            append=True, mode="a",
            data_coumns=["cath_domain", "pdb", "chain", "domain",
                         "cathcode", "class", "arch", "topology", "homology"])

    with open(cath_file) as cath:
        for i, line in enumerate(cath):
            if i < line_start or i > line_end+1:
                continue

            line = line.rstrip()
            if line.startswith("FORMAT"):
                parsing = True
                current_domain = defaultdict(str)
                #if num_domains%50==0: print(num_domains, "of", num_entries)
                num_domains += 1
                continue
            if parsing and line.startswith("//"):
                nseg = 0
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
                elif field == "NSEGMENTS":
                    _current_domain = current_domain.copy()
                    continue
                elif field == "ENDSEG":
                    current_domain["nseg"] = nseg+1
                    append_domain(current_domain)
                    current_domain = _current_domain.copy()
                    nseg += 1
                elif field == "SRANGE":
                    m = srange_re.match(value)
                    if m:
                        current_domain["srange_start"] = m.group(1)
                        current_domain["srange_stop"] = m.group(2)

                else:
                    current_domain[field.lower()] += value
        save_domains()
        print("Finished 100")

def parse_cath(cath_file):
    """
    """
    with open(cath_file) as cath:
        num_entries = sum(1 for l in cath if l.startswith("ENDSEG"))

    indices = []
    with open(cath_file) as cath:
        num_seen = 0
        start_line = 0
        for i, line in enumerate(cath):
            if line.startswith("//"):
                if num_seen == 99:
                    indices.append((start_line, i))
                    start_line = i+1
                else:
                    num_seen += 1

    Parallel(n_jobs=20)(delayed(parse_cath_chunk)(i, cath_file, s, e) for i, (s, e) in enumerate(indices))

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
            row["srange_start"], row["srange_stop"]), axis=1)
        merged.groupby(["from", "to"])["overlap"].idxmax()




if __name__ == "__main__":
    parse_cath("cath-domain-description-file.txt")
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
