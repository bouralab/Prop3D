import os, sys
import urllib.request
import re
from collections import defaultdict
import pandas as pd
import multiprocessing
from joblib import Parallel, delayed

try:
    import h5pyd as h5py
    DISTRIBUTED = True
except ImportError:
    try:
        import h5py
        DISTRIBUTED = False
    except:
        raise ImportError("h5pyd or h5py must be installed")

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

def parse_cath_chunk(chunk, cath_file, cath_full_h5, line_start, line_end, force=False):
    num_domains = 0
    nseg = 0
    current_domain = None
    parsing = False

    srange_re = re.compile("START=([a-zA-Z0-9\-]+)\s+STOP=([a-zA-Z0-9\-]+)")

    store = h5py.File(cath_full_h5, mode="a", use_cache=False)

    def save_domain(domain):
        assert set(domain.keys())==set(cath_desc.keys())

        if domain["nseg"] > 1:
            #Skip other domain segments
            return

        cath_domain = domain.pop("cath_domain")
        cathcode = domain.pop("cathcode")
        domain = pd.Series(domain)
        metadata = domain[["name", "source", "dseqs", "dlength", "nsegments", "verdate", "version"]]
        metadata = metadata.rename({"name":"description"})
        key = f'/{cathcode.replace(".", "/")}/domains/{cath_domain}'
        print(key, cathcode, cath_domain)

        if not force and cath_domain in store[f'/{cathcode.replace(".", "/")}/domains']:
            return

        try:
            group = store.require_group(key)
        except ValueError:
            print("FAILED", key, cathcode, cath_domain)
            raise

        for key, value in metadata.items():
            group.attrs[key] = value

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
                elif field in ["DLENGTH", "SLENGTH"]:
                    current_domain[field.lower()] = int(value)
                elif field == "NSEGMENTS":
                    current_domain[field.lower()] = int(value)
                    _current_domain = current_domain.copy()
                    continue
                elif field == "ENDSEG":
                    current_domain["nseg"] = nseg+1
                    save_domain(current_domain)
                    current_domain = _current_domain.copy()
                    nseg += 1
                elif field == "SRANGE":
                    m = srange_re.match(value)
                    if m:
                        current_domain["srange_start"] = m.group(1)
                        current_domain["srange_stop"] = m.group(2)

                else:
                    current_domain[field.lower()] += value
    store.close()

def parse_cath(cath_file, cath_full_h5):
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
    n_jobs = 10 if DISTRIBUTED else 1
    Parallel(n_jobs=20)(delayed(parse_cath_chunk)(i, cath_file, cath_full_h5, s, e) for i, (s, e) in enumerate(indices) if i>500)

def parse_cath_names_for_group(cath_full_h5, name, group_df):
    store = h5py.File(cath_full_h5, mode="a", use_cache=False)

    print("Start", name)

    for _, row in group_df.iterrows():
        group = store.create_group(row.h5_key)
        group.description = row.description
        group.representativeDomain = row.representative
        if row.cath_code.count(".") == 3:
            store.create_group(f"{row.h5_key}/domains")

    store.close()

    print("End", name)

def parse_cath_names(cath_names_file, cath_full_h5):

    #try:
    #    store = h5py.File(cath_full_h5, mode="r", use_cache=False)
    #Empty file
    #store = h5py.File(cath_full_h5, mode="w", use_cache=False)
    #store.close()

    names = pd.read_csv(cath_names_file, sep="    ", header=None, comment="#",
        names=["cath_code", "representative", "description"])
    names["description"] = names["description"].str[1:]
    names = names.assign(h5_key="/"+names["cath_code"].str.replace(".","/"))
    names = names.assign(group=names["cath_code"].str.split(".", expand=True)[[0,1]].fillna("").agg('.'.join, axis=1))

    root_nodes = names[names["cath_code"].str.split(".").agg(len)==1]
    parse_cath_names_for_group(cath_full_h5, "root", root_nodes)

    groups = names[~names.index.isin(root_nodes.index)].groupby("group")
    Parallel(n_jobs=20)(delayed(parse_cath_names_for_group)(cath_full_h5, name, group) for \
        name, group in groups)


if __name__ == "__main__":
    force = len(sys.argv)>1 and args[1] in ["-f", "--force"]

    if not os.path.isfile("cath-names.txt"):
        urllib.request.urlretrieve(
            "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-names.txt",
            "cath-names.txt")

    if not os.path.isfile("cath-domain-description-file.txt"):
        urllib.request.urlretrieve(
            "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-description-file.txt",
            "cath-domain-description-file.txt")


    cath_full_h5 = "/home/cath-full.h5"
    #parse_cath_names("cath-names.txt", cath_full_h5)
    parse_cath("cath-domain-description-file.txt", cath_full_h5)
