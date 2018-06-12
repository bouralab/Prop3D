import os, sys
sys.path.append("/data/draizene/pdb-tools")
sys.path.append("/data/draizene/molmimic")

import re
import subprocess
import shutil
from collections import defaultdict
from itertools import combinations

import pandas as pd

from molmimic.calculate_features import SwarmJob
from molmimic.map_residues import mmdb_to_pdb_resi
from molmimic.util import atof, natural_keys, to_int, get_interfaces_path

def overlap(min1, max1, min2, max2):
    return max(0, min(max1, max2) - max(min1, min2))

def split_domains_in_pdb(pdb_file, pdb_name, chain, sdi=None, domainNum=None, sdi_pdb_file=None):
    if sdi_pdb_file is None:
        sdi_pdb_file = "/data/draizene/molmimic/data/MMDB/Structural_Domains.csv"

    pdb_by_sdi = pd.read_csv(sdi_pdb_file, dtype = {'pdb': str, 'chain': str})
    pdb_path = os.path.dirname(pdb_file)

    domains = pdb_by_sdi[(pdb_by_sdi["pdbId"]==pdb_name.upper())&((pdb_by_sdi["chnLett"]==chain)|(pdb_by_sdi["chnLettPrefix"]==""))]

    #import pdb; pdb.set_trace()

    for i, domain in domains.iterrows():
        print i
        print domain["sdi"], domain["domNo"]
        print "{}:{}".format(domain["frm"], domain["tu"])
        #Split domain from PDB file

        if sdi is not None and domain["sdi"] != sdi: continue
        if domainNum is not None and domain["domNo"] != domainNum: continue

        domain_file = os.path.join(pdb_path, "{}_{}_sdi{}_d{}.pdb".format(pdb_name, chain, domain["sdi"], domain["domNo"]))
        with open(domain_file, "w") as domain_f:
            subprocess.Popen(["/data/draizene/pdb-tools/pdb_rslice.py",
                "{}:{}".format(domain["frm"], domain["tu"]), pdb_file],
                stdout=domain_f)
        yield domain_file, pdb_name, chain, (domain["sdi"], domain["domNo"])

def add_sdi(dataset_name, cdd_ibis, cdd_sdi, cleanup=False):
    interfaces_path = get_interfaces_path(dataset_name)

    if not os.path.exists(interfaces_path):
        os.makedirs(interfaces_path)

    name_prefix = os.path.splitext(os.path.basename(cdd_ibis))[0]
    inpath = cdd_ibis
    outpath = os.path.join(interfaces_path, "{}.sdi".format(name_prefix))
    observed_outpath = os.path.join(interfaces_path, "{}.observed_sdi".format(name_prefix))
    inferred_outpath = os.path.join(interfaces_path, "{}.inferred_sdi".format(name_prefix))
    extended_outpath = os.path.join(interfaces_path, "{}.extended_sdi".format(name_prefix))

    if not os.path.isfile(cdd_ibis) and not os.path.isfile(cdd_sdi):
        print "Files not found"
        with open(outpath, "w") as f, open(observed_outpath, "w") as f2, open(inferred_outpath, "w") as f3, open(extended_outpath, "w") as f3:
            pass
        return

    cdd_sdi = pd.read_csv(cdd_sdi, dtype = {'pdb': str, 'chain': str})
    cdd_ibis = pd.read_table(cdd_ibis, dtype = {'pdb': str, 'chain': str})

    all_sdi_domains = cdd_sdi.groupby(["pdbId", "chnLett"])

    domains = defaultdict(set)
    observedDomains = defaultdict(set)
    inferredDomains = defaultdict(set)

    extended_sdi_file = open(extended_outpath, "w")
    print >> extended_sdi_file, "pdb\tchain\tsdi\tdomNum\tresi\tresn\tis_multimer\tcdd\tpartner\tpdb_evidence\tobserved"

    for i, row in cdd_ibis.iterrows():
        try:
            sdi_domains = all_sdi_domains.get_group((row["pdb"], row["chain"]))
        except KeyError, IndexError:
            try:
                sdi_domains = all_sdi_domains.get_group((row["pdb"], ""))
            except KeyError, IndexError:
                print "No SDIs for", row["pdb"], row["chain"]
                continue

        resi = row["resi"].split(",")
        resn = row["resn"]

        #Check to see if sdis contain binding site
        if sdi_domains.shape[0] == 1:
            sdi_domain1 = next(sdi_domains.iterrows())[1]
        else:
            for i, sdi_domain1 in sdi_domains.iterrows():
                for j, sdi_domain2 in sdi_domains.iterrows():
                    #Remove sdis that contain multiple sdis, only use smallest
                    if sdi_domain1["sdi"] == sdi_domain2["sdi"]:
                        continue

                    if sdi_domain1["frm"] >= sdi_domain2["frm"] and sdi_domain1["tu"] <= sdi_domain2["tu"]:
                        #sdi2 inside sd1 => skip
                        break
                else:
                    #sd1 not found in any other sdi
                    #Check this domain to see if it contains binding site
                    #if to_int(resi[0]) >= sdi_domain1["frm"] and to_int(resi[-1]) <= sdi_domain1["tu"]:
                    #    break
                    sdi_length = float(sdi_domain1["tu"]-sdi_domain1["frm"])
                    if overlap(to_int(resi[0]), to_int(resi[-1]), sdi_domain1["frm"], sdi_domain1["tu"])/sdi_length >= 0.8:
                        break
            else:
                #Binding site not found in any sdi => skip
                continue

        resi_in_sdi = [r for r in resi if sdi_domain1["frm"]<=to_int(r)<=sdi_domain1["tu"]]

        if len(r) == 0:
            continue
        print >> extended_sdi_file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(row["pdb"], row["chain"], sdi_domain1["sdi"], sdi_domain1["domNo"], ",".join(map(str, sorted(resi_in_sdi, key=natural_keys))), row["resn"], row["is_multimer"], row["cdd"], row["partner"], row["pdb_evidence"], row["is_observed"])

        #Combine resides with all residues (infered and observed)
        domains[(row["pdb"], row["chain"], sdi_domain1["sdi"], sdi_domain1["domNo"])] |= set(resi_in_sdi)

        if bool(row["is_observed"]):
            #Combine resides with observed residues (infered and observed)
            observedDomains[(row["pdb"], row["chain"], sdi_domain1["sdi"], sdi_domain1["domNo"])] |= set(resi_in_sdi)
        else:
            inferredDomains[(row["pdb"], row["chain"], sdi_domain1["sdi"], sdi_domain1["domNo"])] |= set(resi_in_sdi)

    extended_sdi_file.close()

    if len(domains) == 0:
        print "Error no domains found for", name_prefix

    for path, combinded_domains in ((outpath, domains), (observed_outpath, observedDomains), (inferred_outpath, inferredDomains)):
        with open(path, "w") as f:
            print >> f, "pdb\tchain\tsdi\tdomNum\tresi"
            for (pdb, chain, sdi, domNum), resi in combinded_domains.iteritems():
                resi = ",".join(map(str, sorted(resi, key=natural_keys)))
                print >> f, "{}\t{}\t{}\t{}\t{}".format(pdb, chain, sdi, domNum, resi)

    if cleanup:
        os.remove(inpath)

def load_cdd(dataset_name, cdd_dir, cleanup=False, job_name="add_sdi", dependency=None):
    import glob
    sdi_cdd_dir = "/data/draizene/molmimic/data/MMDB/cdd"
    job_name = job_name or os.environ.get("SLURM_JOB_NAME", job_name)
    job = SwarmJob(job_name, walltime="2:00:00", mem="15")

    CDD = pd.read_csv("MMDB/StructDomSfam.csv", usecols=["label"]).drop_duplicates().dropna()
    CDD = sorted(CDD["label"].apply(lambda cdd: cdd.replace("/", "").replace("'", "\'")).tolist())

    for cdd in CDD: #enumerate(glob.glob(cdd_dir+"/*.tsv")):
        ibis_cdd = os.path.join(cdd_dir, "{}.raw".format(cdd))
        sdi_cdd = os.path.join(sdi_cdd_dir, "{}.csv".format(cdd))
        job += "python {} run {} {} {}\n".format(__file__, dataset_name, ibis_cdd, sdi_cdd)

    return job.run(dependency=dependency)

def submit_cdd(cdd_dir, job_name="add_sdi", dependency=None):
    job = SwarmJob(job_name, walltime="96:00:00", individual=True)
    job += "python {} {}\n".format(__file__, cdd_dir)
    job_id = job.submit_individual(dependency=dependency)
    return job_id

if __name__ == "__main__":
    if len(sys.argv) == 3:
        #Path to direcotry of Ibis split CDD
        load_cdd(sys.argv[1], sys.argv[1])
    elif len(sys.argv) == 5 and "run" in sys.argv:
        #Path to sdi by cdd file and ibis split cdd
        add_sdi(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print "Error:", len(sys.argv), sys.argv
