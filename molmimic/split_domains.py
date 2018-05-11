import os, sys
sys.path.append("/data/draizene/pdb-tools")
sys.path.append("/data/draizene/molmimic")

import re
import subprocess
import shutil
from collections import defaultdict

import pandas as pd

from molmimic.calculate_features import SwarmJob
from molmimic.map_residues import map_residues
from molmimic.util import atof, natural_keys, to_int

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

def add_sdi(cdd_sdi, cdd_ibis, cleanup=True):
    name_prefix = os.path.splitext(os.path.basename(cdd_ibis))[0]
    cdd_sdi = pd.read_csv(cdd_sdi, dtype = {'pdb': str, 'chain': str})
    cdd_ibis = pd.read_table(cdd_ibis, header=None, names=["pdb", "chain", "residues", "is_multimer", "domain"], dtype = {'pdb': str, 'chain': str})

    all_sdi_domains = cdd_sdi.groupby(["pdbId", "chnLett"])

    #print all_sdi_domains.groups
    #import pdb; pdb.set_trace()

    domains = defaultdict(set)

    for i, row in cdd_ibis.iterrows():
        try:
            sdi_domains = all_sdi_domains.get_group((row["pdb"], row["chain"]))
        except KeyError, IndexError:
            try:
                sdi_domains = all_sdi_domains.get_group((row["pdb"], ""))
            except KeyError, IndexError:
                print "No SDIs for", row["pdb"], row["chain"]
                continue

        #resi = sorted(row["residues"].split(","), key=natural_keys)

        #sort
        resi = [pdbnum for pdbnum in map_residues(row["pdb"], row["chain"], [row["residues"]]) if pdbnum is not None]

        #Check to see if sdis contain binding site
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
        domains[(row["pdb"], row["chain"], sdi_domain1["sdi"], sdi_domain1["domNo"])] |= set(resi_in_sdi)

    if len(domains) == 0:
        print "Error no domains found for", cdd_sdi, cdd_ibis
        return

    outpath = "{}.sdi.tsv".format(os.path.splitext(cdd_ibis)[0])
    with open(outpath, "w") as f:
        print >> f, "pdb\tchain\tsdi\tdomNum\tresi"
        for (pdb, chain, sdi, domNum), resi in domains.iteritems():
            resi = ",".join(map(str, sorted(resi, key=natural_keys)))
            print >> f, "{}\t{}\t{}\t{}\t{}".format(pdb, chain, sdi, domNum, resi)

    if cleanup:
        os.remove(cdd_ibis)
        shutil.move(outpath, cdd_ibis)

def load_cdd(cdd_dir, cleanup=True, job_name="add_sdi", dependency=None):
    import glob
    sdi_cdd_dir = "/data/draizene/molmimic/data/MMDB/cdd"
    job_name = os.environ.get("SLURM_JOB_NAME", job_name)
    job = SwarmJob(job_name, walltime="2:00:00")

    i = 0
    for i, f in enumerate(glob.glob(cdd_dir+"/*.tsv")):
        print i, f
        cdd = os.path.splitext(os.path.basename(f))[0]
        sdi_cdd = os.path.join(sdi_cdd_dir, "{}.csv".format(cdd))
        job += "python {} run {} {} {}\n".format(__file__, sdi_cdd, f, int(cleanup))

    job.run(dependency=dependency, update_dependencies=True)

def submit_cdd(cdd_dir, job_name="add_sdi", dependency=None):
    job = SwarmJob(job_name, walltime="96:00:00", individual=True)
    job += "python {} {}\n".format(__file__, cdd_dir)
    job_id = job.submit_individual(dependency=dependency)
    return job_id

if __name__ == "__main__":
    if len(sys.argv) == 2:
        #Path to direcotry of Ibis split CDD
        load_cdd(sys.argv[1])
    elif len(sys.argv) == 5 and "run" in sys.argv:
        #Path to sdi by cdd file and ibis split cdd
        add_sdi(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print "Error:", len(sys.argv), sys.argv
