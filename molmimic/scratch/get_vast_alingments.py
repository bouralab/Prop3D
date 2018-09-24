import os, sys
from itertools import groupby
from collections import Counter, defaultdict
import subprocess

import pandas as pd

def run_query_pair(pdb, chain, domain):
    inputfile = "{}_{}_d{}_inputfile.txt".format(pdb, chain, domain)
    with open(inputfile, "w") as f:
        inputline = "{} {} {}".format(pdb, chain, domain)
        print >> f, "{}{}".format(inputline, " "*(8-len(inputline)))

    alignment_file = "/data/draizene/molmimic/data/MMDB/VAST/{}/{}_{}_d{}_alignmentfile".format(pdb.lower()[1:3], pdb, chain, domain)
    print " ".join(["/netopt/structure/bin/queryPair", "1", inputfile, "/dev/null", alignment_file])
    subprocess.call(["/netopt/structure/bin/queryPair", "1", inputfile, "/dev/null", alignment_file])

    os.remove(inputfile)

    return alignment_file

def parse_vast(vast_file):
    for line in vast_file:
        fields = line.rstrip().split("\t")
        fields[0] = tuple(fields[0].split(" "))
        fields[1] = tuple(fields[1].split(" "))
        yield fields

def vast_conservation_score(pdb, chain, sdi, domain, psuedocount=1., cutoff=0.6, save_to_file=False):
    alignment = run_query_pair(pdb, chain, domain)

    scores = defaultdict(lambda: defaultdict(lambda: psuedocount))
    number_alns = Counter()
    with open(alignment) as aln_file:
        for (domainA, domainB), alignments in groupby(parse_vast(aln_file), key=lambda x: x[:2]):
            number_alns[domainA] += 1
            for aln in alignments:
                start, end = aln[3:5]
                for resi in xrange(int(start)-1, int(end)):
                    scores[domainA][resi] += 1

    calc_score_cutoff = lambda domain, resi, total_aln: int(scores[domain][resi]/float(total_aln)>=cutoff)
    calc_score = lambda domain, resi, total_aln: scores[domain][resi]/float(total_aln)

    for domain, residues in scores.iteritems():
        num_domain_alns = number_alns[domain]
        for resi in residues:
            scores[domain][resi] = (
                calc_score(domain, resi, num_domain_alns), 
                calc_score_cutoff(domain, resi, num_domain_alns))

    if save_to_file:
        for (dPdb, dChain, dNum), dScores in scores.iteritems(): 
            with open("/data/draizene/molmimic/data/MMDB/VAST/{}/{}_{}_sdi{}_d{}_scores.tsv".format(dPdb.lower()[1:3], dPdb, dChain, sdi, dNum)) as f:
                print >> f, "index\tscore\tis_conserved"
                dScores = sorted(dScores.iteritems(), key=lambda x: x[0])
                for index, (score, is_conserved) in dScores:
                    print >> f, "{}\t{:.4f}\t{}".format(index, score, is_conserved)

    return scores

def load_sdi_ibis(sdi_ibis_file):
    sdi_ibis = pd.read_table(sdi_ibis_file, header=None, names=["pdb", "chain", "sdi", "domNo", "resi"])
    
    job = SwarmJob("vast_aln")

    pdb_groups = sdi_ibis.groupby(lambda x: sdi_ibis["pdb"].loc[x][1:3])

    i = 0
    for pdb_divided, pdb_group in pdb_groups:
        pdb_path = "/data/draizene/molmimic/data/MMDB/VAST/{}".format(pdb_divided.lower())
        if not os.path.exists(pdb_path):
            os.makedirs(pdb_path)
        for pdb in pdb_group:
            i+=1
            print i, row["pdb"], row["chain"], row["domNo"] 
            job += "python {} run {} {} {} {}\n".format(__file__, row["pdb"], row["chain"], row["sdi"], row["domNo"])

    job.run()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        #Path to direcotry of Ibis split CDD
        load_sdi_ibis(sys.argv[1])
    elif len(sys.argv) == 6 and "run" in sys.argv:
        #Path to sdi by cdd file and ibis split cdd 
        vast_conservation_score(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        print "Error:", len(sys.argv), sys.argv