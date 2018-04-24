from itertools import groupby
from collections import Counter, defaultdict
import subprocess

def run_query_pair(pdb, chain, domain=None):
    inputfile = "{}_{}_{}_inputfile.txt".format(pdb, chain, domain)
    with open(inputfile, "w") as f:
        inputline = "{} {} {}".format(pdb, chain)
        print >> f, "{}{}".format(inputline, " "*(8-len(inputline)))

    alignment_file = "{}_{}_{}_alignmentfile.txt".format(pdb, chain, domain)
    subprocess.call(["/netopt/structure/bin/queryPair", 5, inputfile, "/dev/null", alingment_file])
    return alignment_file

def parse_vast(vast_file):
    for line in vast_file:
        fields = line.rstrip().split("\t")
        fields[0] = tuple(fields[0].split(" "))
        fields[1] = tuple(fields[1].split(" "))
        yield fields

def vast_conservation_score(pdb, chain, domain=0, psuedocount=1., cutoff=None):
    alignment = run_query_pair(pdb, chain, domain=domain)

    scores = defaultdict(lambda: defaultdict(lambda: psuedocount))
    number_alns = Counter()
    with open(alignment) as aln_file:
        for (domainA, domainB), alignments in groupby(parse_vast(aln_file), key=lambda x: x[:2]):
            number_alns[domainA] += 1
            for aln in alignments:
                start, end = aln[3:5]
                for resi in xrange(start-1, end):
                    scores[domainA][resi] += 1

    if cutoff is not None:
        calc_score = lambda domain, resi, total_aln: int(scores[domain][resi]/float(total_aln)>=cutoff)
    else:
        calc_score = lambda domain, resi, total_aln: scores[domain][resi]/float(total_aln)

    for domain, residues in scores.iteritems():
        num_domain_alns = number_alns[domain]
        for resi in residues:
            scores[domain][resi] = calc_score(domain, resi, num_domain_alns)

    return scores
