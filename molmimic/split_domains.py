import re
from collection import defaultdict

def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [ atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text) ]

def to_int(s):
    return int("".join([d for d in s if d.isdigit()]))

def split_domains_pdb(sdi_pdb, pdb, chain):
	pdb_by_sdi = pd.read_csv(sdi_pdb)

	domains = pdb_by_sdi[(pdb_by_sdi["pdbId"]==pdb)&((pdb_by_sdi["chnLettPrefix"]==chain)|(pdb_by_sdi["chnLettPrefix"]==""))]

	for domain in domains:
		#Split domain from PDB file
    	domain_file = os.path.join(pdb_path, "{}_{}_d{}.pdb".format(pdb, chain, domain["domNo"]))
    	with open(domain_file, "w") as domain_f:
        	subprocess.Popen(["/data/draizene/pdb-tools/pdb_rslice.py", 
        		"{}:{}".format(domain["frm"], domain["tu"]), pdb_file], 
        		stdout=delocc)

def add_sdi(cdd_sdi, cdd_ibis):
	cdd_sdi = pd.read_csv(cdd_sdi)
	cdd_ibis = pd.read_csv(cdd_ibis)

	all_sdi_domains = cdd_sdi.groupby(["pdbId", "chnLettPrefix"])

	domains = defaultdict(set)

	for row in cdd_ibis:
		try:
			sdi_domains = all_sdi_domains.get_group([row["pdb"], row["chain"]])
		except KeyError:
			continue

		resi = sorted(row["resi"].split(","), key=natural_keys)
		for sdi_domain in sdi_domains:
			if to_int(resi[0]) >= sdi_domain["frm"] and to_int(resi[-1]) <= sdi_domain["tu"]:
				break
		else:
			continue

		domain[(row["pdb"], row["chain"], sdi_domain["sdi"], sdi_domain["domNo"])].union(resi)

	for (pdb, chain, sdi, domNum), resi in domains.iteritems():
		resi = ",".join(sorted(resi, key=natural_keys))
		print "{}\t{}\t{}\t{}\t{}".format(pdb, chain, sdi, domNum, resi)





