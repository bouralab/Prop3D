import os
import pandas as pd

mmdb_path_prefix = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "MMDB"))

def calculate_mmdb_dataset():
	structural_domains_file = os.path.join(mmdb_path_prefix, "Structural_Domains.csv")
	if os.path.isfile(structural_domains_file):
		structural_domains = pd.read_csv(structural_domains_file)
	else:
		st_domain_file = os.path.join(mmdb_path_prefix, "StDomain.csv")
		st_domain_intvl_file = os.path.join(mmdb_path_prefix, "StDomainIntvl.csv")
		if os.path.isfile(st_domain_file) and os.path.isfile(st_domain_intvl_file):
			st_domain = pd.read_csv(st_domain_file)
			st_domain_intvl = pd.read_csv(st_domain_intvl_file)
			structural_domains = pd.merge(st_domain, st_domain_intvl, on=["sdi"])
			del st_domain
			del st_domain_intvl
		else:
			raise RuntimeError("You don't have the necessary files to build the MMDB dataset: StDomain.csv and StDomain.csv must be in the data directory")

	cdd_file = os.path.join(mmdb_path_prefix, "StructDomSfam.csv")
	cdds = pd.read_csv(cdd_file)

	cdd_groups = cdds.groupby("label")

	cdd_dir = os.path.join(mmdb_path_prefix, "cdd")
	if not os.path.exists(cdd_dir):
		os.makedirs(cdd_dir)

	for cdd, cdd_sdis in cdd_groups:
		cdd_sdis_positions = pd.merge(cdd_sdis, structural_domains, how="left", on="sdi")
		cdd_sdis_positions.to_csv(os.path.join(cdd_dir, "{}.csv".format(cdd)))
		del cdd_sdis_positions

if __name__ == "__main__":
	calculate_mmdb_dataset()


