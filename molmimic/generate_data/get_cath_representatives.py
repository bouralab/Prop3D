import os
import sys
import pandas as pd
from molmimic.parsers.cath import CATHApi


def get_representatives(cathcode, data_dir):
    print("Running", cathcode)
    cath = CATHApi()
    hierarchy = cath.list_children_in_heirarchy(cathcode, 5)
    print(hierarchy)
    representatives = [child["example_domain_id"] for child in hierarchy["children"]]

    data_dir = os.path.join(data_dir, "train_files", *cathcode.split("."))
    all_domains = pd.read_hdf(os.path.join(data_dir, "DomainStructureDataset-full-train.h5"), "table")
    respresentative_domains = all_domains[all_domains["cathDomain"].isin(representatives)]

    respresentative_domains.to_hdf(os.path.join(data_dir, "DomainStructureDataset-representatives.h5"), "table")

if __name__ == "__main__":
    data_dir = sys.argv[1]
    for cathcode in sys.argv[2:]:
        get_representatives(cathcode, data_dir)
