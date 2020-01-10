import pandas as pd
from molmimic.parsers.json import JSONApi

class PDBEApi(JSONApi):
    def __init__(self, store, work_dir=None, download=True, max_attempts=2):
        url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/"
        super(PDBEApi, self).__init__(url, store, work_dir=work_dir,
            download=download, max_attempts=max_attempts)

    def fix_key(self, key):
        return key.lower()

    def get_pdb_residues(self, pdb, chain):
        residues = self.get("residue_listing/{}/chain/{}".format(pdb, chain))

        residue_map = {"residueNumber":[], "pdbResidueNumber":[]}
        for molecule in residues[pdb]["molecules"]:
            for chain in molecule["chains"]:
                for residue in chain["residues"]:
                    resNum = residue["residue_number"]
                    pdbResNum = "{author_residue_number}{author_insertion_code}".format(**residue)
                    residue_map["residueNumber"].append(resNum)
                    residue_map["pdbResidueNumber"].append(pdbResNum)

        return pd.DataFrame(residue_map)
