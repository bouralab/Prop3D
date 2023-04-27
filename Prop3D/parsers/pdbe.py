import pandas as pd
from Prop3D.parsers.json import JSONApi
from Prop3D.generate_data.data_stores import data_stores

class PDBEApi(JSONApi):
    def __init__(self, store=None, work_dir=None, download=True, max_attempts=2, job=None):
        url = "https://www.ebi.ac.uk/pdbe/api/"
        if store is None:
            assert job is not None
            store = data_stores(job).pdbe_store if store is None else store
        super(PDBEApi, self).__init__(url, store, work_dir=work_dir,
            download=download, max_attempts=max_attempts, job=job)

    def get_pdb_residues(self, pdb, chain):
        residues = self.get("pdb/entry/residue_listing/{}/chain/{}".format(pdb.lower(), chain))

        if isinstance(residues, dict) and pdb not in residues:
            print(residues)

        residue_map = {"residueNumber":[], "pdbResidueNumber":[]}
        for molecule in residues[pdb]["molecules"]:
            for chain in molecule["chains"]:
                for residue in chain["residues"]:
                    resNum = residue["residue_number"]
                    pdbResNum = "{author_residue_number}{author_insertion_code}".format(**residue)
                    residue_map["residueNumber"].append(resNum)
                    residue_map["pdbResidueNumber"].append(pdbResNum)

        return pd.DataFrame(residue_map)

    def get_uniprot_mapping(self, accession):
        mappings = self.get("mappings/all_isoforms/{}".format(accession))
        return mappings
