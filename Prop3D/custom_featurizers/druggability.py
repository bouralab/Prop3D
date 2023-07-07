import shutil
from typing import Union
from pathlib import Path

import pandas as pd
from Bio import PDB
from Prop3D.parsers.container import Container
from Prop3D.custom_featurizers import CustomFeaturizer, StructureType

class Druggability(CustomFeaturizer, Container): #Sublcass Container to have the ablity to run docker containers without leaving python
    IMAGE = "docker://fpocket/fpocket"
    ENTRYPOINT = "fpocket"
    PARAMETERS = [
        "-d", #Always output condesned form, easier to parse
        [("pdb_file", "str", ["-f", "{}"]) ]#parameter name (used as var name in python), parameter value type, formatting options
    ]
    CONTAINER_FILE_PREFIX = "/WORKDIR" #Set the working direcotry in the container


    def calculate_prop3D(self, path: str, structure: StructureType, clean: bool = True) -> tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
        """Run fpocket for input structure and save the druggability score of the pocket to each atom in the pocket
        """
        #Call FPocket
        output = next(self(pdb_file=path))

        results_dir = Path(self.work_dir) / f"{Path(path.stem)}_out"
        if not results_dir.is_dir():
            raise RuntimeError("Error running Fpocket")

        #Pocket info, dataframe with cols: cav_id drug_score volume nb_asph inter_chain apol_asph_proportion mean_asph_radius as_density mean_asph_solv_acc mean_loc_hyd_dens flex hydrophobicity_score volume_score charge_score polarity_score a0_apol a0_pol af_apol af_pol n_abpa ala cys asp glu phe gly his ile lys leu met asn pro gln arg ser thr val trp tyr chain_1_type chain_2_type num_res_chain_1 num_res_chain_2 lig_het_tag name_chain_1 name_chain_2
        pocket_info = pd.read_csv(output, sep=" ")[["cav_id", "drug_score"]].set_index("cav_id")

        atom_results = pd.DataFrame(0., index=[a.serial for a in structure.get_atoms], columns=["druggability_score"])
        
        for pocket_file in results_dir.glob("*_atm.pdb"):
            cav_id = pocket_file.stem.split("_").replace("pocket", "")
            drug_score = pocket_info.drug_score[cav_id]
            with pocket_file.open() as f:
                for line in f:
                    if line.startswith("ATOM"):
                        serial_number = int(line[6:11])
                        if atom_results.loc[serial_number].druggability_score < drug_score:
                            atom_results.loc[serial_number, "druggability_score"] = drug_score

        if clean:
            shutil.rmtree(results_dir)

        return atom_results, None