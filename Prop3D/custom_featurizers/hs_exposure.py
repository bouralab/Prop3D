from typing import Union

import pandas as pd
from Bio.PDB import Structure, HSExposure
from Prop3D.custom_featurizers import CustomFeaturizer, StructureType

class HalfSphereExposure(CustomFeaturizer):
    def calculate_prop3D(self, path: str, structure: StructureType) -> tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
        hse = HSExposure()

        # Calculate HSEalpha
        exp_ca = hse.HSExposureCA(structure)

        residue_result = pd.DataFrame(
            [(res_id, exposure) for (chain_id, res_id), exposure in iter(exp_ca)],
            columns=["res_id", "hs_exposure"])
        residue_result = residue_result.set_index("res_id")

        residue_result = residue_result.assign(hs_exposure_norm=
            (residue_result.hs_exposure-residue_result.hs_exposure.mean())/residue_result.hs_exposure.std())

        #No need to calculate is_exposed or is_buried, since they are defined in YAML file with
        # 'from_feature', 'threshold', and 'equality'. But you can if you want.

        return None, residue_result