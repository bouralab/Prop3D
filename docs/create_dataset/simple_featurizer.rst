Create Simple Featurizers
=========================

For whatever reason, you might not even want to use toil following the CATH hierarchy. You can still create a Featurizer calculating all features for a single protein using the following example. You can also use this example in your own workflow system, e.g. nextflow or a simple loop:

.. code-block:: python

    from functools import partial
    from typing import Optional, Union
    
    import pandas as pd
    from Prop3D.common.featurizer import ProteinFeaturizer

    def calculate_features(pdb_file: str, name: str, level: Optional[str] = "atom", work_dir: Optional[str] = None, 
                           input_format: str = "pdb", update_features: Optional[list[str])],
                           to_file: Optional[Union[bool, str]] = None -> Union[str, pd.DataFrame]:
        """Calculate biophysical properties from a single protein

        Parameters
        ----------
        path : str
            Path to local structure file
        name : str
            Name of protein
        level : str
            'atom', 'residue', 'edge' level
        work_dir : str
            Where to save all temporary files from all run software. If None, use cwd.
        input_format : str
            Input format, "pdb", "mmCIF", what ever Bio.PDB understands 
        update_features : list
            List of features names or feature groups to update (while keeping the rest the same). Defualt is None, update all features

        Returns
        -------
        Either the path to the h5 file or a pandas dataframe
        """
        assert level in ['atom', 'residue', 'edges']

        structure = ProteinFeaturizer(
            pdb_file, name, None, work_dir,
            force_feature_calculation=update_features is None,
            update_features=update_features)

        feature_calculator = {
            "atom": structure.calculate_flat_features,
            "residue": structure.calculate_flat_residue_features,
            "edges": partial(structure.calculate_graph, edgelist=True)
        }
        
        df, _ = feature_calculator[level](write=False)

        if level=="edges":
            df["src"] = df["src"].apply(lambda s: "".join(map(str,s[1:])).strip())
            df["dst"] = df["dst"].apply(lambda s: "".join(map(str,s[1:])).strip())
        else:
            df = structure.get_pdb_dataframe(include_features=True, coarse_grained = level=="residue")
        
        if (isinstance(to_file, bool) and to_file) or isinstance(to_file, str):
            if isinstance(to_file, bool):
                to_file = f"{name}.h5"

            df.to_hdf(to_file)

            return to_file 
        
        return df