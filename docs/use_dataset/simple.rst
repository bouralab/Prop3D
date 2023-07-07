Simple 
======

For all of these examples, we will use the Prop3D-20 dataset:

.. code-block:: bash

    export HS_USER=protein
    export HS_PASSWORD=protein
    export HS_ENDPOINT=http://hsds.pods.virginia.edu
    export PROP3D_DATA=/projects/Prop3D/Prop3D-20.h5

Get a single protein from the dataset with python (h5pyd)
---------------------------------------------------------

.. code-block:: python

    import os
    import h5pyd

    with h5pyd.File(os.environ["PROP3D_DATA"], use_cache=False) as f:
        domain = f["2/30/30/100/domains/1kq2A00/atom"][:]

    ... #Code to process domain

Analyze a single domain using Prop3D
------------------------------------

We recommend using `DistributedStructure <https://github.com/bouralab/Prop3D/blob/main/Prop3D/common/DistributedStructure.py>` and/or `DistributedVoxelizedStructure <https://github.com/bouralab/Prop3D/blob/main/Prop3D/common/DistributedVoxelizedStructure.py>` from the `Prop3D GitHub repo <https://github.com/bouralab/Prop3D>``

.. code-block:: python

    import os
    import numpy as np
    from Prop3D.common.DistributedStructure import DistributedStructure
    from Prop3D.common.DistributedVoxelizedStructure import DistributedVoxelizedStructure

    #Load structure and perform actions in atomic coordinate space
    structure = DistributedStructure(
        os.environ["PROP3D_DATA"], 
        key="2/30/30/100", 
        cath_domain_dataset="1kq2A00")

    #Save 5 random rotation sampling from the SO(3) group to pdb files
    for i, (r, M) in enumerate(structure.rotate(num=5)):
        structure.save_pdb(f"1kq2A00_rotation{i}.pdb")

    #Load structure and perform actions in voxelized coordinate space
    structure = DistributedVoxelizedStructure(
        os.environ["PROP3D_DATA"], 
        key="2/30/30/100", 
        cath_domain_dataset="1kq2A00")

    #Save 5 random rotation sampling from the SO(3) group to numpy files
    for i, (r, M) in enumerate(structure.rotate(num=5)):
        coords, feats = structure.map_atoms_to_voxel_space(autoencoder=True)
        np.savez(f"1kq2A00_rotation{i}.npz", coords=coords, feats=feats)
