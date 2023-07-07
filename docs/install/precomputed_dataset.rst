Install Pre-computed Dataset
----------------------------

**Note: If you only want to use the dataset, this is unnecessary. See instructions on how to use the `Pre-computed Dataset <precomputed_datasets>`

*Note: Make sure you install HSDS before running this step*

First, download the data set as a raw .h5 file from Zenodo:

.. code-block:: bash

	export PROP3D_DATA=/home/$USER/Prop3D.h5 #Change to wherever you want to store the Prop3D within HSDS
	wget https://zenodo.org/record/6873024/files/Prop3D-20.h5

Then, import the .h5 file into the HSDS endpoint:

.. code-block:: bash

    python -m Prop3D.generate_data.hsds file:load-Prop3D \
        --maxLocalJobs 96 
        --load Prop3D-20.h5 $PROP3D_DATA #Change last path to whatever you want to name the file in HSDS


Alternatively, you can run the slower non-parallelized version:

.. code-block:: bash
    
    hsload Prop3D.h5 /home/$USER/Prop3D.h5 #Change last path to whatever you want to name the file in HSDS
