===================
Command-line Options
===================

CATH-based workflow
-------------------

Prop3D was initially built to process the CATH Hierarchy. 

First we explain how to create a CATH dataset locally, using the cores present on your system. To build a dataset based on the entire CATH hierarchy, run the following command:

.. code-block:: bash
    
    export TEMP_RUN_NAME="create-cath-dataset" #Change to your desired run name
    export H5_PATH=/path/to/dataset.h5 #Path on HSDS endpoint to the new dataset

    python -m Prop3D.generate_data.main file:$TEMP_RUN_NAME \
        --retryCount 5 \
        --realTimeLogging \
        --hsds_file $H5_PATH

However, if you only want a subset of the dataset, add the argument ``--cathcode`` to specify any level of the hierarchy. It supports multiple entries (separated by space) and different levels such as ``2`` (for All Beta Class), ``2.60`` (for Sandwich Architectures), ``2.40.60`` (for Immunoglobulin-like Topologies), and ``2.40.60.10`` for (domains from the the Immunoglobulin Homologous Superfamily).

.. code-block:: bash
    
    export TEMP_RUN_NAME="create-cath-dataset" #Change to your desired run name
    export H5_PATH=/path/to/dataset.h5 #Path on HSDS endpoint to the new dataset

    python -m Prop3D.generate_data.main file:$TEMP_RUN_NAME \
        --retryCount 5 \
        --realTimeLogging \
        --maxLocalJobs 20 \
        --hsds_file $H5_PATH \
        --cathcode 2.60.40.10 1.10.10.10

You can also specify which cathcodes to exclude with `--skip_cathcode`.

.. code-block:: bash
    
    export TEMP_RUN_NAME="create-cath-dataset" #Change to your desired run name
    export H5_PATH=/path/to/dataset.h5 #Path on HSDS endpoint to the new dataset

    python -m Prop3D.generate_data.main file:$TEMP_RUN_NAME \
        --retryCount 5 \
        --realTimeLogging \
        --maxLocalJobs 20 \
        --hsds_file /projects/Prop3D/PDB.h5 \
        --cathcode 2.60.40.10 1.10.10.10

        
PDB-based workflow
------------------

Another option is build a dataset for the entire PDB, a list of PDB ids, or a local directory of PDB files by specifying the ``--pdb`` option. If no argument values are listed, the entire PDB will be used, otherwise pass in PDB codes or file names as input.

.. code-block:: bash
    
    export TEMP_RUN_NAME="create-cath-dataset" #Change to your desired run name
    export H5_PATH=/path/to/dataset.h5 #Path on HSDS endpoint to the new dataset

    python -m Prop3D.generate_data.main file:$TEMP_RUN_NAME \
        --retryCount 5 \
        --realTimeLogging \
        --hsds_file /projects/Prop3D/PDB.h5 \
        --pdb

Cloud-based Workflows via Toil
------------------------------

All of the previous commands use local parallelization, but because the workflow is based on Toil, it can also run on many cloud-based systems such as AWS, Google Cloud, Kubernetes, Slurm, etc. For more information on running on those environment, please see the `Toil Documentation <https://toil.readthedocs.io/en/latest/running/cloud/cloud.html>`_

Most academics still use Slurm, so we will demonstrate how to run Prop3D on Slurm. Below is an example script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --time=168:00:00
    #SBATCH -A YOUR_ACCOUNT
    #SBATCH -p YOUR_QUEUE
    #SBATCH --mem=64000
    #SBATCH -N 1
    #SBATCH -c 4
    #SBATCH --requeue

    #Load correct python

    mkdir -p ~/tmp3
    export TOIL_SLURM_ARGS="-c 2 -N 1 -A YOUR_ACCOUNT -t 10:00:00 -p YOUR_QUEUE"

    python -m Prop3D.generate_data.main file:create-Prop3D \
        --defaultCores 2 \
        --defaultMemory 32G \
        --retryCount 5 \
        --clean always \
        --realTimeLogging \
        --maxJobs=1000 \
        --workDir=~/tmp3 \
        --statePollingWait 120 \
        --disableAutoDeployment \
        --batchSystem slurm \
        --targetTime 1 \
        --hsds_file /projects/Prop3D/CATH.h5 \
        --restartable

Just save that to a file, such as `run_Prop3D.sh` and submit the job by calling `sbatch run_Prop3D.sh`


