==================
Basic Installation
==================

Before you begin, we recommend creating a new environment for Prop3D:

.. code-block:: bash

    conda create -n Prop3D_env python=3.10
    conda activate Prop3D_env

**WARNING: Prop3D requires Toil and PyTorch, which must be installed prior to Prop3D to allow for you to install with special options**

Prerequisites
-------------

#. **Install Singularity (preferred) or Docker**
Follow these steps:** `<https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps>`_

#. **Install Toil (Allows Building Datasets)**

If you plan on building your own datasets, you must you must install `Toil <https://github.com/DataBiosphere/toil>`_ first, choosing your preferred Toil extras to install, such as cloud-based cluster choice. 

.. code-block:: bash

    pip install toil #Add specific workflow manager you will be use, e.g. toil[aws]

#. **Install PyTorch (Allows Using Datasets)**

If you plan on using the precomputed dataset or one of your own own datasets, it is recommended to install PyTorch. Please follow their install `instructions <https://pytorch.org/get-started/locally/>`_.

An example way to install may be:

.. code-block:: bash

    #Using conda with GPU support:
    conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia

    #Using conda (no GPU support):
    conda install pytorch torchvision torchaudio cpuonly -c pytorch

    #Using pip with GPU support:
    pip3 install torch torchvision torchaudio

    #Using pip (no GPU support):
    pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

Install Prop3D
--------------

Make sure all of the prereuqisites are installed, then run:

.. code-block:: bash

    pip install prop3D

Next, configure HSDS to run with the pre-calculated datasets:

please set the following variables by running ``hsconfig``:

.. code-block::

    # HDFCloud configuration file
    hs_endpoint = http://prop3d-hsds.pods.uvarc.io
    hs_username = None
    hs_password = None
    hs_api_key = None

These values can also be set though the environmental variables, ``HS_ENDPOINT``, ``HS_USERNAME`` and ``HS_PASSWORD``.

The path to each pre-computed dataset on our HSDS endpoint is:

*Prop3D-20* : ``/projects/Prop3D/Prop3D-20.h5``

In Progress:

*PDB* : ``/projects/Prop3D/PDB.h5``

Extras
------

If you are using MinIO underlying HSDS, you should install Prop3D with the ``s3`` extra:

.. code-block:: bash

    pip install prop3D[s3]


