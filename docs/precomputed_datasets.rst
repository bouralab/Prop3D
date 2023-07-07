Pre-computed Datasets
=====================

To use our pre-computed datasets, please set the following variables with ``hsconfig``:

.. code-block::

    # HDFCloud configuration file
    hs_endpoint = http://hsds.pods.uvarc.io
    hs_username = protein
    hs_password = protein
    hs_api_key = None

These values can also be set though the environmental variables, ``HS_ENDPOINT``, ``HS_USERNAME`` and ``HS_PASSWORD``.

The path to each pre-computed dataset on our HSDS endpoint is:

*Prop3D-20* : ``/projects/Prop3D/Prop3D-20.h5``

*PDB* : ``/projects/Prop3D/PDB.h5``