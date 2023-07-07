Create Custom Featurizers
=========================

Prop3D includes the ability to add your own software into the featurization pipeline. All you need is a docker image (if using custom software), a ``Prop3D.meadowlark.CustomFeaturizer`` to run custom software and parse results, and updated yaml file with group and feature definitions.

When creating your own Docker images, please try to follow the `UCSC CGL Docker Lib Philosophy <https://github.com/BD2KGenomics/cgl-docker-lib>`_. However, if the docker image is already built or the ``wrapper.sh`` endpoint does not work, just add an attribute your ``CustomFeaturizer`` subclass specifying the entry point of software inside the Docker image.

.. code-block:: python

    from Prop3D.parsers.contianer import Container
    from Prop3D.custom_featurizers import CustomFeaturizer, StructureType

    class NewFeaturizer(CustomFeaturizer, Container):
        ENTRYPOINT = /path/to/software #Only add if you the UCSC CGL theory didn't work or you don't want to use it
        IMAGE = "docker://USER/image" #Specify the path to docker image
        PARAMETERS = ([
            [("parameter_name", "type", ["--format_option", "{}"]) #parameter name (used as var name in python), parameter value type, formatting options
        ]

        def calculate_prop3D(self, path: str, structure: StructureType) -> tuple[Union[pd.DataFrame, None], Union[pd.DataFrame, None]]:
            ... code to calculate your features ...
            return atom_features, residue_features

Summary
-------
#. Create a Docker image for the new software of interest, uploaded to Docker cloud
#. Subclass ``Prop3D.custom_featurizers.CustomFeaturizer`` and ``Prop3D.parsers.container.Container`` to run the docker image and parse the results. This subclass must contain a ``calculate_prop3D`` method that returns a pandas DataFrame.
#. Modify the custom features YAML file (``Prop3D/custom_featurizers/custom_features.yaml``) by creating a new top-level group it will be in, and then each feature with the specified argument including the new Container you just subclassed. If it exists, make sure it is not hidden with "." in front of its name.

Example with BioPython (Half Sphere Exposure)
---------------------------------------------

First, create a new CustomFeaturizer:

.. literalinclude:: ../../Prop3D/custom_featurizers/hs_exposure.py
  :language: python    

.. code-block:: yaml

    - hs_exposure: # Unique Group Name
        - name: "hs_exposure" # Unique Feature Name
          default: 0.
          aggregate: mean # How to combine for multiple values, e.g. for overlapping atoms in a voxel or if atom feature, how to combine to create a residue level feature
          residue: true #True if calculated at the residue level
          bool: false #True if value is Boolean
          min: -180
          max: 180
          parser: Prop3D.custom_featurizers.HalfSphereExposure
        
        - name: "hs_exposure_norm"
          default: 0.
          aggregate: mean
          residue: true
          bool: false
          min: 0
          max: 1
          parser: Prop3D.custom_featurizers.HalfSphereExposure

        - name: 'is_exposed' #Create new Boolean value based on thresholding previous values
          default: 0
          aggregate: max
          residue: true
          bool: true
          min: 0.
          max: 1.
          threshold: 0.2
          from_feature: "hs_exposure_norm"
          equality: ">="
          parser: Prop3D.custom_featurizers.HalfSphereExposure
        
        - name: 'is_buried'
          default: 0
          aggregate: max
          residue: true
          bool: true
          min: 0.
          max: 1.
          threshold: 0.2
          from_feature: "hs_exposure_norm"
          equality: "<"
          parser: Prop3D.custom_featurizers.HalfSphereExposure

Example with custom software (Fpocket)
--------------------------------------

In this example we will use Fpocket to calculate druggability scores for every atom. A docker image already exists for Fpocket.

First, create a new CustomFeaturizer Container to run the fpocket docker image and parse the results:

.. literalinclude:: ../../Prop3D/custom_featurizers/druggability.py
  :language: python  

Next, we will append the new feature group along with its associated features to the custom features yaml:

.. code-block:: yaml

    - druggability: # Unique Group Name
        - name: "druggability_score" # Unique Feature Name
          default: 0.
          aggregate: mean # How to combine for multiple values, e.g. for overlapping atoms in a voxel or if atom feature, how to combine to create a residue level feature
          residue: false #True if calculated at the residue level
          bool: false #True if value is Boolean
          min: 0.
          max: 1.
          parser: Prop3D.custom_featurizers.Druggability

        - name: 'is_druggable' #Create new Boolean value based on thresholding previous values
          default: 0
          aggregate: max
          residue: true
          bool: true
          min: 0.
          max: 1.
          threshold: 0.5
          from_feature: "druggability_score"
          equality: ">="
          parser: Prop3D.custom_featurizers.Druggability
        
        - name: 'not_druggable'
          default: 0
          aggregate: max
          residue: true
          bool: true
          min: 0.
          max: 1.
          threshold: 0.5
          from_feature: "druggability_score"
          equality: "<"
          parser: Prop3D.custom_featurizers.Druggability