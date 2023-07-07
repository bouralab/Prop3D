"""
Meadowlark
==========
A collection on scripts to process indidual protein structures for use in machine learning tasks. Proteins can be:

1) 'Cleaned' by adding missing residues and atoms;
2) Featurized with atom- and residue-based biophysical prooperties calculated using known structural bioinformatics tool that have been Dockerized (see `Prop3D.ml`). 
3) Convert proteins along with there features into sparse 3D volumes for use in Sparse 3DCNNs
"""

__all__ = ['AbstractStructure', 'DistributedStructure', 'DistributedVoxelizedStructure', 'features', 'featurizer', 'LocalStructure', 'ProteinTables']