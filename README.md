# molmimic

molmimic is a tool to annotate binding sites involved in protein-protein interactions (PPI) given the structure of only one of the interaction partners. This method differs from other binding site prediction programs becuase we treat the problem as a 3D image/volume segmentaion problem. We use Sparse 3D Fully Convolutionary Neural Networks, or Sparse 3D Unet, to make voxel level predictions, where each atom is mapped to one or more voxels.

In the future we hope to report binding partners, binding sites for protein-ligand interactions (PLI), sites shared by PPI an PLI for drug repurposing, and sites that are 'mimicking' known interactions, e.g. host-pathogen PPI.

# Requirements
Python requirements
```
pytorch 0.4
torchvision
pytorchviz
tnt
SparseConvNet
numpy
scikit-learn
Biopython
numba
seaborn
openbabel
freesasa
tqdm
pdb-tools (https://github.com/JoaoRodrigues/pdb-tools)
```
Non python requirements:
```
pdb2pqr
rosetta
snakemake
```
Optional:
```
DSSP
CX
flask
visdom
```
All requirements can be found in the Singularity Container on Singualrity Hub: 
```
https://www.singularity-hub.org/collections/818/
```
It can downloaded using: 
```bash
singularity pull shub://edraizen/SingularityTorch
```

#Data Generation
To create all of the necessary data, please run the Snakemake pipeline in the data directory

```
cd data
snakemake all
```

Feel free to configure a cluster to use with snakemake. Please see the README in the data directory for more information.

![Data Generation Pipeline](figures/data_generation_pipeline.png)

# Running
Training

Inference

# Citations
1. SparseConvNet: Graham, Engelcke, Maaten. [arXiv:1711.10275](https://arxiv.org/abs/1711.10275)
2. 3D U-net: Cicek, Abdulkadir, *et al.* International Conference on Medical Image Computing and Computer Assisted Intervention. 2016. [arXiv](https://arxiv.org/abs/1606.06650)
3. IBIS LUCA: Goncearenco, Shaytan, *et al.* Biophysical J. 2015. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/26213149)
4. IBIS: Shoemaker, Zhang, *et al.* Nucleic Acids Res. 2012. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/22102591)
5. DeepSite: Jiménez, Doerr, *et al.* Bioinformatics. 2017. doi:10.1093/bioinformatics/btx350 [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28575181)
6. 3DCNN: Torng, Altman. BMC Bioinformatics. 2017. doi:10.1186/s12859-017-1702-0. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28615003)
7. EnzyNet: Amidi, *et al.* [arXiv:1707.06017](https://arxiv.org/abs/1707.06017)

