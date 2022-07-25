# Prop3D-20 Dataset #

# Table of contents
1. [Description](#description)
2. [How to use](#use)
    1. [Download](#use-downlaod)
    2. [Use with HSDS](#use-hsds)
        1. [Load data](#use-hsds-load)
        2. [Analyze data with python (h5pyd)](#use-hsds-python)
        3. [Analyze single domain with Prop3D](#use-hsds-prop3d)
        4. [Train ML models with PyTorch Lightning using DeepUrfold](#use-hsds-deepurfold)
    3. [Use raw .h5 file](#use-h5)
3. [Dataset Data Description](#data-description)
    1. [Superfamilies present in dataset](#dataset-description-superfamies)
    2. [Dataset Organization (HDF Groups)](#data-description-organization)
    3. [Domain level Data (HDF Tables)](#data-description-tables)
        1. [Atom Tables](#data-description-tables-atom)
        2. [Residue Tables](#data-description-tables-residues)
        3. [Residue-Residue Contact Tables](#data-description-tables-edges)
    4. [Data Splits](#data-description-splits)
    5. [Representatives](#data-description-representatives)
4. [References](#references)

## 1. Description <a name="description"></a>

Prop3D-20 is a protein structure dataset that combines 3D atomic coordinates with biophysical and evolutionary properties for every atom in every "cleaned" domain structure from <b>20 CATH [1] Homologous Superfamilies.</b> Domain structures are 'cleaned' by adding missing residues with MODELLER [2], missing atoms with SCWRL4 [3], and protonating and energy minimizing (simple de-bump) with PDB2PQR [4]. We follow the CATH hierarchy in a hierarchical data format (HDF) file and include atomic level features, residue level features, residue-residue contact, and pre-calculated train (~80%) / test (~10%) / validation (~10%) splits for each superfamily derived from CATH's sequence identity clusters (e.g. S35 for 35% seq ID).

This dataset was originally stored in the Highly Scalable Data Service ([HSDS](http://www.github.com/HDFGroup/hsds)), and was exported into this raw .h5 file as backup. We recommend loading this data into HSDS for use in [h5pyd](http://www.github.com/HDFGroup/h5pyd), but the .h5 file can opened using [h5py](http://www.github.com/h5py/h5py) as well.

For more information on setting up HSDS and/or recreate this dataset, please see [http://www.github.com/bouralab/Prop3D/README.md](http://www.github.com/bouralab/Prop3D/README.md)

## 2. How to use <a name="use"></a>

### 1. First, download the dataset: <a name="use-download"></a>

```bash
export PROP3D_DATA=/home/$USER/Prop3D.h5 #Change to whereever yous tore Prop3D within HSDS
wget -O Prop3D-20.h5 6873024
```
### 2. If your are using HSDS, follow these steps: <a name="use-hsds"></a>
Follow the instructions at https://gitlab.com/uva-arc/hobo-request/-/blob/main/doc/single-node-k3s-hsds-install.md to configure HSDS, then run:

#### 1. Next, load the dataset into HSDS: <a name="use-hsds-load"></a>

```bash
hsload Prop3D-20.h5 $PROP3D_DATA
```

#### 2. Next, analyze the dataset with python (h5pyd): <a name="use-hsds-python"></a>

```python
import os
import h5pyd

with h5pyd.File(os.environ["PROP3d_DATA"], use_cache=False) as f:
     domain = f["2/30/30/100/domains/1kq2A00/atom"][:]

... #Code to process domain
```

#### 3. For more specific analysis, you can analyze single domains using Prop3D: <a name="use-hsds-prop3d"></a>

We recommend using [DistributedStructure](https://github.com/bouralab/Prop3D/blob/main/Prop3D/common/DistributedStructure.py) and/or [DistributedVoxelizedStructure](https://github.com/bouralab/Prop3D/blob/main/Prop3D/common/DistributedVoxelizedStructure.py) from the [Prop3D github repo](https://github.com/bouralab/Prop3D)

```python
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
```

#### 4. To train ML models in PyTorch Lightning, use the dataloader from DeepUrfold: <a name="use-hsds-deepurfold"></a>

To be written. Please see [DistributedDomainStructureDataModule](https://github.com/bouralab/DeepUrfold/blob/main/DeepUrfold/DataModules/DistributedDomainStructureDataModule.py) and [DistributedDomainStructureDataset](https://github.com/edraizen/DeepUrfold/blob/main/DeepUrfold/Datasets/DistributedDomainStructureDataset.py) from the [DeepUrfold github repo](https://github.com/bouralab/DeepUrfold) for more info.


### 3. If you want to use the raw .h5 locally with only h5py (no HSDS), run: <a name="use-h5"></a>

This is only for testing purposes. Not many of the Prop3D analyss functionality will work with raw .h5 files.

```bash
import os
import h5py

with h5py.File("/path/to/Prop3D-20.h5") as f: #Replace with the path the .h5 file present locallt on your machine, not inside HSDS
     domain = f["2/30/30/100/domains/1kq2A00/atom"] #Warning: might take a while!

... #Code to process domain
```

## 3. Dataset Data Description <a name="#data-description"></a>

### 1. Superfamilies present in dataset <a name="#dataset-description-superfamies"></a>
We include the following 20 superfamilies from CATH v4.3:

- Winged helix-like DNA-binding domain (1.10.10.10)
- EF-hand (1.10.238.10)
- Globins (1.10.490.10)
- Transferase (1.10.510.10)
- Ferritin (1.20.1260.10)
- SH3 (2.30.30.100)
- OB (2.40.50.140)
- Immunoglobulin (2.60.40.10)
- Beta-grasp (3.10.20.30)
- Ribosomal Protein S5 (3.30.230.10)
- KH (3.30.300.20)
- Sm-like (domain 2) (3.30.310.60)
- Gyrase A (3.30.1360.40)
- K Homology domain, type 1 (3.30.1370.10)
- Hedgehog domain (3.30.1380.10)
- P-loop NTPases (3.40.50.300)
- Rossmann-like (3.40.50.720)
- Ribonuclease Inhibitor (3.80.10.10)
- NTPase (3.90.79.10)
- Oxidoreductase (3.90.420.10)

### 2. Dataset Organization (HDF Groups) <a name="#data-description-organization"></a>

We organize each superfamily by their CATH hierarchy within the HDF file:

```
/1
    /10
        /10
            /10
                /domains
                /data_splits
                /representatives
        /238
            /10
                /domains
                /data_splits
                /representatives
        /490
            /10
                /domains
                /data_splits
                /representatives
        /510
            /10
                /domains
                /data_splits
                /representatives
    /20
        /1260
            /10
                /domains
                /data_splits
                /representatives
/2
    /30
        /30
            /100
                /domains
                /data_splits
                /representatives
    /40
        /50
            /140
                /domains
                /data_splits
                /representatives
    /60
        /40
            /10
                /domains
                /data_splits
                /representatives
/3
    /10
        /20
            /30
                /domains
                /data_splits
                /representatives
    /30
        /230
            /10
                /domains
                /data_splits
                /representatives
        /300
            /20
                /domains
                /data_splits
                /representatives
        /310
            /60
                /domains
                /data_splits
                /representatives
        /1360
            /40
                /domains
                /data_splits
                /representatives
        /1370
            /10
                /domains
                /data_splits
                /representatives
        /1380
            /10
                /domains
                /data_splits
                /representatives
    /40
        /50
            /300
                /domains
                /data_splits
                /representatives
            /720
                /domains
                /data_splits
                /representatives
    /80
        /10
            /10
                /domains
                /data_splits
                /representatives
    /90
        /79
            /10
                /domains
                /data_splits
                /representatives
        /420
            /10
                /domains
                /data_splits
                /representatives
```

### 3. Domain level Data (HDF Tables) <a name="#data-description-tables"></a>

Each "domain" group within a superfamily contains another group for each domain in that superfamily named by its CATH domain name e.g. "1kq2A00". This lower contains three datasets:

- 1. atom: atomic level properties for each atom in the domain, with the following columns: <a name="#data-description-tables"></a>

 ```python
 dtype([('serial_number', '<i8'), ('atom_name', 'S5'), ('residue_id', 'S8'), ('residue_name', 'S8'), ('chain', 'S2'), ('bfactor', '<f8'), ('X', '<f8'), ('Y', '<f8'), ('Z', '<f8'), ('H', '<f8'), ('HD', '<f8'), ('HS', '<f8'), ('C', '<f8'), ('A', '<f8'), ('N', '<f8'), ('NA', '<f8'), ('NS', '<f8'), ('OA', '<f8'), ('OS', '<f8'), ('F', '<f8'), ('MG', '<f8'), ('P', '<f8'), ('SA', '<f8'), ('S', '<f8'), ('CL', '<f8'), ('CA', '<f8'), ('MN', '<f8'), ('FE', '<f8'), ('ZN', '<f8'), ('BR', '<f8'), ('I', '<f8'), ('Unk_atom', '<f8'), ('C_elem', '<f8'), ('N_elem', '<f8'), ('O_elem', '<f8'), ('S_elem', '<f8'), ('H_elem', '<f8'), ('F_elem', '<f8'), ('MG_elem', '<f8'), ('P_elem', '<f8'), ('CL_elem', '<f8'), ('CA_elem', '<f8'), ('MN_elem', '<f8'), ('FE_elem', '<f8'), ('ZN_elem', '<f8'), ('BR_elem', '<f8'), ('I_elem', '<f8'), ('Unk_elem', '<f8'), ('vdw', '<f8'), ('charge', '<f8'), ('neg_charge', '<f8'), ('pos_charge', '<f8'), ('electrostatic_potential', '<f8'), ('is_electronegative', '<f8'), ('cx', '<f8'), ('is_concave', '<f8'), ('hydrophobicity', '<f8'), ('is_hydrophobic', '<f8'), ('biological_hydrophobicity', '<f8'), ('octanal_hydrophobicity', '<f8'), ('atom_asa', '<f8'), ('residue_rasa', '<f8'), ('residue_buried', '<f8'), ('ALA', '<f8'), ('CYS', '<f8'), ('ASP', '<f8'), ('GLU', '<f8'), ('PHE', '<f8'), ('GLY', '<f8'), ('HIS', '<f8'), ('ILE', '<f8'), ('LYS', '<f8'), ('LEU', '<f8'), ('MET', '<f8'), ('ASN', '<f8'), ('PRO', '<f8'), ('GLN', '<f8'), ('ARG', '<f8'), ('SER', '<f8'), ('THR', '<f8'), ('VAL', '<f8'), ('TRP', '<f8'), ('TYR', '<f8'), ('Unk_residue', '<f8'), ('phi', '<f8'), ('phi_sin', '<f8'), ('phi_cos', '<f8'), ('psi', '<f8'), ('psi_sin', '<f8'), ('psi_cos', '<f8'), ('is_helix', '<f8'), ('is_sheet', '<f8'), ('Unk_SS', '<f8'), ('is_regular_helix', '<f8'), ('is_beta_bridge', '<f8'), ('is_extended_strand', '<f8'), ('is_310_helix', '<f8'), ('is_pi_helix', '<f8'), ('is_hbond_turn', '<f8'), ('is_bend', '<f8'), ('no_ss', '<f8'), ('hydrophobic_atom', '<f8'), ('aromatic_atom', '<f8'), ('hbond_acceptor', '<f8'), ('hbond_donor', '<f8'), ('metal', '<f8'), ('eppic_entropy', '<f8'), ('is_conserved', '<f8'), ('vdw_radii', '<f8'), ('Unk_element', '<f8')])
 ```

- 2. residue: residue level properties for each residue in the domain, with the following columns: <a name="##data-description-tables-residues"></a>

```python
dtype([('residue_id', 'S8'), ('residue_name', 'S8'), ('chain', 'S2'), ('bfactor', '<f8'), ('X', '<f8'), ('Y', '<f8'), ('Z', '<f8'), ('vdw', '<f8'), ('charge', '<f8'), ('neg_charge', '<f8'), ('pos_charge', '<f8'), ('electrostatic_potential', '<f8'), ('is_electronegative', '<f8'), ('cx', '<f8'), ('is_concave', '<f8'), ('hydrophobicity', '<f8'), ('is_hydrophobic', '<f8'), ('biological_hydrophobicity', '<f8'), ('octanal_hydrophobicity', '<f8'), ('residue_rasa', '<f8'), ('residue_buried', '<f8'), ('ALA', '<f8'), ('CYS', '<f8'), ('ASP', '<f8'), ('GLU', '<f8'), ('PHE', '<f8'), ('GLY', '<f8'), ('HIS', '<f8'), ('ILE', '<f8'), ('LYS', '<f8'), ('LEU', '<f8'), ('MET', '<f8'), ('ASN', '<f8'), ('PRO', '<f8'), ('GLN', '<f8'), ('ARG', '<f8'), ('SER', '<f8'), ('THR', '<f8'), ('VAL', '<f8'), ('TRP', '<f8'), ('TYR', '<f8'), ('Unk_residue', '<f8'), ('phi', '<f8'), ('phi_sin', '<f8'), ('phi_cos', '<f8'), ('psi', '<f8'), ('psi_sin', '<f8'), ('psi_cos', '<f8'), ('is_helix', '<f8'), ('is_sheet', '<f8'), ('Unk_SS', '<f8'), ('is_regular_helix', '<f8'), ('is_beta_bridge', '<f8'), ('is_extended_strand', '<f8'), ('is_310_helix', '<f8'), ('is_pi_helix', '<f8'), ('is_hbond_turn', '<f8'), ('is_bend', '<f8'), ('no_ss', '<f8'), ('eppic_entropy', '<f8'), ('is_conserved', '<f8'), ('vdw_radii', '<f8'), ('Unk_element', '<f8')])
```

- 3. edges: residue-residue contact for each residue in the domain, with the following columns: <a name="#data-description-tables"></a>

```python
dtype([('src', 'S8'), ('dst', 'S8'), ('distance', '<f8'), ('angle', '<f8'), ('omega', '<f8'), ('theta', '<f8')])
```

### 4. Data Splits <a name="#data-description-splits"></a>
Next, each superfamily group (H level) also contains a 'data_split' group, which is split into 4 different sequence IDs to base splitting off of: S35, S60, S95, S100. Each of those sequence IDs is further split into train (~80%), train (~10%), and validation (~10%) groups, which map back to domains in the 'domains' group via soft links.

### 5. Representatives <a name="#data-description-representatives"></a>
Finally, each superfamily group (H level) also contains a 'representatives' group that contains representative domains for that superfamily computed by CATH, which map back to domains in the 'domains' group via soft links.

## 4. References

1. Sillitoe I, Bordin N, Dawson N, Waman VP, Ashford P, Scholes HM, Pang CSM, Woodridge L, Rauer C, Sen N, Abbasian M, Le Cornu S, Lam SD, Berka K, Varekova IH, Svobodova R, Lees J, Orengo CA. CATH: increased structural coverage of functional space. Nucleic Acids Res. 2021 Jan 8;49(D1):D266-D273. doi: 10.1093/nar/gkaa1079. PMID: 33237325; PMCID: PMC7778904.

2. Webb B, Sali A. Comparative Protein Structure Modeling Using MODELLER. Curr Protoc Bioinformatics. 2016 Jun 20;54:5.6.1-5.6.37. doi: 10.1002/cpbi.3. PMID: 27322406; PMCID: PMC5031415.

3. Krivov GG, Shapovalov MV, Dunbrack RL Jr. Improved prediction of protein side-chain conformations with SCWRL4. Proteins. 2009 Dec;77(4):778-95. doi: 10.1002/prot.22488. PMID: 19603484; PMCID: PMC2885146.

4. Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei GW, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Sci. 2018 Jan;27(1):112-128. doi: 10.1002/pro.3280. Epub 2017 Oct 24. PMID: 28836357; PMCID: PMC5734301.

