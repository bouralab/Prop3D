import os

try:
    import torch
except ImportError:
    raise ImportError("In order to the Prop3D datasets, you must install pytorch")

import numpy as np
from scipy.stats import special_ortho_group

from Prop3D.ml.datasets.DistributedDataset import DistributedDataset
from Prop3D.common.DistributedVoxelizedStructure import DistributedVoxelizedStructure

class DistributedDomainStructureDataset(DistributedDataset):
    hierarchy = ["C", "A", "T", "H", "S35", "S60", "S95", "S100"]
    rvs = np.eye(3)

    def __init__(self, path, key, use_features=None, predict_features=None,
      truth_key=None, cluster_level="S35", test=False, validation=False, 
      representatives=False, domains=None, all_domains=False, file_mode="r", 
      dataset_group_name=None, use_keys=None, ignore_keys=None, 
      return_structure=False, label_encoder_classes=None):
        assert [validation, test].count(True)<2, "Can only select none or one at a time"
        self.use_features = use_features
        self.cluster_level = cluster_level
        self.use_features = use_features
        self.predict_features = predict_features
        self.truth_key = truth_key
        self.test = test
        self.validation = validation
        self.domains = domains
        self.all_domains = all_domains
        self.representatives = representatives
        self.return_structure = return_structure

        if truth_key is not None:
            if isinstance(predict_features, (list, tuple)) and len(predict_features)>0:
                print("ignoring predict_features")
        elif isinstance(predict_features, (list, tuple)) and len(predict_features)>0:
            assert isinstance(use_features, (list, tuple)) and len(use_features) > 0, \
                "Cannot train on all features and predict a feature in that set"
            assert len(set(predict_features).intersection(set(use_features)))==0, \
                "Cannot train on and predict the same features"

        assert [isinstance(domains,(str,list,tuple)), all_domains, representatives].count(True)<2

        if dataset_group_name is None:
            if all_domains or isinstance(domains, (list,tuple)):
                dataset_group_name = "domains"
            elif representatives:
                dataset_group_name = "representatives"
            elif cluster_level is not None:
                dataset_group_name = f"data_splits/{cluster_level}/"
                if self.validation:
                    dataset_group_name += "validation"
                elif self.test:
                    dataset_group_name += "test"
                else:
                    dataset_group_name += "train"
            else:
                dataset_group_name = "domains"
        self.dataset_group_name = dataset_group_name

        if key is None:
            key = "/"
            self.dataset_group_name = ""

        super().__init__(path, key, test=self.test, dataset_group_name=dataset_group_name,
            use_keys=self.domains, file_mode=file_mode, label_encoder_classes=label_encoder_classes)


    @staticmethod
    def encode_labels(labels, from_fitted=False):
        from sklearn.preprocessing import LabelEncoder
        if not from_fitted:
            embedding = LabelEncoder().fit([l.split("/")[-1] for l in labels])
        else:
            embedding = LabelEncoder()
            embedding.classes_ = labels
        return embedding

    def __getitem__(self, index):
        #xyz[resn][resa][atom] = np.array([x,y,z])
        # all_atoms = np.array(t['xyz'][res,])[0,] #[L, 14, 3]
        # coords_dict_chain['N_chain_'+letter]=all_atoms[:,0,:].tolist()
        # coords_dict_chain['CA_chain_'+letter]=all_atoms[:,1,:].tolist()
        # coords_dict_chain['C_chain_'+letter]=all_atoms[:,2,:].tolist()
        # coords_dict_chain['O_chain_'+letter]=all_atoms[:,3,:].tolist()
        # my_dict['coords_chain_'+letter]=coords_dict_chain

        voxelizer = DistributedVoxelizedStructure(
            self.path, key, cath_domain_dataset, 
            use_features=self.use_features, predict_features=self.predict_features,
            replace_na=True)
        
        coords_dict_chain = {
            "N_chain_{voxelizer.chain}}": [],
            "CA_chain_{voxelizer.chain}}": [],
            "C_chain_{voxelizer.chain}}": [],
            "O_chain_{voxelizer.chain}}": [],
        }
        feats = []
        for r in voxelizer.get_residues():
            for atom_type in ("N", "CA", "C", "O"):
                for a in r:
                    if a[atom_type]:
                        coords_dict_chain["{atom_type}_chain_{voxelizer.chain}}"].append(a[["X", "Y", "Z"]].tolist())
                        feats.append(a[self.predict_features].tolist())
                    else:
                        coords_dict_chain["{atom_type}_chain_{voxelizer.chain}}"].append([np.nan]*3)
                        feats.append([np.nan]*len(self.predict_features))
        
        seq = voxelizer.get_sequence()
        return {
            "coords_chain_{voxelizer.chain}": coords_dict_chain,
            "masked_list": np.zeros(len(seq)),
            "visible_list": np.ones(len(seq)),
            "num_of_chains": 1,
            "seq": seq,
            "prop3d_features": feats
        }