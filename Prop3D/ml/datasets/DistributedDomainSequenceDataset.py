import os

try:
    import torch
except ImportError:
    raise ImportError("In order to the Prop3D datasets, you must install pytorch")

import numpy as np
from scipy.stats import special_ortho_group

from Prop3D.ml.datasets.DistributedDataset import DistributedDataset
from Prop3D.common.DistributedVoxelizedStructure import DistributedVoxelizedStructure

class DistributedDomainSequenceDataset(DistributedDataset):
    hierarchy = ["C", "A", "T", "H", "S35", "S60", "S95", "S100"]
    rvs = np.eye(3)

    def __init__(self, path, key, use_features=None, predict_features=None,
      truth_key=None, cluster_level="S35", nClasses=1,
      test=False, validation=False, representatives=False, domains=None, all_domains=False,
      file_mode="r", dataset_group_name=None, use_keys=None, ignore_keys=None, 
      return_structure=False, label_encoder_classes=None):
        assert [validation, test].count(True)<2, "Can only select none or one at a time"
        self.use_features = use_features
        self.cluster_level = cluster_level
        self.use_features = use_features
        self.predict_features = predict_features
        self.truth_key = truth_key
        self.nClasses = nClasses
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
            assert (use_features is None) or (isinstance(use_features, (list, tuple)) and len(use_features) > 0), \
                "Cannot train on all features and predict a feature in that set"
            assert (use_features is None) or (isinstance(use_features, (list, tuple)) and len(set(predict_features).intersection(set(use_features)))==0), \
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
        cath_domain_dataset = super().__getitem__(index)
        key = self.order[index]

        voxelizer = DistributedVoxelizedStructure(
            self.path, key, cath_domain_dataset, rotate=False,
            use_features=self.use_features, predict_features=self.predict_features,
            replace_na=True, coarse_grained=True)

        if self.return_structure:
            return voxelizer

        sequence = voxelizer.get_sequence()
        predict_features = voxelizer.features[self.predict_features].tolist()

        if self.use_features is not None:
            use_features = voxelizer.features[self.use_features]
            return (sequence, use_features), predict_features
        else:
            return sequence, predict_features


