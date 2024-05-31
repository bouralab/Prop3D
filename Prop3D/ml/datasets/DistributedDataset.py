import os
import time
from numbers import Number

import numpy as np
import pandas as pd

try:
    from torch.utils.data import Dataset
    from torch.utils.data import ConcatDataset as _ConcatDataset
except ImportError:
    raise ImportError("In order to the Prop3D datasets, you must install pytorch")

import h5pyd

class ConcatDataset(_ConcatDataset):
    def __getattr__(self, name):
        if hasattr(self.datasets[0], name):
            if callable(getattr(self.datasets[0], name)):
                def method(*args, **kwds):
                    r = [getattr(d, name)(*args, **kwds) for d in self.datasets]
                    return r
                return method
            else:
                return [getattr(d, name) for d in self.datasets]
        else:
            raise AttributeError(name)

class DistributedDataset(Dataset):
    """Read dataset from h5 file. If key specifies a dataset, each row is an
    independent sample. If kay specifies a group, each dataset is an independent
    sample.

    Paramaters
    ----------
    path : str
        The name of the full h5 file with all groups and datasets
    key : str
        The key to the dataset or group, specifying all intermediate groups
    test : bool
        Set mode to testing. This saves all rows or datasets IDs as an embedding
        to compare against
    dataset_group_name : str
        Name of group to split datasets and data_splits. Default is 'datasets'
    file_mode : str
        Open h5 file for reading and or writing. Writing should only be used if
        creating data_splits
    """

    def __init__(self, path, key, test=False, dataset_group_name="datasets", use_keys=None, ignore_keys=None, file_mode="r", use_cache=False, retries=100, label_encoder_classes=None):
        self.path = path
        self.key = key
        self.test = test
        self.file_mode = file_mode
        self.use_cache = use_cache
        self.retries = retries
        self._f = None

        try:
            f = h5pyd.File(self.path, self.file_mode, use_cache=self.use_cache, retries=self.retries)
        except Exception as e:
            import traceback as tb
            with open(f"{key.replace('/', '_')}.error", "w") as f:
                print(tb.format_exc(), file=f)
            raise

        data = f[key]

        self.embedding = None
        if not isinstance(data, h5pyd.Dataset):
            self.is_dataset = False
            if dataset_group_name in data.keys():
                self.key = os.path.join(key, dataset_group_name)
                data = data[dataset_group_name]
            self.order = sorted(data.keys())


            if ignore_keys is None:
                ignore_keys = []
            elif not isinstance(ignore_keys, (list, tuple)):
                ignore_keys = [ignore_keys]

            if use_keys is not None:
                #Preserve order of use_keys and allow multiple of same key to support oversampling
                if not isinstance(use_keys, (list, tuple)):
                    use_keys = [use_keys]
                #self.order = sorted(list(set(self.order).intersection(set(use_keys))))
                self.order = [k for k in use_keys if k in self.order]

            self.order = [os.path.join(self.key, k) for k in self.order if k not in ignore_keys]

            if self.test:
                if label_encoder_classes is not None:
                    self.embedding = self.encode_labels(label_encoder_classes, from_fitted=True)
                else:
                    self.embedding = self.encode_labels(self.order)
        else:
            self.order = list(range(len(data)))
            self.is_dataset = True

    @classmethod
    def create_multiple(cls, path, keys, balance=None, **kwds):
        if not isinstance(keys, (list, tuple)):
            keys = [keys]

        datasets = [cls(path, key, **kwds) for key in keys]

        if len(datasets) == 1:
            return datasets[0]

        print("Created", len(datasets), "datates")

        if balance in ["oversample", "undersample"]:
            if balance == "oversample":
                from imblearn.over_sampling import RandomOverSampler as Sampler
            else:
                from imblearn.under_sampling import RandomUnderSampler as Sampler

            print(f"Balancing dataset with {len(datasets)} classes with each with "
                + f"{[len(d.order) for d in datasets]} samples")

            X, y = zip(*[(d, name) for name, dataset in zip(keys, datasets) \
                for d in dataset.order])

            sampler = Sampler(random_state=0)
            X_resampled, y_resampled = sampler.fit_resample(np.array(X).reshape(-1, 1), y)
            X = pd.DataFrame({"X":X_resampled.flatten(),"y":y_resampled})
            X.to_csv(f"resampled_dataset_{balance}_{np.random.randint(4)}.csv")
            datasets = [cls(path, key, use_keys=df["X"], **kwds) for key, df in X.groupby("y")]

        return ConcatDataset(datasets)

    @staticmethod
    def encode_labels(labels, from_fitted=False):
        from sklearn.preprocessing import LabelEncoder
        if not from_fitted:
            embedding = LabelEncoder().fit(labels)
        else:
            embedding = LabelEncoder()
            embedding.classes_ = labels
        return embedding

    def __len__(self):
        return len(self.order)

    def __getitem__(self, index):
        for i in range(self.retries):
            try:
                if not self.is_dataset:
                    return self.f[self.order[index]]
                else:
                    return self.f[self.key][self.order[index]]
            except Exception:
                pass
            
            time.sleep(1) #Issue with multiple requests? wait a second

    @property
    def f(self):
        if self._f is None:
            self._f = h5pyd.File(self.path, self.file_mode, use_cache=self.use_cache, retries=self.retries)
        return self._f

    @f.deleter
    def f(self):
        if self._f is not None:
            self._f.close()


class DistributedDatasetSplit(DistributedDataset):
    def __init__(self, path, key, split_key):
        f = h5pyd(path, 'r')
        self.group = f[key]
        assert hasattr(self.group, "data_split")
        self.data = self.group.data_split[split_key]
        self.order = sorted(self.data.keys())

    def split(self, split_size={"train":0.8, "validation":0.1, "test":0.1}, ):
        if self.file_mode == "r":
            if isinstance(self.data, h5pyd.Dataset):
                return {self.data.attrs["data_split"][split_name] \
                    for split_name in split_size.items()}
            else:
                name = self.data.visit(
                    lambda n: n if f'data_splits/{split_name}' in n else None)
                return f[name]
        else:
            raise NotImplementedError("You must subclass DistributedDataset in order to split your data")
