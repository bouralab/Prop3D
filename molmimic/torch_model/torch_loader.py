import os
import traceback
from itertools import groupby

try:
    import torch
    from torch.utils.data.dataset import Dataset
    from torch.utils.data import DataLoader
except ImportError:
    torch = None
    Dataset = object
    DataLoader = None
    
import pandas as pd
import numpy as np

try:
    import sparseconvnet as scn
except:
    #Allow module to load without scn to read data
    scn = None

from molmimic.biopdbtools import Structure, InvalidPDB

def dense_collate(data, batch_size=1):
    batch = {"data":None, "truth":None}

    for i, d in enumerate(data):
        if batch["data"] is None:
            if len(d["data"].shape) == 3:
                batch["data"] = torch.FloatTensor(batch_size, 1, d["data"].shape[0], d["data"].shape[1], d["data"].shape[2])
                batch["truth"] = torch.FloatTensor(batch_size, 1, d["truth"].shape[0], d["truth"].shape[1], d["truth"].shape[2])
            else:
                batch["data"] = torch.FloatTensor(batch_size, d["data"].shape[3], d["data"].shape[0], d["data"].shape[1], d["data"].shape[2])
                batch["truth"] = torch.FloatTensor(batch_size, d["truth"].shape[3], d["truth"].shape[0], d["truth"].shape[1], d["truth"].shape[2])

        batch["data"][i, ...] = torch.from_numpy(d["data"])
        batch["truth"][i, ...] = torch.from_numpy(d["truth"])

    return batch

def sparse_collate(data, input_shape=(256,256,256), create_tensor=False):
    if scn is None or not create_tensor:
        batch = {
            "indices": [],
            "data": [],
            "truth": []
            }

        def add_sample(indices, features, truth):
            batch["indices"].append(indices)
            batch["data"].append(features)
            batch["truth"].append(truth)

    else:
        inputSpatialSize = torch.LongTensor(input_shape)
        batch = {
            "data": scn.InputBatch(3, inputSpatialSize),
            "truth": scn.InputBatch(3, inputSpatialSize)
        }

        def add_sample(indices, features, truth):
            batch["data"].addSample()
            batch["truth"].addSample()
            indices = torch.LongTensor(indices)
            try:
                batch["data"].setLocations(indices, torch.FloatTensor(features), 0) #Use 1 to remove duplicate coords?
                batch["truth"].setLocations(indices, torch.FloatTensor(truth), 0)
            except AssertionError:
                #PDB didn't fit in grid?
                pass
            #del features
            #del truth

    for d in data:
    	if d["data"] is None: continue
        add_sample(d["indices"], d["data"], d["truth"])

    #print "Made batch"
    if create_tensor:
        batch["data"].precomputeMetadata(1)

    return batch

class IBISDataset(Dataset):
    def __init__(self, ibis_data, transform=True, input_shape=(96, 96, 96), tax_glob_group="A_eukaryota", num_representatives=2, only_aa=False, start_index=0, end_index=None, train=True):
        self.transform = transform
        self.only_aa = only_aa
        self.input_shape = input_shape
        self.epoch = None
        self.batch = None

        # Open and load text file including the whole training data
        data = pd.read_table(ibis_data, sep="\t")
        self.data = data.loc[(data["tax_glob_group"] == tax_glob_group) & (data["n"] >= num_representatives)]

        try:
            skip = "{}_skip.tab".format(os.path.splitext(ibis_data)[0])
            with open(skip) as skip_f:
                skip_ids = [int(l.rstrip()) for l in skip_f if l]
            osize = self.data.shape[0]
            self.data= self.data.loc[~self.data["unique_obs_int"].isin(skip_ids)]
            print "{}: Skipping {} of {}, {} remaining of {}".format("Train" if train else "Validate", len(skip_ids), osize, self.data.shape[0], osize)
        except IOError:
            print "No Skip ID file"
            pass

        if 0 < start_index < 1:
            start_index *= self.data.shape[0]

        if end_index is None:
            end_index = self.data.shape[0]
        elif end_index < 1:
            end_index *= self.data.shape[0]

        start_index = int(start_index)
        end_index = int(end_index)

        if start_index != 0 or end_index != self.data.shape[0]:
            self.data = self.data.iloc[start_index:end_index]

        self.n_samples = self.data.shape[0]
        self.train = train

    @classmethod
    def get_training_and_validation(cls, ibis_data, transform=True, input_shape=(96, 96, 96), tax_glob_group="A_eukaryota", num_representatives=2, data_split=0.8, only_aa=False, train_full=False, validate_full=False):
        print "Train full", train_full, "Validate full", validate_full
        train = cls(
            ibis_data,
            transform=transform,
            input_shape=input_shape,
            tax_glob_group=tax_glob_group,
            num_representatives=num_representatives,
            end_index=data_split,
            only_aa=only_aa,
            train=train_full)
        validate = cls(
            ibis_data,
            transform=transform,
            input_shape=input_shape,
            tax_glob_group=tax_glob_group,
            num_representatives=num_representatives,
            start_index=data_split,
            only_aa=only_aa,
            train=validate_full)
        return {"train":train, "val":validate}

    def get_data_loader(self, batch_size, shuffle, num_workers):
        return DataLoader(self, batch_size=batch_size, \
            shuffle=shuffle, num_workers=num_workers, collate_fn=sparse_collate)

    def __getitem__(self, index, verbose=True):
        datum = self.data.iloc[index]
        #print "Running {} ({}.{}): {}".format(datum["unique_obs_int"], datum["pdb"], datum["chain"], ",".join(["{}{}".format(i,n) for i, n in zip(datum["resi"].split(","), datum["resn"].split(","))]))
        try:
            inidices, data, truth = Structure.features_from_string(
                datum["pdb"],
                datum["chain"],
                datum["resi"],
                id=datum["unique_obs_int"],
                input_shape=self.input_shape,
                rotate=self.transform,
                only_aa=self.only_aa,
                expand_features=False,
                include_full_protein=self.train)
        except (KeyboardInterrupt, SystemExit):
            raise
        except InvalidPDB:
            trace = traceback.format_exc()
            print "Error:", trace
            with open("{}_{}_{}_{}_{}.error".format(datum["pdb"], datum["chain"], datum["unique_obs_int"], self.epoch, self.batch), "w") as ef:
                print >> ef, trace

            #return
            return {"indices": None,
                    "data": None,
                    "truth": None
                    }
        except:
            trace = traceback.format_exc()
            print "Error:", trace
            with open("{}_{}_{}_{}_{}.error".format(datum["pdb"], datum["chain"], datum["unique_obs_int"], self.epoch, self.batch), "w") as ef:
                print >> ef, trace
            raise

        #print "Got data for", index

        data = np.nan_to_num(data)

        #print truth
        #print truth.shape

        grid_indices = []
        grid_features = []
        grid_truth = []
        for grid_index, sample_indices in groupby(sorted(enumerate(indices.tolist()), key=lambda x:x[1]), key=lambda x:x[1]):
            sample_indices, _ = zip(*sample_indices)
            sample_indices = list(sample_indices)
            # print "SAME GRIDS"
            #print data[sample_indices].shape
            # features = np.sum(data[sample_indices], axis=0)
            # #assert np.argmax(features) < data[sample_indices].shape[1]
            # features = np.eye(1, 21, np.argmax(features))[0].tolist() #One-hot
            # grid_indices.append(grid_index)
            # grid_features.append(features)
            # print "NEW FEATURES"
            # print features

            if verbose and len(sample_indices) > 1:
                print "More than one atom per voxel in", datum["pdb"], datum["chain"], datum["unique_obs_int"]
                print "   ", len(sample_indices), "at", grid_index, "but contain the same residue" if data[sample_indices[0]].tolist()==data[sample_indices[1]].tolist() else "but do not contain the same residue"
                print "   Using the first found atom"

            grid_indices.append(grid_index)
            grid_features.append(data[sample_indices[0]].tolist())
            grid_truth.append(truth[sample_indices[0]].tolist())

        sample = {
            "indices": grid_indices,
            "data": grid_features,
            "truth": grid_truth #np.ones((len(grid_indices), 1), dtype=int).tolist()
            }

        return sample

    # Override to give PyTorch size of dataset
    def __len__(self):
        return self.data.shape[0]

class SphereDataset(Dataset):
    def __init__(self, shape, cnt=5, r_min=10, r_max=30, border=10, sigma=20, n_samples=1000, train=True):
        self.shape = shape
        self.cnt = cnt
        self.r_min = r_min
        self.r_max = r_max
        self.border = border
        self.sigma = sigma
        self.n_samples = n_samples
        self.train = train

    @classmethod
    def get_training_and_validation(cls, shape, cnt=5, r_min=10, r_max=30, border=10, n_samples=1000, data_split=0.8):
        train = cls(shape, cnt=cnt, r_min=r_min, r_max=r_max, border=border, n_samples=n_samples*data_split, train=True)
        validate = cls(shape, cnt=cnt, r_min=r_min, r_max=r_max, border=border, n_samples=n_samples*(1-data_split), train=False)
        return {"train":train, "val":validate}

    def get_data_loader(self, batch_size, shuffle, num_workers):
        return DataLoader(self, batch_size=batch_size, \
            shuffle=shuffle, num_workers=num_workers, collate_fn=sparse_collate if self.train else dense_collate)

    # Override to give PyTorch size of dataset
    def __len__(self):
        return self.n_samples

    def __getitem__(self, index):
        indices, data, truth = create_spheres(self.cnt, self.shape, self.border, self.r_min, self.r_max, train=self.train)

        sample = {
            "indices": indices,
            "data": data,
            "truth": truth
            }

        return sample

def create_spheres(num_spheres, shape=(96, 96, 96), border=10, min_r=5, max_r=15, hallow=True, train=True):
    """Create randomly placed and randomy sized spheres inside of a grid
    """
    if train:
        class sphere:
            indices = []
            volume = []
            labels = []
        def add_indices(_indices, value, truth=1):
            ind = _indices.T.tolist()
            sphere.indices += ind
            sphere.volume += [[value]]*len(ind)
            sphere.labels += [[truth]]*len(ind)
    else:
        class sphere:
            indices = []
            volume = np.random.random(list(shape))
            labels = np.zeros(list(shape))
        def add_indices(_indices, value, truth=1):
            sphere.indices += _indices.T.tolist()
            sphere.volume[_indices] = value
            sphere.labels[_indices] = truth

    for i in xrange(num_spheres):
        #Define random center of sphere and radius
        center = [np.random.randint(border, edge-border) for edge in shape]
        r = np.random.randint(min_r, max_r)
        color = np.random.random()

        y, x, z = np.ogrid[-center[0]:shape[0]-center[0], -center[1]:shape[1]-center[1], -center[2]:shape[2]-center[2]]
        m = x*x + y*y + z*z < r*r

        if hallow:
            n = x*x + y*y + z*z > (r-1)*(r-1)
            sphere_indices = np.array(np.where((m==True)&(n==True)))
        else:
            sphere_indices = np.array(np.where(m==True))

        add_indices(sphere_indices, color)

    return sphere.indices, sphere.volume, sphere.labels
