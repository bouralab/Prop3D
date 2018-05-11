import sys
sys.path.append("/data/draizene/molmimic")

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
from scipy.spatial import cKDTree, distance

try:
    import sparseconvnet as scn
except:
    #Allow module to load without scn to read data
    scn = None

from molmimic.biopdbtools import Structure, InvalidPDB
from itertools import product, izip

def dense_collate(data, batch_size=1):
    batch = {"data":None, "truth":None, "scaling":1.0}

    num_true = 0
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
        num_true += np.sum(d["truth"])

    batch["weight"] = float(num_true)/len(data) #(256*256*256)

    return batch

def sparse_collate(data, input_shape=(256,256,256), create_tensor=False):
    if scn is None or not create_tensor:
        batch = {
            "indices": [],
            "data": [],
            "truth": [],
            "id": []
            }

        def add_sample(indices, features, truth, id):
            batch["indices"].append(indices)
            batch["data"].append(features)
            batch["truth"].append(truth)
            batch["id"].append(id)

    else:
        inputSpatialSize = torch.LongTensor(input_shape)
        batch = {
            "data": scn.InputBatch(3, inputSpatialSize),
            "truth": scn.InputBatch(3, inputSpatialSize),
            "id": []
        }

        def add_sample(indices, features, truth, id):
            batch["data"].addSample()
            batch["truth"].addSample()
            indices = torch.from_numpy(indices).long()
            try:
                batch["data"].setLocations(indices, torch.from_numpy(features).float(), 0) #Use 1 to remove duplicate coords?
                batch["truth"].setLocations(indices, torch.from_numpy(truth).float(), 0)
                batch["id"].append(id)
            except AssertionError:
                #PDB didn't fit in grid?
                pass
            #del features
            #del truth

    sample_weights = []
    batch_weight = None
    num_data = 0.0
    for i, d in enumerate(data):
        if d["data"] is None: continue
        if batch_weight is None:
            batch_weight = 0.0 if d["truth"].shape[1] == 2 else np.zeros(data[0]["truth"].shape[1])

        add_sample(d["indices"], d["data"], d["truth"], d["id"])
        if d["truth"].shape[1] == 2:
            num_true = np.sum(d["truth"][:, 0])
            true_prob = num_true/float(d["truth"].shape[0])
            sample_weights.append(np.array((1-true_prob, true_prob)))
            batch_weight += num_true
            num_data += d["truth"].shape[0]
        else:
            num_true = np.sum(d["truth"], axis=0)
            batch_weight += num_true
            sample_weights.append(num_true/float(d["truth"].shape[0]))
            num_data += d["truth"].shape[0]

    batch_weight /= float(num_data)

    #print "Made batch"
    if create_tensor:
        batch["data"].precomputeMetadata(1)

    if isinstance(batch_weight, float):
        batch["sample_weights"] = np.array(sample_weights)
        batch["weight"] = np.array([1-batch_weight, batch_weight]) #None #1.-float(num_true)/len(data) #(256*256*256)
    else:
        batch["sample_weights"] = np.array(sample_weights)
        batch["weight"] = batch_weight
    return batch

class IBISDataset(Dataset):
    def __init__(self, ibis_data, transform=True, input_shape=(264,264,264), tax_glob_group="A_eukaryota", num_representatives=2, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, course_grained=False, expand_atom=False, start_index=0, end_index=None, undersample=False, oversample=False, remove_multimers=False, cellular_organisms=False, random_features=None, train=True):
        self.transform = transform
        self.only_aa = only_aa
        self.only_atom = only_atom
        self.non_geom_features = non_geom_features
        self.use_deepsite_features = use_deepsite_features
        self.oversample = oversample
        self.undersample = undersample
        self.train = train
        self.course_grained = course_grained
        self.expand_atom = expand_atom
        self.input_shape = input_shape
        self.epoch = None
        self.batch = None

        self.random_features = random_features
        if self.random_features is not None:
            if random_features[1]:
                self.rfeats = np.array(list(product([0,1], repeat=random_features[0])))
            else:
                self.rfeats = np.eye(random_features[0])

            #Return blue (0,0,1) if nFeatures=3, else always 3rd feature is used or last feature if nFeature < 3
            if len(random_features) >= 3:
                self.r_bs = np.array(map(int, random_features[2]))
            else:
                self.r_bs = self.rfeats[min(2, random_features[0])]



        # Open and load text file including the whole training data
        dtype = {'pdb': str, 'chain': str}
        sep = "\t"
        data = pd.read_table(ibis_data, sep="\t", dtype=dtype)
        if len(data.columns) == 1:
            sep = ","
            data = pd.read_table(ibis_data, sep=",", dtype=dtype)
        if data.columns[0] != "pdb":
            if len(data.columns) == 5:
                #Full protein
                names = ["pdb", "chain", "resi", "is_multimer", "cdd"]
            else:
                #Individual chain
                names = ["pdb", "chain", "sdi", "domNo", "resi"]
            data = pd.read_table(ibis_data, sep=sep, header=None, names=names, dtype=dtype)

        data.columns = [col.lower() for col in data.columns]

        try:
            if cellular_organisms:
                self.full_data = data.loc[(data["tax_glob_group"] == tax_glob_group) | (data["tax_glob_group"] == "cellular")]
                self.full_data = self.full_data.loc[self.full_data["n"] >= num_representatives]
                print "cellular_organisms"
            else:
                self.full_data = data.loc[(data["tax_glob_group"] == tax_glob_group) & (data["n"] >= num_representatives)]
        except KeyError:
            self.full_data = data

        try:
            skip = "{}_skip.tab".format(os.path.splitext(ibis_data)[0])
            with open(skip) as skip_f:
                skip_ids = [int(l.rstrip()) for l in skip_f if l if not l.startswith("#") and len(l.rstrip())>0]
            osize = self.full_data.shape[0]
            skip_size = len(skip_ids)
            self.full_data = self.full_data.loc[~self.full_data["unique_obs_int"].isin(skip_ids)]
            self.full_data = self.full_data.loc[~self.full_data["gi"].isin(skip_ids)]
            if remove_multimers:
                multimers = [648375,660939,660938,668915,659919,513280,666570,535768,17338,654986,631439,656314,532931,665852,545875,513834,640379,511746,606239,509681,618577,527589,530024,630491,360637,511758,660132,660138,652064,651818,639149,630492,630493,634831,634826,647996,647997,643192,614273,614267,647151,617702,616098,516702,478916,461717,640385,661187,547577,541160,493657,493658,641974,541196,541453,447450,428933,514064,372904,372906,514011,514060,545252,401694,549817,299853,514062,514066,661242,635739,642241,614531,549815,549803,483520,483519,338893,174920,661197,554878,661250,622610,92851,563706,562371,632146,592870,658319,615442,542548,573967,123213,575823,621426,376703,576582,608295,579900,579888,533721,473097,504981,642845,647951,327972,642859,647950,267698,617700,627019,642808,642798,635782,227110,372997,284288,231335,231345,600901,231343,476220,476224,237050,243070,497904,453293,404236,404238,384495,395891]
                skip_size += len(multimers)
                self.full_data = self.full_data.loc[~self.full_data["unique_obs_int"].isin(multimers)]
            print "{}: Skipping {} of {}, {} remaining of {}".format("Train" if train else "Validate", skip_size, osize, self.full_data.shape[0], osize)
        except IOError:
            print "No Skip ID file"
            pass

        if 0 < start_index < 1:
            start_index *= self.full_data.shape[0]

        if end_index is None:
            end_index = self.full_data.shape[0]
        elif end_index < 1:
            end_index *= self.full_data.shape[0]

        start_index = int(start_index)
        end_index = int(end_index)

        if start_index != 0 or end_index != self.full_data.shape[0]:
            self.full_data = self.full_data.iloc[start_index:end_index]

        #MMDB splits separate molecules in one chain, get original chain back by removing "_"
        self.full_data.loc[:,["chain"]] = self.full_data["chain"].apply(lambda chain: chain.split("_")[0] if isinstance(chain, str) else chain)
        self.full_data.loc[:, ["pdb"]] = self.full_data["pdb"].apply(lambda pdb: pdb.replace("e+", "E") if len(pdb)==5 else pdb)

        if "sdi" not in self.full_data:
            self.data = self.full_data[["pdb", "chain"]].drop_duplicates()
        else:
            self.data = self.full_data
        self.data.loc[:, "include_negatives"] = True
        self.n_samples = self.data.shape[0]

        if self.oversample:
            d = self.data.copy()
            d.loc[:, "include_negatives"] = False
            d = pd.concat([d]*5, ignore_index=True)
            self.data = pd.concat((self.data, d), ignore_index=True)
            self.n_samples = self.data.shape[0]

    @classmethod
    def get_training_and_validation(cls, ibis_data, transform=True, input_shape=(264,264,264), tax_glob_group="A_eukaryota", num_representatives=2, data_split=0.8, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False, course_grained=False, expand_atom=False, oversample=False, undersample=False, random_features=None, cellular_organisms=False):
        train = cls(
            ibis_data,
            transform=transform,
            input_shape=input_shape,
            tax_glob_group=tax_glob_group,
            num_representatives=num_representatives,
            end_index=data_split,
            only_aa=only_aa,
            only_atom=only_atom,
            non_geom_features=non_geom_features,
            use_deepsite_features=use_deepsite_features,
            course_grained=course_grained,
            expand_atom=expand_atom,
            undersample=undersample,
            oversample=oversample,
            random_features=random_features,
            cellular_organisms=cellular_organisms,
            train=True)
        validate = cls(
            ibis_data,
            transform=transform,
            input_shape=input_shape,
            tax_glob_group=tax_glob_group,
            num_representatives=num_representatives,
            start_index=data_split,
            only_aa=only_aa,
            non_geom_features=non_geom_features,
            use_deepsite_features=use_deepsite_features,
            course_grained=course_grained,
            expand_atom=expand_atom,
            only_atom=only_atom,
            undersample=False,
            oversample=False,
            random_features=random_features,
            cellular_organisms=cellular_organisms,
            train=False)
        return {"train":train, "val":validate}

    def get_data_loader(self, batch_size, shuffle, num_workers):
        return DataLoader(self, batch_size=batch_size, \
            shuffle=shuffle, num_workers=num_workers, collate_fn=sparse_collate)

    def get_number_of_features(self):
        if self.random_features is not None:
            return self.random_features[0]
        return Structure.number_of_features(
            only_aa=self.only_aa,
            only_atom=self.only_atom,
            non_geom_features=self.non_geom_features,
            use_deepsite_features=self.use_deepsite_features,
            course_grained=self.course_grained
            )

    def __getitem__(self, index, verbose=True):
        pdb_chain = self.data.iloc[index]

        if "sdi" not in self.full_data:
            binding_sites = self.full_data.loc[(self.full_data["pdb"]==pdb_chain["pdb"])&(self.full_data["chain"]==pdb_chain["chain"])]
            binding_site_residues = ",".join([binding_site["resi"] for _, binding_site in binding_sites.iterrows()])
            gi = "{}.{}".format(pdb_chain["pdb"], pdb_chain["chain"]) #binding_sites["gi"].iloc[0] #"{}.{}".format(pdb_chain["pdb"], pdb_chain["chain"]) if self.course_grained else binding_sites["unique_obs_int"].iloc[0]
            sdi = None
            domain = None
        else:
            binding_site_residues = pdb_chain["resi"]
            gi = "{}.{}.d{}".format(pdb_chain["pdb"], pdb_chain["chain"], pdb_chain["domnum"])
            sdi = pdb_chain["sdi"]
            domain = pdb_chain["domnum"]

        try:
            indices, data, truth = Structure.features_from_string(
                pdb_chain["pdb"],
                pdb_chain["chain"],
                sdi = sdi,
                domain = domain,
                resi=binding_site_residues,
                #id=gi,
                input_shape=self.input_shape,
                rotate=self.transform,
                only_aa=self.only_aa,
                only_atom=self.only_atom,
                non_geom_features=self.non_geom_features,
                use_deepsite_features=self.use_deepsite_features,
                course_grained=self.course_grained,
                expand_atom=self.expand_atom,
                include_full_protein=pdb_chain["include_negatives"],
                undersample=self.undersample)
        except (KeyboardInterrupt, SystemExit):
            raise
        except InvalidPDB, RuntimeError:
            trace = traceback.format_exc()
            with open("{}_{}_{}_{}_{}.error".format(pdb_chain["pdb"], pdb_chain["chain"], gi, self.epoch, self.batch), "w") as ef:
                print trace
                print >> ef, trace

            #return
            return {"indices": None,
                    "data": None,
                    "truth": None,
                    "id": gi
                    }
        except:
            trace = traceback.format_exc()
            print "Error:", trace
            with open("{}_{}_{}_{}_{}.error".format(pdb_chain["pdb"], pdb_chain["chain"], gi, self.epoch, self.batch), "w") as ef:
                print >> ef, trace
            raise

        if self.random_features is not None:
            #Use random features from SphereDataset
            data = self.rfeats[np.random.choice(len(self.rfeats), indices.shape[0])]
            truth_indices = np.where(truth[:, 1]==1)
            data[truth_indices] = self.r_bs

        sample = {
            "indices": indices,
            "data": np.nan_to_num(data),
            "truth": truth,
            "id": gi
            }

        return sample

    # Override to give PyTorch size of dataset
    def __len__(self):
        return self.n_samples

class DenseSphereDataset(Dataset):
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
        np.random.seed()
        indices, data, truth = create_spheres(self.cnt, self.shape, self.border, self.r_min, self.r_max, train=self.train)

        sample = {
            "indices": indices,
            "data": data,
            "truth": truth
            }

        return sample

class SphereDataset(Dataset):
    def __init__(self, shape, cnt=3, r_min=15, r_max=30, border=10, sigma=20, n_samples=1000, nFeatures=3, allow_feature_combos=False, bs_feature=None, bs_feature2=None, bs_features=None, train=True, stripes=False):
        assert nFeatures > 0
        self.shape = np.array(shape)
        self.cnt = cnt
        self.r_min = r_min
        self.r_max = r_max
        self.border = border
        self.sigma = sigma
        self.n_samples = n_samples
        self.nFeatures = nFeatures
        self.train = train

        x = np.arange(0, self.shape[0])
        y = np.arange(0, self.shape[1])
        z = np.arange(0, self.shape[2])
        mx, my, mz = np.meshgrid(x, y, z)
        self.voxel_tree = cKDTree(zip(mx.ravel(), my.ravel(), mz.ravel()))

        if allow_feature_combos:
            self.features = np.array(list(product([0,1], repeat=nFeatures)))
        else:
            self.features = np.eye(self.nFeatures)

        if bs_features is not None:
            if not isinstance(bs_features, (list,tuple)):
                raise RuntimeError("Invalid bs_features")

            bs_features = [[np.array(map(int, bs_feat)) for bs_feat in bs_feature.split(",")] for bs_feature in bs_features]
            try:
                self.bs_color, self.bs_color2 = zip(*bs_features)
            except ValueError:
                self.bs_color = [f[0] for f in bs_features]

            self.n_classes = len(self.bs_color)+1
        else:
            #Return blue (0,0,1) if nFeatures=3, else always 3rd feature is used or last feature if nFeature < 3
            if isinstance(bs_feature, str):
                self.bs_color = np.array(map(int, bs_feature))
            elif isinstance(bs_feature, np.ndarray):
                self.bs_color = bs_feature
            else:
                self.bs_color = self.features[min(2, self.nFeatures-1)]

            if isinstance(bs_feature2, str):
                self.bs_color2 = np.array(map(int, bs_feature2))
            elif isinstance(bs_feature2, np.ndarray):
                self.bs_color2 = bs_feature2
            else:
                self.bs_color2 = None

            self.n_classes = 2

        if stripes:
            print "Making stripes"
            self.color_patch = make_stripes
        else:
            self.color_patch = alternate_colors

    @classmethod
    def get_training_and_validation(cls, shape, cnt=3, r_min=15, r_max=30, border=10, n_samples=1000, nFeatures=3, allow_feature_combos=False, bs_feature=None, bs_feature2=None, bs_features=None, stripes=False, data_split=0.8):
        train = cls(shape, cnt=cnt, r_min=r_min, r_max=r_max, border=border, n_samples=n_samples*data_split, nFeatures=nFeatures, allow_feature_combos=allow_feature_combos, bs_feature=bs_feature, bs_feature2=bs_feature2, bs_features=bs_features, stripes=False, train=True)
        validate = cls(shape, cnt=cnt, r_min=r_min, r_max=r_max, border=border, n_samples=n_samples*(1-data_split), nFeatures=nFeatures, allow_feature_combos=allow_feature_combos, bs_feature=bs_feature, bs_feature2=bs_feature2, bs_features=bs_features, stripes=False, train=False)
        return {"train":train, "val":validate}

    def get_data_loader(self, batch_size, shuffle, num_workers):
        return DataLoader(self, batch_size=batch_size, \
            shuffle=shuffle, num_workers=num_workers, collate_fn=sparse_collate)

    # Override to give PyTorch size of dataset
    def __len__(self):
        return self.n_samples

    def __getitem__(self, index):
        np.random.seed()

        #center = [np.random.randint(self.border+self.r_max, edge-border-self.r_max) for edge in shape]
        r = np.random.randint(self.r_min, self.r_max)
        center = np.round(self.shape/2.).astype(int)

        outer_sphere_points = self.voxel_tree.query_ball_point(center, r)
        inner_sphere_points = self.voxel_tree.query_ball_point(center, r-1)
        sphere_points = np.setdiff1d(outer_sphere_points, inner_sphere_points)
        indices = np.array([self.voxel_tree.data[i] for i in sphere_points])

        cnt = np.random.randint(1, self.cnt+1)

        tree = cKDTree(indices)

        colors = self.features[np.random.choice(len(self.features), indices.shape[0])] #self.color_patch(indices, self.bs_color, self.bs_color2) #
        truth = np.zeros((indices.shape[0], self.n_classes))
        truth[:, 0] = 1.

        used_points = None
        distances = None
        for bs_id in xrange(cnt):
            num_search = 0
            while True:
                idx = np.random.randint(0, indices.shape[0])
                bs_position = indices[idx]
                if used_points is None: break
                distances = distance.cdist(used_points, bs_position[None], 'euclidean')
                if len(np.where(distances <= 2)[0]) == 0:
                    break
                num_search += 1
            size = np.random.randint(5, 8)
            ball_indices = list(tree.query_ball_point(bs_position, r=size))
            points = [tree.data[idx] for idx in ball_indices]

            if self.bs_color is not None and isinstance(self.bs_color, (tuple,list)):
                class_num = np.random.choice(len(self.bs_color))
                bs_color = self.bs_color[class_num]
                if self.bs_color2 is not None:
                    bs_color2 = self.bs_color2[class_num]
                    colors[ball_indices, :] = self.color_patch(points, bs_color, bs_color2)
                else:
                    colors[ball_indices, :] = bs_color
                truth[ball_indices, class_num] = 1.
                truth[ball_indices, 0] = 0.

            else:
                if self.bs_color2 is not None:
                    colors[ball_indices, :] = self.color_patch(points, self.bs_color, self.bs_color2)
                else:
                    colors[ball_indices, :] = self.bs_color
                truth[ball_indices, :] = np.array([0., 1.])

            points = np.array(points)
            if used_points is None:
                used_points = points
            else:
                used_points = np.concatenate((used_points, points))
        del used_points
        del tree
        del distances

        sample = {
            "indices": indices,
            "data": colors,
            "truth": truth,
            "id":index
            }

        return sample

def alternate_colors(patch, color1, color2):
    """Make sure no voxel in the patch is the same color as its 8 neighest neighbors
    """
    tree = cKDTree(patch)
    features = np.tile(color1, (len(patch),1))
    for i, pnt in enumerate(patch):
        #Max distance should be sqrt(3) to account for voxels in one the 27 surrounding voxels,
        #but since it is circle like, there will only be 8 neighbors
        neighbors = tree.query(pnt, k=8, distance_upper_bound=1.42)
        for d, n in izip(*neighbors):
            if d == np.inf or tree.data[n].tolist() == tree.data[i].tolist(): continue
            if features[n].tolist() == features[i].tolist():
                features[n] = color2
    return features

def make_stripes(patch, color1, color2, gap_dist=np.pi/60.):
    """Make stripes on patch
    """
    #Start with all features of the patch as color1
    features = np.tile(color1, (len(patch),1))
    patch = np.array(patch)
    patch_spherical = to_spherical(patch-np.array((96,96,96)))
    psis = patch_spherical[:,2]
    psi_range = np.arange(np.min(psis), np.max(psis), gap_dist)

    for psi in psi_range:
        stripe_points = np.where((patch_spherical[:,2]>=psi-0.01)&(patch_spherical[:,2]<=psi+0.01))
        stripe = patch[stripe_points]
        features[stripe_points] = color2

    return features

def to_spherical(xyz):
    """(theta, phi, r)
    """
    ptsnew = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,1] = np.arctan2(xyz[:,1], xyz[:,0])
    ptsnew[:,2] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    return ptsnew

def to_cartesian(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

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
