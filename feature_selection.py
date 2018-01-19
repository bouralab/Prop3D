import sys
sys.path.append("/data/draizene/3DUnetCNN")
sys.path.append("/usr/share/pdb2pqr")

import os
import argparse

import itertools as it

import numpy as np

import h5py

from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LassoCV, SGDClassifier
from sklearn.decomposition import PCA
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel

from biopdbtools import Structure, InvalidPDB
from pdb_generator import IBISGenerator

from Bio import PDB

from map_residues import map_residues

def linear_selection():
    lsvc = LinearSVC(C=0.01, penalty="l1", dual=False)
    # sgd = SGDRegressor(alpha=0.0001, average=False, epsilon=0.1, eta0=0.01,
 #       fit_intercept=True, l1_ratio=0.15, learning_rate='invscaling',
 #       loss='squared_loss', max_iter=None, n_iter=None, penalty='l1',
 #       power_t=0.25, random_state=None, shuffle=True, tol=None,
 #       verbose=0, warm_start=False)
    model = SelectFromModel(lsvc)
    return model, True

def variance_selection():
    sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
    return sel, False

def lasso_selection():
    # We use the base estimator LassoCV since the L1 norm promotes sparsity of features.
    clf = LassoCV()
    sfm = SelectFromModel(clf, threshold=0.25)
    return sfm, False

def inc_pca():
    return PCA(), False

def feature_selection(ibis_data, batch_size=100):
    data = IBISGenerator(ibis_data, input_shape=(96,96,96,59))
    feature_names = Structure.get_feature_names()

    models = {
         "linear":linear_selection(),
         "variance":variance_selection(),
         "lasso":lasso_selection(),
         "pca":inc_pca()
    }
    if False and os.path.isfile("collpased_features.npy"):
        X = np.memmap("collpased_features.npy", mode="r")
    else:
        #X = np.require(np.memmap("collpased_features.npy", mode="w+", shape=(data.data.shape[0], Structure.nFeatures)), requirements=['O'])
        X = []
        for i, batch_start in enumerate(xrange(0, data.data.shape[0], batch_size)):
            for j, index in enumerate(xrange(batch_start, min(batch_start+batch_size, data.data.shape[0]))):
                row = data.data.iloc[index]
                id = "{}_{}_{}".format(row["pdb"].lower(), row["chain"].split("_")[0], row["unique_obs_int"])
                precalc_features_path = os.path.join(
                    os.path.join(os.environ.get("MOLMIMIC_FEATURES", "/data/draizene/molmimic/features"),
                    "{}.h5".format(id)))
                try:
                    precalc_features = h5py.File(precalc_features_path, 'r')["features"][:]
                except IOError:
                    continue
                X += np.nan_to_num(precalc_features).tolist()
        X = np.array(X)

    for name, (model, should_create) in models.iteritems():
    	try:
        	model.fit(X.tolist())
        except ValueError as e:
        	print name
        	print e
        	print X
        	print
        	continue
        if should_create:
        	model = SelectFromModel(model, prefit=True)
        if name != "pca":
            print name, "=", [feature_names[i] for i, s in enumerate(model.get_support()) if s]

    cumsum = np.cumsum(models["pca"][0].explained_variance_ratio_)
    d = np.argmax(cumsum > 0.95) + 1
    print "pca, n-components: {}; explained variance: {}".format(d, models["pca"][0].explained_variance_ratio_)

def test_generator(ibis_data, max_iter=100):
    data = IBISGenerator(ibis_data, input_shape=(96,96,96))
    for i, (X, y) in enumerate(data.generate()):
        if i>=max_iter: break
        print i
        print "X=", X
        print "y=", y
        print

def test_features(ibis_data):
    data = IBISGenerator(ibis_data, input_shape=(96,96,96))
    for i, datum in data.data.loc[data.data["pdb"]=="3J7A"].iterrows():
        print "Running {} ({}.{}): {}".format(datum["unique_obs_int"], datum["pdb"], datum["chain"], ",".join(["{}{}".format(i,n) for i, n in zip(datum["resi"].split(","), datum["resn"].split(","))]))
        Structure.features_from_string(datum["pdb"], datum["chain"], datum["resi"], id=datum["unique_obs_int"], input_shape=(96,96,96))

def match_resnames(ibis_data):
    data = IBISGenerator(ibis_data, input_shape=(96,96,96,52))
    for i, datum in data.data.loc[data.data["pdb"]=="3B13"].iterrows():
        print datum["pdb"], datum["chain"], datum["resi"], datum["resn"]
        for (resi, resn, pdb_resi, ncbi_resi), true_resn in it.izip(map_residues(datum["pdb"], datum["chain"], [datum["resi"]]), datum["resn"].split(",")):
            print "   ", resi, "{}={}".format(resn, true_resn), pdb_resi, ncbi_resi
        print


        # s = Structure(datum["pdb"], datum["chain"]).extract_chain()
        # resn = []
        # try:
        #     for r in s.align_seq_to_struc(datum["resi"], return_residue=True):
        #         try:
        #             resn.append(PDB.Polypeptide.three_to_one(r.get_resname()))
        #         except KeyError:
        #             resn.append("X")
        # except KeyError:
        #     continue

        # if ",".join(resn) != datum["resn"]:
        #     print "Error in {}.{}:{}".format(datum["pdb"], datum["chain"], datum["resi"])
        #     print "    ", s.starting_index_seq, s.starting_index_struc
        #     print "    ", ",".join(resn), datum["resn"]
        # else:
        #     "Match for {}.{}:{} = {}".format(datum["pdb"], datum["chain"], datum["resi"], datum["resn"])

def parse_args():
    parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")

    parser.add_argument("--test_gen", default=False, action="store_true")

    parser.add_argument(
        "ibis_data")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if args.test_gen:
        test_generator(args.ibis_data)
    else:
        feature_selection(args.ibis_data)
