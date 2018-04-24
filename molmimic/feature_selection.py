import sys
sys.path.append("/data/draizene/molmimic")

import os
import argparse
import itertools as it

import numpy as np

from sklearn.feature_selection import VarianceThreshold
from sklearn.linear_model import LassoCV, SGDClassifier
from sklearn.decomposition import PCA
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel

from molmimic.biopdbtools import Structure, InvalidPDB
from molmimic.torch_model.torch_loader import IBISDataset

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

def feature_selection(ibis_data):
    dataset = IBISDataset(ibis_data)
    feature_names = Structure.get_feature_names()

    models = {
         "linear":linear_selection(),
         "variance":variance_selection(),
         "lasso":lasso_selection(),
         "pca":inc_pca()
    }

    X = np.array((len(dataset), len(feature_names)))
    for i in xrange(len(dataset)):
        features = dataset[i]
        X[i] = np.nan_to_num(features["data"])

    for name, (model, should_create) in models.iteritems():
    	try:
        	model.fit(X)
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

def parse_args():
    parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")

    parser.add_argument(
        "ibis_data")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    feature_selection(args.ibis_data)
