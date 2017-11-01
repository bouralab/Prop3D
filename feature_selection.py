from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel

from biopdbtools import Structure
from pdb_generator import IBISGenerator

def linear_selection(X, y):
	lsvc = LinearSVC(C=0.01, penalty="l1", dual=False).fit(X, y)
	model = SelectFromModel(lsvc, prefit=True)
	model_features = model.get_support()
	return model, model_features

def variance_selection(X, y):
	sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
	sel.fit(X)
	sel_features = sel.get_support()
	return sel, model_features

def lasso_selection(X, y):
	# We use the base estimator LassoCV since the L1 norm promotes sparsity of features.
	clf = LassoCV()
	sfm = SelectFromModel(clf, threshold=0.25)
	sfm.fit(X, y)
	sfm_features = sfm.get_support()
	return sfm, sfm_features

def feature_selection():
	data = IBISGenerator(ibis_data, input_shape=(96,96,96,52))
	X = [Stucture.feature_from_string(datum["pdb"], datum["chain"], datum["resi"]) for i, row in data.data.iterrows()]
	y = [1]*len(features)

	feature_names = Structure.get_feature_names()

def parse_args():
	parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")

	parser.add_argument(
		"ibis_data")

	return parser.parse_args()

