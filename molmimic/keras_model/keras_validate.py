import sys
sys.path.append("/data/draizene/3DUnetCNN")
sys.path.append("/data/draizene/molmimic")
sys.path.append("/usr/share/pdb2pqr")

import argparse
from keras import backend as K
K.set_image_dim_ordering('tf')

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, matthews_corrcoef

from matplotlib.backends.backend_pdf import PdfPages

from unet3d.training import load_old_model

from pdb_generator import IBISGenerator

def validate(model_file, ibis_data, input_shape=(96,96,96), batch_size=2, data_split=0.8, num_gpus=None, only_aa=False):
    model = load_old_model(model_file)
    validate_data = IBISGenerator(ibis_data, input_shape=input_shape, batch_size=batch_size, only_aa=only_aa, start_index=data_split, shuffle=False)

    data_gen = validate_data.generate()

    y_true = []
    y_pred = []
    num_batches = validate_data.data.shape[0]/batch_size
    for i in xrange(num_batches):
        print i, "of", num_batches-1
        X, y = data_gen.next()
        y_true_sample = y.flatten().astype(int)
        y_pred_sample = model.predict_on_batch(X).flatten().astype(int)

        y_true += y_true_sample.tolist()
        y_pred += y_pred_sample.tolist()

        print np.where(y_true_sample == 1)
        print np.where(y_pred_sample == 1)

        tpr, fpr, _ = roc_curve(y_true_sample, y_pred_sample)
        roc_auc = auc(fpr, tpr)
        print "Batch ROCAUC:", roc_auc

    tpr, fpr, _ = roc_curve(y_true_sample, y_pred_sample.tolist())
    roc_auc = auc(fpr, tpr)
    print "Total TPR, FPR:", tpr, fpr

    model_name = os.path.splitext(os.path.basename(model_file))[0]
    pp = PdfPages("model_evaluation-{}.pdf".format(model_name))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(fpr,tpr, lw=3., label="ROC (AUC: {})".format(roc_auc))
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_xlim([0, 1.0])
    ax.set_ylim([0.0, 1.05])

    fig.suptitle("{} Model Evaluation".format(model_name), fontsize=20)

    pp.savefig()
    pp.close()

    return roc_auc


def parse_args():
    parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")
    parser.add_argument(
        "-s",
        "--shape",
        nargs=3,
        default=(96, 96, 96))
    parser.add_argument(
        "--batch",
        type=int,
        default=2)
    parser.add_argument(
        "--split",
        type=float,
        default=0.8)
    parser.add_argument(
        "--only-aa",
        default=False,
        action="store_true",
        help="Only use one feature: aa (20 features since aa is one hot encoded. Else use all 59 features.")

    gpus = parser.add_mutually_exclusive_group()
    gpus.add_argument(
        "--num_gpus",
        type=int,
        default=1)
    gpus.add_argument(
        "--all_gpus",
        action="store_true",
        default=False)

    parser.add_argument(
        "model_file")

    parser.add_argument(
        "ibis_data")

    args = parser.parse_args()

    if args.all_gpus:
        args.num_gpus = len(get_available_gpus())

    return args

if __name__ == "__main__":
    args = parse_args()
    validate(args.model_file, args.ibis_data, args.shape, args.batch, args.split, args.num_gpus, args.only_aa)