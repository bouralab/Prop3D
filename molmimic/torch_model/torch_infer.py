import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/data/draizene/tqdm")
#sys.path.append("/usr/share/pdb2pqr")
#sys.path.append("/data/draizene/seaborn")
#sys.path.append("/data/draizene/pytorchviz")
#sys.path.append("/data/draizene/torchbiomed")

import os
import time
import math
import multiprocessing
from datetime import datetime
from itertools import izip

import numpy as np
import gc
import torch
from torch.optim import Adam, SGD
from torch.autograd import Variable
from torch.optim.lr_scheduler import StepLR, LambdaLR
from torchnet.logger import MeterLogger

import sparseconvnet as scn

from tqdm import tqdm

from molmimic.torch_model.Loss import DiceLoss
import molmimic.torch_model.Loss as LossFunctions
from molmimic.torch_model.torch_model import UNet3D, ResNetUNet
from molmimic.torch_model.torch_loader import IBISDataset, IBISUnclusteredDataset
from molmimic.torch_model.ModelStats import ModelStats, add_to_logger, format_meter, graph_logger

from molmimic.biopdbtools import Structure

from torchviz import dot
#import torchbiomed.loss as bioloss

def load_model(model_file, no_batch_norm=False, use_resnet_unet=True, dropout_depth=False, dropout_width=False, dropout_p=0.5, unclustered=False, nFeatures=None, only_aa=False, only_atom=False, non_geom_features=False, use_deepsite_features=False):
    model_prefix = "./molmimic_model_{}".format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))

    dtype = 'torch.cuda.FloatTensor' if torch.cuda.is_available() else 'torch.FloatTensor'

    if nFeatures is None or nFeatures<0:
        nFeatures = Structure.number_of_features(only_aa=only_aa, only_atom=only_atom, non_geom_features=non_geom_features, use_deepsite_features=use_deepsite_features)

    if use_resnet_unet:
        model = ResNetUNet(nFeatures, 2, dropout_depth=dropout_depth, dropout_width=dropout_width, dropout_p=dropout_p)
    else:
        model = UNet3D(nFeatures, 2, batchnorm=not no_batch_norm)

    states =  torch.load(model_file)
    model.load_state_dict(states)
    model.type(dtype)
    model.train(False)  # Set model to evaluate mode

    return model

def infer(model, data, input_shape=(256,256,256), use_gpu=True, course_grained=False):
    model_prefix = "./molmimic_model_{}".format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))

    dtype = 'torch.cuda.FloatTensor' if torch.cuda.is_available() else 'torch.FloatTensor'
    inputSpatialSize = torch.LongTensor(input_shape)

    labels = None

    mlog = MeterLogger(nclass=2, title="Sparse 3D UNet Inference")

    phase = "test"
    epoch = 0

    batch_weight = data.get("weight", None)
    if batch_weight is not None:
        batch_weight = torch.from_numpy(batch_weight).float().cuda()
    sample_weights = data.get("sample_weights", None)
    if sample_weights is not None:
        sample_weights = torch.from_numpy(sample_weights).float().cuda()
        use_size_average = False
        weight = sample_weights if use_size_average else batch_weight

    if data["data"].__class__.__name__ == "InputBatch":
        sparse_input = True
        inputs = data["data"]
        labels = data["truth"] if "truth" in data else None
        if use_gpu:
            inputs = inputs.cuda().to_variable(requires_grad=False)
            labels = labels.cuda().to_variable()
        else:
            inputs = inputs.to_variable(requires_grad=False)
            labels = labels.to_variable()

    elif isinstance(data["data"], (list, tuple)):
        sparse_input = True
        inputs = scn.InputBatch(3, inputSpatialSize)
        labels = scn.InputBatch(3, inputSpatialSize) if "truth" in data else None

        if isinstance(data["data"][0], np.ndarray):
            long_tensor = lambda arr: torch.from_numpy(arr).long()
            float_tensor = lambda arr: torch.from_numpy(arr).float()
        elif isinstance(data["data"][0], (list, tuple)):
            long_tensor = lambda arr: torch.LongTensor(arr)
            float_tensor = lambda arr: torch.FloatTensor(arr)
        else:
            raise RuntimeError("invalid datatype")

        for sample, (indices, features, truth) in enumerate(izip(data["indices"], data["data"], data["truth"])):
            inputs.addSample()
            if labels is not None: 
                labels.addSample()

            try:
                indices = long_tensor(indices)
                features = float_tensor(features)
                if labels is not None:
                    truth = float_tensor(truth)
            except RuntimeError as e:
                print e
                continue

            inputs.setLocations(indices, features, 0) #Use 1 to remove duplicate coords?
            if labels is not None: 
                labels.setLocations(indices, truth, 0)

        del data
        del indices
        del truth

        inputs.precomputeMetadata(1)

        if use_gpu:
            inputs = inputs.cuda()
            labels = labels.cuda()

        inputs = inputs.to_variable(requires_grad=True)
        labels = labels.to_variable()

    elif isinstance(data["data"], torch.FloatTensor):
        #Input is dense
        print "Input is Dense"
        sparse_input = False
        if use_gpu:
            inputs = inputs.cuda()
            labels = labels.cuda()
        inputs = Variable(data["data"], requires_grad=True)
        inputs = scn.DenseToSparse(3)(inputs)
        try:
            inputs = inputs.cuda().to_variable(requires_grad=True)
        except:
            pass
        labels = Variable(data["truth"])

    else:
        raise RuntimeError("Invalid data from dataset")

    # forward
    try:
        outputs = model(inputs)
    except AssertionError:
        print nFeatures, inputs
        raise

    if labels is None:
        return outputs
    else:
        loss_fn = torch.nn.CrossEntropyLoss(weight=weight)

        loss = loss_fn(outputs, torch.max(labels.features, 1)[1])

        mlog.update_loss(loss, meter='loss')
        mlog.update_meter(outputs, torch.max(labels.features, 1)[1], meters={'accuracy', 'map'})
        add_to_logger(mlog,  "Train" if phase=="train" else "Test", epoch, outputs, labels.features, batch_weight)

        del inputs
        del labels
        del loss
        del loss_fn
        del batch_weight
        del sample_weights

        return outputs, mlog

def create_args(parser=None, default_shape=(256,256,256), add_model=True, add_data=True):
    if parser is None:
        import argparse
        parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")

    #Data specific
    parser.add_argument(
        "-s",
        "--shape",
        type=int,
        nargs=3,
        default=default_shape)
    parser.add_argument(
        "--only-aa",
        default=False,
        action="store_true",
        help="Only use one feature: aa (20 features since aa is one hot encoded). Else use all 59 features.")
    parser.add_argument(
        "--only-atom",
        default=False,
        action="store_true",
        help="Only use one feature: atom type (5 features since atom is one hot encoded). Else use all 59 features.")
    parser.add_argument(
        "--non-geom-features",
        default=False,
        action="store_true",
        help="Only use non geometric features")
    parser.add_argument(
        "--use_deepsite_features",
        default=False,
        action="store_true",
        help="Only use DeepSite features")
    parser.add_argument(
        "--course-grained",
        default=False,
        action="store_true",
        help="Validate the network using full protein rather than just the binding site"
    )
    parser.add_argument(
        "--unclustered",
        default=False,
        action="store_true",
    )

    #Model specific
    parser.add_argument(
        "--no-batch-norm",
        default=False,
        action="store_true",
        help="Do not use BatchNorm after each conv layer"
    )
    parser.add_argument(
        "--use-resnet-unet",
        default=False,
        action="store_true",
        help="Do not use BatchNorm after each conv layer"
    )

    gpus = parser.add_mutually_exclusive_group()
    gpus.add_argument(
        "--num_gpus",
        type=int,
        default=1)
    gpus.add_argument(
        "--all_gpus",
        action="store_true",
        default=False)

    # cpus = parser.add_mutually_exclusive_group()
    # cpus.add_argument(
    #     "--num_cpus",
    #     type=int,
    #     default=1)
    # cpus.add_argument(
    #     "--all_cpus",
    #     action="store_true",
    #     default=False)

    if add_model:
        parser.add_argument("model")

    if add_data:
        parser.add_argument("data")

    return parser

def parse_args():
    args = parser.parse_args()

    # if args.all_gpus:
    #     args.num_gpus = len(get_available_gpus())

    # if args.all_cpus:
    #     args.num_cpus = None

    return args

if __name__ == "__main__":
    args = parse_args()

    model = load(args.model, no_batch_norm=args.no_batch_norm, use_resnet_unet=args.use_resnet_unet)

    train(
        args.ibis_data,
        input_shape           = args.shape,
        model_prefix          = args.prefix,
        only_aa               = args.only_aa,
        only_atom             = args.only_atom,
        non_geom_features     = args.non_geom_features,
        use_deepsite_features = args.use_deepsite_features,
        num_workers           = args.num_cpus,
        use_gpu               = args.num_gpus > 0,
        course_grained        = args.course_grained,
        no_batch_norm         = args.no_batch_norm,
        use_resnet_unet       = args.use_resnet_unet,
        unclustered           = args.unclustered,
    )
