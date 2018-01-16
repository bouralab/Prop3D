import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/usr/share/pdb2pqr")
sys.path.append("/data/draizene/seaborn")

import os
import time
import multiprocessing
import math
from datetime import datetime
from itertools import izip

import matplotlib
matplotlib.use("Agg")

import numpy as np
from sklearn import metrics

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

from matplotlib.backends.backend_pdf import PdfPages

from torch.optim import Adam, SGD
import torch
from torch.autograd import Variable
from torch.nn.modules.loss import _Loss
from torch.optim.lr_scheduler import StepLR, LambdaLR
import sparseconvnet as scn

from molmimic.torch_model.torch_model import UNet3D
from molmimic.torch_model.torch_loader import IBISDataset

import subprocess

def get_gpu_memory_map():
    """Get the current gpu usage.

    Returns
    -------
    usage: dict
        Keys are device ids as integers.
        Values are memory usage as integers in MB.
    """
    result = subprocess.check_output(
        [
            'nvidia-smi', '--query-gpu=memory.used',
            '--format=csv,nounits,noheader'
        ])
    # Convert lines into a dictionary
    gpu_memory = [int(x) for x in result.strip().split('\n')]
    gpu_memory_map = dict(zip(range(len(gpu_memory)), gpu_memory))
    return gpu_memory_map
    #return {torch.cuda.getMemoryUsage(i) for i in xrange(torch.cuda.device_count())}

class ModelStats(object):
    all_dice = {"train":[], "val":[]}
    all_accuracies = {"train":[], "val":[]}
    all_fpr = {"train":[], "val":[]}
    all_tpr = {"train":[], "val":[]}

    def __init__(self):
        self.top1 = 0
        self.top5 = 0
        self.n = 0
        self.nll = 0
        self.running_corrects = None
        self.true_corrects = None
        self.running_corrects_num = 0 #TP+TN
        #self.confusion_matrix = torch.FloatTensor(shape).zero_()

        self.accuracies = []
        self.losses = []

    def update(self, output, target, loss):
        batchSize = output.size(0)
        #import pdb; pdb.set_trace()

        self.n += batchSize
        self.nll += loss * batchSize

        predicted_corrects = output >=0.9
        predicted_corrects = predicted_corrects.float()
        predicted_corrects_num = predicted_corrects.eq(target).sum()
        self.running_corrects_num += predicted_corrects_num #TP+TN
        #self.running_incorrects += batchSize-predicted_corrects #FP+FN
        if self.running_corrects is None:
            self.running_corrects = predicted_corrects.cpu().view(-1)
            self.true_corrects = target.cpu().view(-1)
        else:
            self.running_corrects = torch.cat((self.running_corrects, predicted_corrects.cpu().view(-1)), 0)
            self.true_corrects = torch.cat((self.true_corrects, target.cpu().view(-1)), 0)

        self.accuracies.append(predicted_corrects_num/float(batchSize))
        self.losses.append(loss)


    def save(self, phase, epoch):
        ModelStats.all_dice[phase] += self.losses
        ModelStats.all_accuracies[phase].append(self.running_corrects_num/float(self.n))
        self.plot(phase, epoch)

    def top1pct(self):
        return 100 * (1 - 1.0 * self.top1 / float(self.n))

    def top5pct(self):
        return 100 * (1 - 1.0 * self.top5 / float(self.n))

    def nllpct(self):
        return 100*self.nll/float(self.n)

    def correctpct(self):
        return 100*self.running_corrects_num/float(self.n)

    def plot_accuracy(self, ax, final_phase=None):
        ax.plot(self.accuracies if final_phase is None else ModelStats.all_accuracies[final_phase])
        ax.set_title("Accuracy Increases per {}".format("Batch" if final_phase is None else "Epoch"))
        ax.set_xlabel("Batch #" if final_phase is None else "Epoch #")
        ax.set_ylabel("Accuracy")
        ax.set_ylim([0, 1])

    def plot_loss(self, ax, final_phase=None):
        ax.plot(self.losses if final_phase is None else ModelStats.all_dice[final_phase])
        ax.set_title("Dice Coefficent Decreases per {}".format("Batch" if final_phase is None else "Epoch"))
        ax.set_xlabel("Batch #" if final_phase is None else "Epoch #")
        ax.set_ylabel("Loss")
        ax.set_ylim([-1, 0])

    def plot_roc(self, ax, phase, final=False):
        if not final:
            #import pdb; pdb.set_trace()
            fpr, tpr, _ = metrics.roc_curve(
                self.running_corrects.numpy(), self.true_corrects.numpy(), pos_label=1.)
            ModelStats.all_fpr[phase] += fpr.tolist()
            ModelStats.all_tpr[phase] += tpr.tolist()
        else:
            fpr = ModelStats.all_fpr[phase]
            tpr = ModelStats.all_tpr[phase]

        rocauc = metrics.auc(fpr, tpr, reorder=True)

        ax.plot(fpr, tpr, "-", label="AUC: {:4f}".format(rocauc))
        ax.set_title('ROC')
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_ylim([0, 1])
        legend = ax.legend(loc='lower right')

    def plot(self, phase, epoch):
        pp = PdfPages('epoch{}_{}_statistics.pdf'.format(epoch, phase))
        f, axes = plt.subplots(1, 2, figsize=(16,6))
        f.suptitle("{} Epoch {} Statistics".format(phase.title(), epoch), fontsize=14)
        self.plot_accuracy(axes[0])
        self.plot_loss(axes[1])
        #self.plot_roc(axes[2], phase)
        f.subplots_adjust(wspace=.4)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)

    def plot_final(self):
        pp = PdfPages('final_train_statistics.pdf')
        f, axes = plt.subplots(1, 2)
        f.suptitle("Final Train Statistics", fontsize=14)
        self.plot_accuracy(axes[0], final_phase="train")
        self.plot_loss(axes[1], final_phase="train")
        #self.plot_roc(axes[2], phase="train", final=True)
        f.subplots_adjust(wspace=.4)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)

        pp = PdfPages('final_validation_statistics.pdf')
        f, axes = plt.subplots(1, 2)
        f.suptitle("Final Validation Statistics", fontsize=14)
        self.plot_accuracy(axes[0], final_phase="val")
        self.plot_loss(axes[1], final_phase="val")
        #self.plot_roc(axes[2], phase="val", final=True)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)

def train(ibis_data, input_shape=(96,96,96), model_prefix=None, check_point=True, save_final=True, only_aa=False, num_workers=None, num_epochs=30, batch_size=20, shuffle=True, use_gpu=True, initial_learning_rate=0.0001, learning_rate_drop=0.5, learning_rate_epochs=10, lr_decay=4e-2, data_split=0.8, train_full=False, validate_full=False, no_batch_norm=False):
    if model_prefix is None:
        model_prefix = "./molmimic_model_{}".format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))

    if num_workers is None:
        num_workers = multiprocessing.cpu_count()-1
    print "Using {} workers".format(num_workers)

    since = time.time()

    if ibis_data == "spheres":
        from torch_loader import SphereDataset
        datasets = SphereDataset.get_training_and_validation(input_shape, cnt=1, n_samples=1000, data_split=0.99)
        nFeatures = 1
        validation_batch_size = 1
        input_shape = (96, 96, 96)
    elif os.path.isfile(ibis_data):
        datasets = IBISDataset.get_training_and_validation(
            ibis_data,
            input_shape=input_shape,
            only_aa=only_aa,
            data_split=data_split,
            train_full=train_full,
            validate_full=validate_full)
        nFeatures = 21 if only_aa else 59
        validation_batch_size = batch_size
    else:
        raise RuntimeError("Invalid training data")

    dataloaders = {name:dataset.get_data_loader(
        batch_size if dataset.train else validation_batch_size,
        shuffle,
        num_workers) \
        for name, dataset in datasets.iteritems()}

    dtype = 'torch.cuda.FloatTensor' if torch.cuda.is_available() else 'torch.FloatTensor'

    model = UNet3D(nFeatures, 1, batchnorm=not no_batch_norm)
    model.type(dtype)
    #optimizer = Adam(model.parameters(), lr=initial_learning_rate)
    optimizer = SGD(model.parameters(),
        lr = initial_learning_rate,
        momentum = 0.9,
        weight_decay=1e-4,
        nesterov=True)

    if False:
        scheduler = StepLR(optimizer, step_size=1, gamma=learning_rate_drop)
    else:
        scheduler = LambdaLR(optimizer, lambda epoch: math.exp((1 - epoch) * lr_decay))

    check_point_model_file = "{}_checkpoint_model.pth".format(model_prefix)
    check_point_epoch_file = "{}_checkpoint_epoch.pth".format(model_prefix)
    if check_point and os.path.isfile(check_point_model_file) and os.path.isfile(check_point_epoch_file):
        start_epoch = torch.load(check_point_epoch_file)
        print "Restarting at epoch {} from {}".format(start_epoch+1, check_point_file)
        model.load_state_dict(torch.load(check_point_model_file))
    else:
        start_epoch = 0

    criterion = DiceLoss()

    inputSpatialSize = torch.LongTensor(input_shape)

    for epoch in xrange(start_epoch, num_epochs):
        print "Epoch {}/{}".format(epoch, num_epochs - 1)
        print "-" * 10

        #print get_gpu_memory_map()

        # Each epoch has a training and validation phase
        for phase in ['train', 'val']:
            datasets[phase].epoch = epoch
            num_batches = int(np.ceil(len(datasets[phase])/float(batch_size if phase == "train" else validation_batch_size)))

            if phase == 'train':
                scheduler.step()
                model.train(True)  # Set model to training mode
            else:
                model.train(False)  # Set model to evaluate mode

            stats = ModelStats()

            # Iterate over data.
            for data_iter_num, data in enumerate(dataloaders[phase]):
                #Input will be turned SparseConvNetTensor
                print "{} Batch: {} of {}".format(phase.title(), data_iter_num, num_batches),
                datasets[phase].batch = data_iter_num
                #print type(data["data"]), data["data"].__class__.__name__, data["data"].__class__.__name__ == "InputBatch"
                if data["data"].__class__.__name__ == "InputBatch":
                    sparse_input = True
                    inputs = data["data"]
                    labels = data["truth"]
                    if use_gpu:
                        inputs = inputs.cuda().to_variable(requires_grad=True)
                        labels = labels.cuda().to_variable()
                    else:
                        inputs = inputs.to_variable(requires_grad=True)
                        labels = labels.to_variable()

                elif isinstance(data["data"], (list, tuple)):
                    sparse_input = True
                    inputs = scn.InputBatch(3, inputSpatialSize)
                    labels = scn.InputBatch(3, inputSpatialSize)

                    for sample, (indices, features, truth) in enumerate(izip(data["indices"], data["data"], data["truth"])):
                        inputs.addSample()
                        labels.addSample()

                        indices = torch.LongTensor(indices)
                        features = torch.FloatTensor(features)
                        truth = torch.FloatTensor(truth)

                        try:
                            inputs.setLocations(indices, features, 0) #Use 1 to remove duplicate coords?
                            labels.setLocations(indices, truth, 0)
                        except AssertionError:
                            #PDB didn't fit in grid?
                            continue

                    del data

                    inputs.precomputeMetadata(1)

                    if use_gpu:
                        inputs = inputs.cuda().to_variable(requires_grad=True)
                        labels = labels.cuda().to_variable()
                    else:
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

                # zero the parameter gradients
                optimizer.zero_grad()

                # forward
                #print inputs
                outputs = model(inputs)

                if sparse_input:
                    loss = criterion(outputs.features, labels.features)

                    if math.isnan(loss.data[0]):
                        print "Loss is Nan?"
                        import pdb; pdb. set_trace()

                    stats.update(outputs.features.data, labels.features.data, loss.data[0])
                else:
                    outputs = scn.SparseToDense(3, 1)(outputs)
                    loss = criterion(outputs.cpu(), labels.cpu())
                    stats.update(outputs.data.cpu().view(-1), labels.data.cpu().view(-1), loss.data[0])

                # backward + optimize only if in training phase
                if phase == 'train':
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()

                print "Epoch {} {}: corrects:{:.2f}% nll:{:.2f}% dice:{:.4f}% time:{:.1f}s".format(
                    epoch, phase, stats.correctpct(), stats.nllpct(), loss.data[0]*-100, time.time() - since)

            stats.save(phase, epoch)

            del inputs
            del labels
            del outputs

            if epoch < num_epochs-1:
                del stats

            if check_point:
                torch.save(epoch, check_point_epoch_file)
                torch.save(model.state_dict(), check_point_model_file)

            model.set_log_level(1)

            #print get_gpu_memory_map()

    stats.plot_final()

    time_elapsed = time.time() - since
    print 'Training complete in {:.0f}m {:.0f}s'.format(time_elapsed/60, time_elapsed % 60)

    if save_final:
        torch.save(model.state_dict(), "{}.pth".format(model_prefix))

    return model

class DiceLoss(_Loss):
    def __init__(self, size_average=True, smooth=1.):
        super(DiceLoss, self).__init__(size_average)
        self.smooth = smooth

    def forward(self, input, target):
        return -self.dice_coef(input, target)

    def dice_coef(self, input, target):
        iflat = input.view(-1)
        tflat = target.view(-1)
        intersection = (iflat * tflat).sum()

        #print iflat
        #print tflat
        #print intersection

        dice = ((2. * intersection + self.smooth) / ((iflat.sum() + tflat.sum() + self.smooth)))



        return dice

class IoULoss(_Loss):
    def __init__(self, size_average=True, smooth=1.):
        super(IoULoss, self).__init__(size_average)
        self.smooth = smooth

    def forward(self, input, target):
        return -self.IoU(input, target)

    def IoU(self, input, target):
        # y_pred_f = input.view(input.numel())
        # y_true_f = target.view(target.numel())
        # intersection = torch.sum(y_true_f*y_pred_f)
        # dice = (2. * intersection + self.smooth)/(torch.sum(y_true_f) + torch.sum(y_pred_f) + self.smooth)
        # print dice
        # return dice
        iflat = input.view(-1)
        tflat = target.view(-1)
        intersection = (iflat * tflat).sum()

        return ((intersection + self.smooth) / ((iflat.sum() + tflat.sum() + intersection + self.smooth)))

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")
    parser.add_argument(
        "-s",
        "--shape",
        type=int,
        nargs=3,
        default=(512,512,512)) #(256,256,256))
    parser.add_argument(
        "--batch-size",
        type=int,
        default=20)
    parser.add_argument(
        "--epochs",
        type=int,
        default=30)
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=0.1)
    parser.add_argument(
        "--learning-rate-drop",
        type=float,
        default=0.5)
    parser.add_argument(
        "--learning-rate-epochs",
        type=int,
        default=5)
    parser.add_argument(
        "--data-split",
        type=float,
        default=0.8)
    parser.add_argument(
        "--no-shuffle",
        default=False,
        action="store_true",
        help="Do not shuffle data")
    parser.add_argument(
        "--only-aa",
        default=False,
        action="store_true",
        help="Only use one feature: aa (20 features since aa is one hot encoded. Else use all 59 features.")
    parser.add_argument(
        "--train-full",
        default=False,
        action="store_true",
        help="Train the network using full protein rather than just the binding site"
    )
    parser.add_argument(
        "--validate-full",
        default=False,
        action="store_true",
        help="Validate the network using full protein rather than just the binding site"
    )
    parser.add_argument(
        "--no-batch-norm",
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

    cpus = parser.add_mutually_exclusive_group()
    cpus.add_argument(
        "--num_cpus",
        type=int,
        default=1)
    cpus.add_argument(
        "--all_cpus",
        action="store_true",
        default=False)

    parser.add_argument(
        "ibis_data")

    args = parser.parse_args()

    if args.all_gpus:
        args.num_gpus = len(get_available_gpus())

    if args.all_cpus:
        args.num_cpus = None

    return args

if __name__ == "__main__":
    args = parse_args()

    train(
        args.ibis_data,
        input_shape           = args.shape,
        only_aa               = args.only_aa,
        num_workers           = args.num_cpus,
        num_epochs            = args.epochs,
        batch_size            = args.batch_size,
        shuffle               = not args.no_shuffle,
        use_gpu               = args.num_gpus > 0,
        initial_learning_rate = args.learning_rate,
        learning_rate_drop    = args.learning_rate_drop,
        learning_rate_epochs  = args.learning_rate_epochs,
        data_split            = args.data_split,
        train_full            = args.train_full,
        validate_full         = args.validate_full,
        no_batch_norm         = args.no_batch_norm
    )
