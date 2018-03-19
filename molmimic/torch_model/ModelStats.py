import time
from itertools import groupby

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

from matplotlib.backends.backend_pdf import PdfPages

from sklearn import metrics

import numpy as np

import torch

values = {}
epoch_values = {}
def add_to_logger(logger, phase, epoch, output, target, weight, locations=None, n_classes=2, smooth=1.0):
    weight_sum = weight.sum()

    dice_combined = 0.0
    weighted_dice_combined = 0.0
    mcc_combined = 0.0
    precision_combined = 0.0
    kappa_combined = 0.0
    f1_combined = 0.0

    for c in range(n_classes):
        iflat = output[:,c].cpu().view(-1)
        tflat = target[:,c].cpu().view(-1)
        intersection = (iflat * tflat).sum()

        dice = ((2. * intersection + smooth) / ((iflat.sum() + tflat.sum() + smooth))).data[0]
        weighted_dice = weight[c]*dice

        iflat = (iflat > .7).long()
        tflat = tflat.long()
        mcc = metrics.matthews_corrcoef(tflat.data.numpy(), iflat.data.numpy())
        precision = metrics.precision_score(tflat.data.numpy(), iflat.data.numpy())
        kappa = metrics.cohen_kappa_score(tflat.data.numpy(), iflat.data.numpy())
        f1 =  metrics.f1_score(tflat.data.numpy(), iflat.data.numpy())

        logger.update_loss(np.array([dice]),          meter='dice_class{}'.format(c))
        logger.update_loss(np.array([weighted_dice]), meter='weighted_dice_class{}'.format(c))
        logger.update_loss(np.array([mcc]),           meter='mcc_class{}'.format(c))
        logger.update_loss(np.array([precision]),     meter='precision_class{}'.format(c))
        logger.update_loss(np.array([kappa]),         meter='kappa_class{}'.format(c))
        logger.update_loss(np.array([f1]),            meter='f1_class{}'.format(c))

        dice_combined += dice
        weighted_dice_combined += weighted_dice
        mcc_combined += mcc
        precision_combined += precision
        kappa_combined += kappa
        f1_combined += f1

        del iflat
        del tflat
        del intersection
        del mcc
        del precision
        del kappa
        del f1

    logger.update_loss(np.array([dice_combined]),          meter='dice_sum')
    logger.update_loss(np.array([weighted_dice_combined]), meter='weighted_dice_sum')

    weighted_dice_combined /= weight_sum
    logger.update_loss(np.array([weighted_dice_combined]), meter='weighted_dice_wavg')

    dice_combined /= n_classes
    mcc_combined /= n_classes
    precision_combined /= n_classes
    kappa_combined /= n_classes
    f1_combined /= n_classes

    logger.update_loss(np.array([dice_combined]),      meter='dice_avg')
    logger.update_loss(np.array([mcc_combined]),       meter='mcc_avg')
    logger.update_loss(np.array([precision_combined]), meter='precision_avg')
    logger.update_loss(np.array([kappa_combined]),     meter='kappa_avg')
    logger.update_loss(np.array([f1_combined]),        meter='f1_avg')

    if n_classes == 2:
        result = (output>0.7).cpu().float()
        _, results_flat = result.max(dim=1)
        _, targets_flat = result.max(dim=1)

        iflat = results_flat.view(-1)
        tflat = targets_flat.view(-1)
        intersection = (iflat * tflat).sum()

        dice = ((2. * intersection + smooth) / ((iflat.sum() + tflat.sum() + smooth))).data[0]
        weighted_dice = weight_sum*dice
        mcc = metrics.matthews_corrcoef(tflat.data.numpy(), iflat.data.numpy())
        precision = metrics.precision_score(tflat.data.numpy(), iflat.data.numpy())
        kappa = metrics.cohen_kappa_score(tflat.data.numpy(), iflat.data.numpy())
        f1 = metrics.f1_score(tflat.data.numpy(), iflat.data.numpy())

        logger.update_loss(np.array([dice]),          meter='dice_flat')
        logger.update_loss(np.array([weighted_dice]), meter='weighted_dice_flat')
        logger.update_loss(np.array([mcc]),           meter='mcc_flat')
        logger.update_loss(np.array([precision]),     meter='precision_flat')
        logger.update_loss(np.array([kappa]),         meter='kappa_flat')
        logger.update_loss(np.array([f1]),            meter='f1_flat')

        del result
        del results_flat
        del targets_flat
        del weight_sum
        del iflat
        del tflat
        del intersection

    global values
    for meter_name, value in logger.meter.iteritems():
        val = value.value()
        val = val[0] if isinstance(val, tuple) and val[0] is not None else val
        val = val if val is not None else 0.0

        try:
            values[meter_name].append(val)
        except KeyError:
            values[meter_name] = [val]

        try:
            epoch_val = value.mean
        except AttributeError:
            ep_val = value.value()
            ep_val = ep_val[0] if isinstance(val, tuple) else ep_val
            ev_val = ep_val if ep_val is not None else 0.0
            epoch_val = ep_val[0] if isinstance(val, tuple) else ep_val
        try:
            epoch_values[meter_name][epoch] = epoch_val
        except KeyError:
            epoch_values[meter_name] = {epoch: epoch_val}

def format_meter(logger, mode, iepoch=0, ibatch=1, totalbatch=1, meterlist=None, prefix=True):
    pstr = "{}:\t[{}][{}/{}] ".format(mode, iepoch, ibatch, totalbatch) if prefix else ""
    space = len(pstr)

    if meterlist is None:
        meterlist = ["loss"]+[m for m in logger.meter.keys() if m != "loss"]
    for i, meter in enumerate(meterlist):
        if i > 0:
            pstr += " "*space

        if meter in ['confusion', 'histogram', 'image']:
            continue
        if meter == 'accuracy':
            pstr += "Acc@1 {:.2f}% \t Acc@{} {:.2f}% \n".format(logger.meter[meter].value()[0], logger.topk, logger.meter[meter].value()[1])
        elif meter == 'map':
            pstr += "mAP {:.3f} \n".format(logger.meter[meter].value())
        elif meter == 'auc':
            pstr += "AUC {:.3f} \n".format(logger.meter[meter].value())
        else:
            pstr += "{} {:.3f} ({:.3f}) \n".format(meter, logger.meter[meter].val, logger.meter[meter].mean)

    return pstr

def graph_logger(logger, phase, epoch, final=False, meterlist=None):
    if meterlist is None:
        meterlist = logger.meter.keys()

    if not final:
        stats_fname = "Sparse3DUnet_statistics_{}.tsv".format(phase)
        with open(stats_fname, "a") as stats_file:
            print ("epoch {} {}".format(epoch, format_meter(logger, phase)), file=stats_file)

    for meter_name in meterlist:
        if meter_name in ['confusion', 'histogram', 'image']:
            continue
        pp = PdfPages('{}_epoch{}_{}.pdf'.format(phase if not final else "final", epoch, meter_name))
        f, ax = fig, ax = plt.subplots(figsize=(8.5,11))
        if final:
            f.suptitle("Sparse 3D Unet {} Final - {}".format(phase.title(), meter_name), fontsize=14)
        else:
            f.suptitle("Sparse 3D Unet {} Epoch {} - {}".format(phase.title(), epoch, meter_name), fontsize=14)
        ax.set_xlabel("Batch #" if not final else "Epoch #")
        ax.set_ylabel("")
        ax.plot(values[meter_name] if not final else [epoch_values[meter_name][epoch] \
            for epoch in sorted(epoch_values[meter_name].keys())])
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)
        plt.clf()
        plt.cla()
        del f
        del ax

    global values
    del values
    values = {}

class ModelStats(object):
    all_dice = {"train":[], "val":[]}
    all_accuracies = {"train":[], "val":[]}
    all_fpr = {"train":[], "val":[]}
    all_tpr = {"train":[], "val":[]}

    def __init__(self, phase):
        self.phase = phase
        self.top1 = 0
        self.top5 = 0
        self.n = 0
        self.nll = 0
        self.running_corrects_num = 0 #TP+TN
        self.fpr = []
        self.tpr = []
        self.mcc = []
        self.precision = []
        self.kappa = []
        self.f1 = []

        self.accuracies = []
        self.losses = []

    def update(self, output, target, loss, locations, epoch, since, loss_function="dice", smooth=1.0):
        batchSize = output.size(0)
        #import pdb; pdb.set_trace()



        self.n += batchSize
        # self.nll += loss * batchSize
        # _target = target.cpu()
        predicted_corrects_raw = output>0.7
        predicted_corrects = predicted_corrects_raw.float()
        predicted_corrects_num = predicted_corrects.eq(target).sum().data[0]
        self.running_corrects_num += predicted_corrects_num #TP+TN
        predicted_corrects_num /= batchSize

        iflat = predicted_corrects.view(-1)
        tflat = target.view(-1)
        intersection = (iflat * tflat).sum()

        dice = ((2. * intersection + smooth) / ((iflat.sum() + tflat.sum() + smooth))).data[0]
        mcc = metrics.matthews_corrcoef(target.view(-1).data.numpy(), predicted_corrects.view(-1).data.numpy())
        precision = metrics.precision_score(target.view(-1).data.numpy(), predicted_corrects.view(-1).data.numpy())
        kappa = metrics.cohen_kappa_score(target.view(-1).data.numpy(), predicted_corrects.view(-1).data.numpy())
        f1 = metrics.f1_score(target.view(-1).data.numpy(), predicted_corrects.view(-1).data.numpy())

        del predicted_corrects_raw
        del predicted_corrects

        previous_row = 0
        predicted_corrects_batch = 0
        mcc_batch = 0
        precision_batch = 0
        kappa_batch = 0
        f1_batch= 0
        dice_batch = 0

        samples = locations[:, 3]
        num_samples = float(samples[-1]+1)

        for i, sample in groupby(enumerate(samples), key=lambda x:x[1]):
            for voxel_end in sample: pass
            batch_predictions = output[previous_row:voxel_end[0]+1]
            target_values = target[previous_row:voxel_end[0]+1]


            iflat = (batch_predictions>0.7).float().view(-1)
            tflat = target_values.view(-1)
            intersection = (iflat * tflat).sum()

            dice_val = ((2. * intersection + smooth) / ((iflat.sum() + tflat.sum() + smooth))).data[0]
            dice_batch += dice_val

            sample_size = float(voxel_end[0]-previous_row+1)

            #predicted_corrects_raw_sample = batch_predictions>=0.8
            predicted_corrects_sample = batch_predictions.cpu().float()

            print "target", tflat.sum().data[0], "of", sample_size
            print "output", iflat.sum().data[0], "of", sample_size
            predicted_corrects_num_sample = iflat.eq(tflat)
            print "equal", predicted_corrects_num_sample.sum().data[0], "of", sample_size
            #predicted_corrects_num_sample = predicted_corrects_num_sample
            #print "sum", predicted_corrects_num_sample
            #predicted_corrects_num_sample /= sample_size
            predicted_corrects_batch += predicted_corrects_num_sample.sum().data[0]

            prediction = (predicted_corrects_sample>0.7).int().cpu().data.numpy()
            mcc_batch += metrics.matthews_corrcoef(tflat.data.numpy().astype(int), prediction)
            precision_batch += metrics.precision_score(tflat.data.numpy(), prediction)
            kappa_batch += metrics.cohen_kappa_score(tflat.data.numpy(), prediction)
            f1_batch += metrics.f1_score(tflat.data.numpy(), prediction)

            previous_row = voxel_end[0]
            #del predicted_corrects_raw_sample
            del predicted_corrects_sample
            del predicted_corrects_num_sample

        mcc_batch /= num_samples
        precision_batch /= num_samples
        kappa_batch /= num_samples
        f1_batch /= num_samples
        dice_batch /= num_samples
        predicted_corrects_batch /= float(batchSize)

        print """Epoch {} {}:
    loss:      {:.4f}
    dice:      {:.4f}%    {:.4f}%
    accuracy:  {:.2f}%    {:.2f}%
    mcc:       {:.4f}%    {:.4f}%
    precision: {:.4f}%    {:.4f}%
    kappa:     {:.4f}%    {:.4f}%
    f1:        {:.4f}%    {:.4f}%
    time:      {:.1f}s""".format(
            epoch,
            self.phase,
            loss,
            dice*100, dice_batch*100,
            predicted_corrects_num*100,  predicted_corrects_batch*100,
            mcc*100, mcc_batch*100,
            precision*100, precision_batch*100,
            kappa*100, kappa_batch*100,
            f1*100, f1_batch*100,
            time.time() - since)

        # try:
        #     fpr, tpr, _ = metrics.roc_curve(_target.view(-1).numpy(), predicted_corrects.view(-1).numpy(), pos_label=1.)
        #     self.tpr += tpr.tolist()
        #     self.fpr += fpr.tolist()
        #     ModelStats.all_fpr[phase] += self.fpr
        #     ModelStats.all_tpr[phase] += self.tpr
        # except (KeyboardInterrupt, SystemExit):
        #     raise
        # except:
        #     pass



        # self.mcc += mcc
        # self.precision += precision
        # self.kappa += kappa
        # self.f1 += f1
        # self.n += n
        # print "Epoch {} {}: corrects:{}% dice:{}% mcc:{}% precision:{}% kappa:{}% f1:{}% time:{}s".format(
        #     epoch, self.phase, accuracy*100, dice*100, mcc*100, precision*100, kappa*100, f1*100, time.time() - since)

        # print "Epoch {} {}: accuracy:{:.2f}% avg-acc:{:.2f}% dice:{:.4f}% mcc:{:.4f}% precision:{:.4f}% kappa:{:.4f}% f1:{:.4f}% time:{:.1f}s".format(
        #     epoch, self.phase, accuracy*100, running_corrects_num*100, dice*100, mcc*100, precision*100, kappa*100, f1*100, time.time() - since)

        #self.accuracies.append(predicted_corrects_num/float(batchSize))
        #self.losses.append(loss)
        #del predicted_corrects_raw
        #del predicted_corrects
        #del running_corrects
        #del true_corrects
        #del _target

    def save(self, phase, epoch):
        ModelStats.all_dice[phase] += self.losses
        ModelStats.all_accuracies[phase].append(self.running_corrects_num/self.n)
        self.plot(phase, epoch)

    def top1pct(self):
        return 100 * (1 - 1.0 * self.top1 / self.n)

    def top5pct(self):
        return 100 * (1 - 1.0 * self.top5 / self.n)

    def nllpct(self):
        return 100*self.nll/self.n

    def correctpct(self):
        return 100*self.running_corrects_num/self.n

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
            fpr, tpr = self.fpr, self.tpr
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
        f, axes = plt.subplots(1, 1, figsize=(8.5,11))
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
        self.plot_roc(axes[2], phase="train", final=True)
        f.subplots_adjust(wspace=.4)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)

        pp = PdfPages('final_validation_statistics.pdf')
        f, axes = plt.subplots(1, 2)
        f.suptitle("Final Validation Statistics", fontsize=14)
        self.plot_accuracy(axes[0], final_phase="val")
        self.plot_loss(axes[1], final_phase="val")
        self.plot_roc(axes[2], phase="val", final=True)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)
