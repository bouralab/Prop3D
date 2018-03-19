from itertools import groupby
from torch.autograd import Function
from torch.nn.modules.loss import _Loss

class DiceLoss(_Loss):
    def __init__(self, size_average=True, smooth=1.):
        super(DiceLoss, self).__init__(size_average)
        self.smooth = smooth

    def forward(self, input, target, locations, weight=None):
        if self.size_average:
            return -self.dice_coef_samples(input, target, locations, weight)

        return -self.dice_coef_batch(input, target, weight)

    def dice_coef_batch(self, input, target, weight=None):
        iflat = (input>0.7).float().view(-1)
        tflat = target.view(-1)
        intersection = (iflat * tflat).sum()

        #Do per batch
        dice = ((2. * intersection + self.smooth) / ((iflat.sum() + tflat.sum() + self.smooth)))

        if weight is not None:
            dice *= weight

        return dice

    def dice_coef_samples(self, input, target, locations, weight=None):
        total_size = input.size(0)
        samples = locations[:, 3]
        num_samples = samples[-1]+1
        previous_row = 0
        dice = None
        total = input.size(0)

        if weight is not None:
            use_sample_weights = weight.shape[0]>1
            if use_sample_weights:
                assert use_sample_weights and weight.shape[0] == num_samples, "{} {}".format(weight.shape[0], num_samples)
                total = weight.sum()
            else:
                total = weight.data[0]

        for i, sample in groupby(enumerate(samples), key=lambda x:x[1]):
            for voxel_end in sample: pass

            print i, ":", previous_row, "to", voxel_end[0], "of", total_size
            batch_predictions = input[previous_row:voxel_end[0]+1]
            target_values = target[previous_row:voxel_end[0]+1]
            previous_row = voxel_end[0]

            iflat = (batch_predictions > 0.7).float().view(-1)
            tflat = target_values.view(-1)
            intersection = (iflat * tflat).sum()

            dice_val = ((2. * intersection + self.smooth) / ((iflat.sum() + tflat.sum() + self.smooth)))
            print dice_val

            if use_sample_weights:
                dice_val *= weight[i]

            if dice is None:
                dice = dice_val
            else:
                dice += dice_val

        if weight is not None and not use_sample_weights:
            dice *= weight

        return dice/total if total > 0 else self.smooth

class IoULoss(_Loss):
    """A intersect B / A union B = A intersect B / (area(A) + area(B) - A intersect B)
    """
    def __init__(self, size_average=True, smooth=1.):
        super(IoULoss, self).__init__(size_average)
        self.smooth = smooth

    def forward(self, input, target, weight=None):
        if self.size_average:
            return -self.dice_coef_samples(input, target, locations, weights)

        return -self.IoU_batch(input, target, weight=None)

    def IoU_batch(self, input, target):
        iflat = input.view(-1)
        tflat = target.view(-1)
        intersection = (iflat * tflat).sum()

        return ((intersection + self.smooth) / ((iflat.sum() + tflat.sum() - intersection + self.smooth)))

    def IoU_samples(self, input, target, locations, weight=None):
        samples = locations[:, 3]
        num_samples = samples[-1]+1
        previous_row = 0
        IoU = None

        if weight is not None and isinstance(weight, (list, tuple)):
            assert len(weight) == num_samples
            use_sample_weights = True

        for i, sample in groupby(enumerate(samples), key=lambda x:x[1]):
            for voxel_end in sample: pass

            batch_predictions = input[previous_row:voxel_end[1]+1]
            target_values = target[previous_row:voxel_end[1]+1]
            previous_row = voxel_end[1]

            iflat = batch_predictions.view(-1)
            tflat = target_values.view(-1)
            intersection = (iflat * tflat).sum()

            IoU_val = ((intersection + self.smooth) / ((iflat.sum() + tflat.sum() - intersection + self.smooth)))

            if dice is None:
                IoU = IoU_val
            else:
                IoU += IoU_val

        if weight is not None and not use_sample_weights:
            IoU *= weights

        return IoU/num_samples

def dice_loss(input, target, smooth=1.0, weight=None):
    return DiceLoss(smooth=smooth)(input, target)
    iflat = input.view(-1)
    tflat = target.view(-1)
    intersection = (iflat * tflat).sum()

    #Do per batch
    dice = ((2.0 * intersection + smooth) / ((iflat.sum() + tflat.sum() + smooth)))

    if weight is not None:
        dice *= weight

    return -1.0*dice

class DiceLossTorchBioMed(Function):
    """From torchbiomed: https://github.com/mattmacy/torchbiomed"""
    def __init__(self, *args, **kwargs):
        pass

    def forward(self, _input, _target, save=True):

        input = _input.view(-1)
        target = _target.view(-1)

        if save:
            self.save_for_backward(input, target)

        eps = 0.000001
        _, result_ = input.max(0)
        result_ = torch.squeeze(result_)
        if input.is_cuda:
            result = torch.cuda.FloatTensor(result_.size())
            self.target_ = torch.cuda.FloatTensor(target.size())
        else:
            result = torch.FloatTensor(result_.size())
            self.target_ = torch.FloatTensor(target.size())
        result.copy_(result_)
        self.target_.copy_(target)
        target = self.target_
        intersect = torch.dot(result, target)
        # binary values so sum the same as sum of squares
        result_sum = torch.sum(result)
        target_sum = torch.sum(target)
        union = result_sum + target_sum + (2*eps)

        # the target volume can be empty - so we still want to
        # end up with a score of 1 if the result is 0/0
        IoU = intersect / union
        print('union: {:.3f}\t intersect: {:.6f}\t target_sum: {:.0f} IoU: result_sum: {:.0f} IoU {:.7f}'.format(
            union, intersect, target_sum, result_sum, 2*IoU))
        out = torch.FloatTensor(1).fill_(2*IoU)
        self.intersect, self.union = intersect, union
        return out

    def backward(self, grad_output):
        input, _ = self.saved_tensors
        intersect, union = self.intersect, self.union
        target = self.target_
        gt = torch.div(target, union)
        IoU2 = intersect/(union*union)
        pred = torch.mul(input[:, 1], IoU2)
        dDice = torch.add(torch.mul(gt, 2), torch.mul(pred, -4))
        grad_input = torch.cat((torch.mul(dDice, -grad_output[0]),
                                torch.mul(dDice, grad_output[0])), 0)
        return grad_input , None

class WeightedDiceLossVNet(_Loss):
    """
    https://mattmacy.github.io/vnet.pytorch
    Not functional. Port from Lasagna in progress.
    """
    def __init__(self, size_average=True, smooth=1.):
        super(WeightedDiceLossVNet, self).__init__(size_average)
        self.smooth = smooth

    def forward(self, input, target):
        n = input.size(0)
        dice = torch.FloatTensor(n).zero_()
        self.union = torch.FloatTensor(n).zero_()
        self.intersection = torch.FloatTensor(n).zero_()

        self.result = np.reshape(np.squeeze(np.argmax(bottom[0].data[...],axis=1)),[bottom[0].data.shape[0],bottom[0].data.shape[2]])
        self.gt = np.reshape(np.squeeze(bottom[1].data[...]),[bottom[1].data.shape[0],bottom[1].data.shape[2]])

        self.gt = (self.gt > 0.5).astype(dtype=np.float32)
        self.result = self.result.astype(dtype=np.float32)

        for i in xrange(0,n):
            # compute dice
            CurrResult = (self.result[i,:]).astype(dtype=np.float32)
            CurrGT = (self.gt[i,:]).astype(dtype=np.float32)

            self.union[i] = torch.sum(CurrResult)+torch.sum(CurrGT)
            self.intersection[i] = torch.sum(CurrResult * CurrGT)

            dice[i] = 2 * self.intersection[i] / (self.union[i]+0.00001)
            print dice[i]

        top[0].data[0]=np.sum(dice)

    def backward(self, top, propagate_down, bottom):
        for btm in [0]:
            prob = bottom[0].data[...]
            bottom[btm].diff[...] = np.zeros(bottom[btm].diff.shape, dtype=np.float32)
            for i in range(0, bottom[btm].diff.shape[0]):

                bottom[btm].diff[i, 0, :] += 2.0 * (
                    (self.gt[i, :] * self.union[i]) / ((self.union[i]) ** 2) - 2.0*prob[i,1,:]*(self.intersection[i]) / (
                    (self.union[i]) ** 2))
                bottom[btm].diff[i, 1, :] -= 2.0 * (
                    (self.gt[i, :] * self.union[i]) / ((self.union[i]) ** 2) - 2.0*prob[i,1,:]*(self.intersection[i]) / (
                    (self.union[i]) ** 2))
