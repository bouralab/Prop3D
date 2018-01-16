import torch
import torch.nn as nn
import sparseconvnet as scn

class UNet3D(nn.Module):
    def __init__(self, in_channel, n_classes, batchnorm=True):
        self.in_channel = in_channel
        self.n_classes = n_classes
        super(UNet3D, self).__init__()

        self.conv1_1 = self.encoder(in_channel, 32, bias=False, batchnorm=batchnorm)
        self.conv1_2 = self.encoder(32, 64, bias=False, batchnorm=batchnorm)
        self.pool1 = scn.MaxPooling(3, 2, 2)

        self.conv2_1 = self.encoder(64, 64, bias=False, batchnorm=batchnorm)
        self.conv2_2 = self.encoder(64, 128, bias=False, batchnorm=batchnorm)
        self.pool2 = scn.MaxPooling(3, 2, 2)

        self.conv3_1 = self.encoder(128, 128, bias=False, filter_stride=1, filter_size=3, batchnorm=batchnorm)
        self.conv3_2 = self.encoder(128, 256, bias=False, filter_stride=1, filter_size=3, batchnorm=batchnorm)
        self.pool3 = scn.MaxPooling(3, 2, 2)

        self.conv4_1 = self.encoder(256, 256, bias=False, batchnorm=batchnorm)
        self.conv4_2 = self.encoder(256, 512, bias=False, batchnorm=batchnorm)

        self.up5_1 = self.decoder(512, 512, filter_size=2, filter_stride=2, bias=False)
        self.up5_2 = scn.JoinTable()
        self.conv5_1 = self.encoder(256+512, 256, bias=False, batchnorm=batchnorm)
        self.conv5_2 = self.encoder(256, 256, bias=False, batchnorm=batchnorm)

        self.up6_1 = self.decoder(256, 256, filter_size=2, filter_stride=2, bias=False)
        self.up6_2 = scn.JoinTable()
        self.conv6_1 = self.encoder(128 + 256, 128, bias=False)
        self.conv6_2 = self.encoder(128, 128, bias=False)

        self.up7_1 = self.decoder(128, 128, filter_size=2, filter_stride=2, bias=False)
        self.up7_2 = scn.JoinTable()
        self.conv7_1 = self.encoder(64 + 128, 64, bias=False, batchnorm=batchnorm)
        self.conv7_2 = self.encoder(64, 64, bias=False, batchnorm=batchnorm)

        self.conv8 = self.encoder(64, n_classes, filter_size=1, bias=False, batchnorm=batchnorm)
        self.act = scn.Sigmoid()

        self.log_level = 0

    def input_spatial_size(self, out_size):
        return out_size

    def set_log_level(self, level=None):
        self.log_level = level or 0

    def encoder(self, in_channels, out_channels, filter_size=3, filter_stride=1, bias=True, batchnorm=True, submanifold=True):
        layer = scn.Sequential(
            scn.SubmanifoldConvolution(3, in_channels, out_channels, filter_size, bias) if submanifold \
                else scn.Convolution(3, in_channels, out_channels, filter_size, filter_stride, bias),
            scn.BatchNormReLU(out_channels) if batchnorm else scn.ReLU())
        return layer

    def decoder(self, in_channels, out_channels, filter_size, filter_stride=1, bias=True):
        layer = scn.Sequential(
            scn.Deconvolution(3, in_channels, out_channels, filter_size, filter_stride, bias),
            scn.ReLU())
        return layer

    def forward(self, x):
        verbose = self.log_level
        if verbose > 0: print "input", x, x.spatial_size.tolist(), x.features.size()
        if verbose > 1: print "   ", x.features.view(-1).cpu().data.numpy().tolist()
        conv1 = self.conv1_1(x)
        if verbose > 0: print "conv1_1", conv1, conv1.spatial_size.tolist(), conv1.features, conv1.features.size()
        if verbose > 1: print "   ", conv1.features
        conv1 = self.conv1_2(conv1)
        if verbose > 0: print "conv1_2", conv1, conv1.spatial_size.tolist(), conv1.features, conv1.features.size()
        if verbose > 1: print "   ", conv1.features
        pool1 = self.pool1(conv1)
        if verbose > 0: print "pool1", pool1.spatial_size.tolist(), pool1.features.size()
        if verbose > 1: print "   ", pool1.features

        conv2 = self.conv2_1(pool1)
        if verbose > 0: print "conv2_1", conv2.spatial_size.tolist(), conv2.features.size()
        if verbose > 1: print "   ", conv2.features
        conv2 = self.conv2_2(conv2)
        if verbose > 0: print "conv2_2", conv2.spatial_size.tolist(), conv2.features.size()
        if verbose > 1: print "   ", conv2.features
        pool2 = self.pool2(conv2)
        if verbose > 0: print "pool2", pool2.spatial_size.tolist(), pool2.features.size()
        if verbose > 1: print "   ", pool2.features

        conv3 = self.conv3_1(pool2)
        if verbose > 0: print "conv3_1", conv3.spatial_size.tolist(), conv3.features.size()
        if verbose > 1: print "   ", conv3.features
        conv3 = self.conv3_2(conv3)
        if verbose > 0: print "conv3_3", conv3.spatial_size.tolist(), conv3.features.size()
        if verbose > 1: print "   ", conv3.features
        pool3 = self.pool3(conv3)
        if verbose > 0: print "pool3", pool3.spatial_size.tolist(), pool3.features.size()
        if verbose > 1: print "   ", pool3.features

        conv4 = self.conv4_1(pool3)
        if verbose > 0: print "conv4_1", conv4.spatial_size.tolist(), conv4.features.size()
        if verbose > 1: print "   ", conv4.features
        conv4 = self.conv4_2(conv4)
        if verbose > 0: print "conv4_2", conv4.spatial_size.tolist(), conv4.features.size()
        if verbose > 1: print "   ", conv4.features

        up5 = self.up5_1(conv4)
        if verbose > 0: print "up5_1", up5.spatial_size.tolist(), up5.features.size()
        if verbose > 1: print "   ", up5.features
        up5 = self.up5_2((up5, conv3))
        if verbose > 0: print "up5_2", up5.spatial_size.tolist(), up5.features.size()
        if verbose > 1: print "   ", up5.features
        conv5 = self.conv5_1(up5)
        if verbose > 0: print "conv5_1", conv5.spatial_size.tolist(), conv5.features.size()
        if verbose > 1: print "   ", conv5.features
        conv5 = self.conv5_2(conv5)
        if verbose > 0: print "conv5_2", conv5.spatial_size.tolist(), conv5.features.size()
        if verbose > 1: print "   ", conv5.features

        up6 = self.up6_1(conv5)
        if verbose > 0: print "up6_1", up6.spatial_size.tolist(), up6.features.size()
        if verbose > 1: print "   ", up6.features
        up6 = self.up6_2((up6, conv2))
        if verbose > 0: print "up6_2", up6.spatial_size.tolist(), up6.features.size()
        if verbose > 1: print "   ", up6.features
        conv6 = self.conv6_1(up6)
        if verbose > 0: print "conv6_1", conv6.spatial_size.tolist(), conv6.features.size()
        if verbose > 1: print "   ", conv6.features
        conv6 = self.conv6_2(conv6)
        if verbose > 0: print "conv6_2", conv6.spatial_size.tolist(), conv6.features.size()
        if verbose > 1: print "   ", conv6.features

        up7 = self.up7_1(conv6)
        if verbose > 0: print "up7_1", up7.spatial_size.tolist(), up7.features.size()
        if verbose > 1: print "   ", up7.features
        up7 = self.up7_2((up7, conv1))
        if verbose > 0: print "up7_2", up7.spatial_size.tolist(), up7.features.size()
        if verbose > 1: print "   ", up7.features
        conv7 = self.conv7_1(up7)
        if verbose > 0: print "conv7_1", conv7.spatial_size.tolist(), conv7.features.size()
        if verbose > 1: print "   ", conv7.features
        conv7 = self.conv7_2(conv7)
        if verbose > 0: print "conv7_2", conv7.spatial_size.tolist(), conv7.features.size()
        if verbose > 1: print "   ", conv7.features

        conv8 = self.conv8(conv7)
        if verbose > 0: print "conv8", conv8.spatial_size.tolist(), conv8.features.size()
        if verbose > 1: print "   ", conv8.features
        act = self.act(conv8)

        return act
