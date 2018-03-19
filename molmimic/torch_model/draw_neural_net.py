"""
Copyright (c) 2017, Gavin Weiguang Ding
All rights reserved.
Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

input [256, 256, 256]
conv1_1 [256, 256, 256]
conv1_2 [256, 256, 256]
pool1 [128, 128, 128]
conv2_1 [128, 128, 128]
conv2_2 [128, 128, 128]
pool2 [64, 64, 64]
conv3_1 [64, 64, 64]
conv3_3 [64, 64, 64]
pool3 [32, 32, 32]
conv4_1 [32, 32, 32]
conv4_2 [32, 32, 32]
up5_1 [64, 64, 64]
up5_2 [64, 64, 64]
conv5_1 [64, 64, 64]
conv5_2 [64, 64, 64]
up6_1 [128, 128, 128]
up6_2 [128, 128, 128]
conv6_1 [128, 128, 128]
conv6_2 [128, 128, 128]
up7_1 [256, 256, 256]
up7_2 [256, 256, 256]
conv7_1 [256, 256, 256]
conv7_2 [256, 256, 256]
conv8 [256, 256, 256]
"""


import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


NumConvMax = 8
NumFcMax = 20
White = 1.
Light = 0.7
Medium = 0.5
Dark = 0.3
Darker = 0.15
Black = 0.

class DrawCNN(object):
    def __init__(self, model, fc_unit_size=2, layer_width=40, ax=None, fig=None):
        self.model = model
        self.fc_unit_size = fc_unit_size
        self.layer_width = layer_width
        self.layers = []
        self.patches = []

        if ax is None:
            self.fig, self.ax = plt.subplots()
        else:
            self.fig = fig or plt.gcf()
            self.ax = ax

        self.xy = np.array((0., 800.))
        self.loc_diff = np.array([-3, -3])
        self.previous_depth = None

        self.depth_xy = {}
        self.depths = {}

    def save(self, prefix="", show=True):
        plt.tight_layout()
        plt.axis('equal')
        #plt.axis('off')
        self.fig.set_size_inches(12, 12)
        if show:
            plt.show()


        fig.savefig(os.path.join(prefix, "{}_conv_net.png".format(self.model)),
                    bbox_inches='tight', pad_inches=0)

    def _label(self, xy, text, xy_off=[0, 4]):
        plt.text(xy[0] + xy_off[0], xy[1] + xy_off[1], text, family='sans-serif', size=8)

    def add_conv(self, dim, draw_size, depth=1, size=64, in_channels=1, out_channels=32, num=3):
        if self.previous_depth == depth:
            self.xy += np.array((self.layers[-1].get_width()+self.layer_width, 0), dtype=float)
        elif self.previous_depth is not None and self.previous_depth < depth:
            self.xy += np.array((self.layers[-1].get_width()-draw_size, -(self.layers[-1].get_height()+self.layer_width)), dtype=float)
        elif self.previous_depth is not None and self.previous_depth > depth:
            self.xy = np.array((self.layers[-1].get_x()+6, self.depth_xy[depth][1]), dtype=float)
        elif self.previous_depth is None:
            #Start
            self.xy -= np.array((self.layer_width, 0), dtype=float)

        self.depth_xy[depth] = self.xy.copy()
        self.depths[len(self.layers)] = depth
        conv_type = "Inputs" if len(self.layers)==0 else "Channels"
        name = "{}\n{}@{}".format(conv_type, in_channels, "x".join(map(str, [size]*dim)))
        self._label(self.xy+np.array([0,draw_size]), name)
        for ind in xrange(num):
            #print ind * self.loc_diff
            patch = Rectangle(self.xy + ind * self.loc_diff, draw_size, draw_size)
            color = Medium if ind % 2 else Light
            patch.set_color(color * np.ones(3))
            patch.set_edgecolor(Black * np.ones(3))
            self.ax.add_patch(patch)
            if ind == num-1:
                self.layers.append(patch)

        self.previous_depth = depth
        return len(self.layers)-1

         # conv layers
        # size_list = [(32, 32), (18, 18), (10, 10), (6, 6), (4, 4)]
        # num_list = [3, 32, 32, 48, 48]
        # x_diff_list = [0, layer_width, layer_width, layer_width, layer_width]
        # text_list = ['Inputs'] + ['Feature\nmaps'] * (len(size_list) - 1)
        # loc_diff_list = [[3, -3]] * len(size_list)

        # num_show_list = list(map(min, num_list, [NumConvMax] * len(num_list)))
        # top_left_list = np.c_[np.cumsum(x_diff_list), np.zeros(len(x_diff_list))]

        # for ind in range(len(size_list)):
        #     add_layer(patches, colors, size=size_list[ind],
        #               num=num_show_list[ind],
        #               top_left=top_left_list[ind], loc_diff=loc_diff_list[ind])
        #     label(top_left_list[ind], '{}\n{}@{}x{}'.format(text_list[ind],
        #         num_list[ind], size_list[ind][0], size_list[ind][1]))

    def add_conv_patch(self, conv1, conv2, patch_size=5, stride=1, size=3, ratio=None, direction=1):
        if ratio is None:
            ratio = [0.7, 0.5]

        start_loc = self.layers[conv1].get_xy()+np.array((self.layers[conv1].get_width()*ratio[0], self.layers[conv1].get_height()*ratio[1]))
        end_loc = self.layers[conv2].get_xy()+np.array((self.layers[conv2].get_width()*ratio[0], self.layers[conv2].get_height()*ratio[1]))

        # start_loc = top_left_list[ind_bgn] \
        #     + (num_show_list[ind_bgn] - 1) * np.array(loc_diff_list[ind_bgn]) \
        #     + np.array([start_ratio[0] * (size_list[ind_bgn][1] - patch_size[1]),
        #                 - start_ratio[1] * (size_list[ind_bgn][0] - patch_size[0])]
        #                )

        # end_loc = top_left_list[ind_bgn + 1] \
        #     + (num_show_list[ind_bgn + 1] - 1) * np.array(
        #         loc_diff_list[ind_bgn + 1]) \
        #     + np.array([end_ratio[0] * size_list[ind_bgn + 1][1],
        #                 - end_ratio[1] * size_list[ind_bgn + 1][0]])

        patch = Rectangle(start_loc, patch_size, -patch_size)
        patch.set_color(Dark * np.ones(3))
        patch.set_edgecolor(Black * np.ones(3))
        self.ax.add_patch(patch)

        line1 = Line2D([start_loc[0], end_loc[0]], [start_loc[1], end_loc[1]])
        line1.set_color(Darker * np.ones(3))
        self.ax.add_line(line1)

        line2 = Line2D([start_loc[0] + patch_size, end_loc[0]], [start_loc[1], end_loc[1]])
        line2.set_color(Darker * np.ones(3))
        self.ax.add_line(line2)

        line3 = Line2D([start_loc[0], end_loc[0]], [start_loc[1] - patch_size, end_loc[1]])
        line3.set_color(Darker * np.ones(3))
        self.ax.add_line(line3)

        line4 = Line2D([start_loc[0] + patch_size, end_loc[0]], [start_loc[1] - patch_size, end_loc[1]])
        line4.set_color(Darker * np.ones(3))
        self.ax.add_line(line4)


        # in between layers
        # start_ratio_list = [[0.4, 0.5], [0.4, 0.8], [0.4, 0.5], [0.4, 0.8]]
        # end_ratio_list = [[0.4, 0.5], [0.4, 0.8], [0.4, 0.5], [0.4, 0.8]]
        # patch_size_list = [(5, 5), (2, 2), (5, 5), (2, 2)]
        # ind_bgn_list = range(len(patch_size_list))
        # text_list = ['Convolution', 'Max-pooling', 'Convolution', 'Max-pooling']

        # for ind in range(len(patch_size_list)):
        #     add_mapping(
        #         patches, colors, start_ratio_list[ind], end_ratio_list[ind],
        #         patch_size_list[ind], ind,
        #         top_left_list, loc_diff_list, num_show_list, size_list)
        #     label(top_left_list[ind], text_list[ind] + '\n{}x{} kernel'.format(
        #         patch_size_list[ind][0], patch_size_list[ind][1]), xy_off=[26, -65]
        #     )

    def add_max_pool_patch(self, conv1, conv2, pool=2):
        pass

    def add_deconv_patch(self, conv1, conv2, stride=1, size=3):
        pass

    def add_fully_connected(self):
        pass
        # fully connected layers
        # size_list = [(fc_unit_size, fc_unit_size)] * 3
        # num_list = [768, 500, 2]
        # num_show_list = list(map(min, num_list, [NumFcMax] * len(num_list)))
        # x_diff_list = [sum(x_diff_list) + layer_width, layer_width, layer_width]
        # top_left_list = np.c_[np.cumsum(x_diff_list), np.zeros(len(x_diff_list))]
        # loc_diff_list = [[fc_unit_size, -fc_unit_size]] * len(top_left_list)
        # text_list = ['Hidden\nunits'] * (len(size_list) - 1) + ['Outputs']

        # for ind in range(len(size_list)):
        #     add_layer(patches, colors, size=size_list[ind], num=num_show_list[ind],
        #               top_left=top_left_list[ind], loc_diff=loc_diff_list[ind])
        #     label(top_left_list[ind], text_list[ind] + '\n{}'.format(
        #         num_list[ind]))

        # text_list = ['Flatten\n', 'Fully\nconnected', 'Fully\nconnected']

        # for ind in range(len(size_list)):
        #     label(top_left_list[ind], text_list[ind], xy_off=[-10, -65])


if __name__ == '__main__':
    graph = DrawCNN("Molmimic 3DCNN")
    conv1_1 = graph.add_conv(3, 64, depth=0, size=264, in_channels=21, out_channels=32)
    conv1_2 = graph.add_conv(3, 64, depth=0, size=264, in_channels=32, out_channels=64)
    graph.add_conv_patch(conv1_1, conv1_2, stride=1, size=3)
    conv1_3 = graph.add_conv(3, 64, depth=0, size=264, in_channels=32, out_channels=64)
    graph.add_conv_patch(conv1_2, conv1_3, stride=1, size=3, lines=False)

    conv2_1 = graph.add_conv(3, 32, depth=1, size=132, in_channels=64, out_channels=64)
    graph.add_conv_patch(conv1_2, conv2_1, ratio=[0.3, 0.2])
    conv2_2 = graph.add_conv(3, 32, depth=1, size=132, in_channels=64, out_channels=128)
    graph.add_conv_patch(conv2_1, conv2_2, patch_size=3, stride=1, size=3)
    conv2_3 = graph.add_conv(3, 32, depth=1, size=132, in_channels=64, out_channels=128)
    graph.add_conv_patch(conv2_2, conv2_3, patch_size=3, stride=1, size=3, lines=False)

    conv3_1 = graph.add_conv(3, 16, depth=2, size=66, in_channels=128, out_channels=128)
    graph.add_conv_patch(conv2_2, conv3_1, patch_size=3, ratio=[0.3, 0.2])
    conv3_2 = graph.add_conv(3, 16, depth=2, size=66, in_channels=128, out_channels=256)
    graph.add_conv_patch(conv3_1, conv3_2, patch_size=3, stride=1, size=3, ratio=(0.5, 0.7))

    conv4_1 = graph.add_conv(3, 8, depth=3, size=33, in_channels=256, out_channels=256)
    graph.add_conv_patch(conv3_2, conv4_1, patch_size=3, ratio=[0.3, 0.3])
    conv4_2 = graph.add_conv(3, 8, depth=3, size=33, in_channels=256, out_channels=512)
    graph.add_conv_patch(conv4_1, conv4_2, patch_size=2, stride=1, size=3, ratio=[0.5, 0.4])

    conv5_1 = graph.add_conv(3, 16, depth=2, size=66, in_channels=512, out_channels=32)
    graph.add_conv_patch(conv4_2, conv5_1, size=3, patch_size=2, ratio=(0.2, 0.8))
    conv5_2 = graph.add_conv(3, 16, depth=2, size=66, in_channels=512, out_channels=32)
    graph.add_conv_patch(conv5_1, conv5_2, patch_size=3, stride=1, size=3, ratio=(0.5, 0.3))

    conv6_1 = graph.add_conv(3, 32, depth=1, size=132, in_channels=512, out_channels=32)
    graph.add_conv_patch(conv5_2, conv6_1, size=3, patch_size=3, ratio=(0.1, 0.9))
    conv6_2 = graph.add_conv(3, 32, depth=1, size=132, in_channels=512, out_channels=32)
    graph.add_conv_patch(conv6_1, conv6_2, patch_size=3, stride=1, size=3)

    conv7_1 = graph.add_conv(3, 64, depth=0, size=264, in_channels=512, out_channels=32)
    graph.add_conv_patch(conv6_2, conv7_1, size=3, patch_size=3, ratio=(0.1, 0.9))
    conv7_2 = graph.add_conv(3, 64, depth=0, size=264, in_channels=512, out_channels=32)
    graph.add_conv_patch(conv7_1, conv7_2, stride=1, size=3)

    conv8 = graph.add_conv(3, 64, depth=0, size=264, in_channels=64, out_channels=1)
    graph.add_conv_patch(conv7_2, conv8, stride=1, ratio=(0.7, 0.3))

    graph.save()
