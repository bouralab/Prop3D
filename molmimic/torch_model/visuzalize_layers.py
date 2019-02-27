"""Originally by @author: Utku Ozbulak - github.com/utkuozbulak
"""

import os
import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/data/draizene/molmimic/molmimic/visualize")

import matplotlib
matplotlib.use("Agg")
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)

import argparse
from itertools import product

import numpy as np
import torch
from torch.optim import SGD

import sparseconvnet as scn

import molmimic.torch_model.torch_loader as loader
import mayavi_vieiwer as mv

from molmimic.torch_model import torch_infer as infer

class Sparse3DCNNLayerVisualization():
    """
        Produces an image that minimizes the loss of a convolution
        operation for a specific layer and filter
    """
    def __init__(self, model, selected_layer, selected_filter, use_gpu=True):
        self.model = model
        self.model.eval()
        self.selected_layer = selected_layer
        self.selected_filter = selected_filter
        self.conv_output = 0

        # Generate a 3D random volume, fills entire volume, but uses sparse matrices
        dim = np.arange(0, 96)
        x, y, z = np.meshgrid(dim, dim, dim)
        points = list(zip(x.ravel(), y.ravel(), z.ravel()))
        features = np.array(list(product([0,1], repeat=3)))
        features = features[np.random.choice(8, 96*96*96)]

        self.inputSpatialSize = torch.LongTensor((96,96,96))
        self.created_image = scn.InputBatch(3, self.inputSpatialSize)
        self.created_image.addSample()
        indices = torch.LongTensor(points)
        labels = torch.from_numpy(features).float()
        self.created_image.setLocations(indices, labels, 0)

        if use_gpu:
            self.created_image = self.created_image.cuda()

        self.created_image = self.created_image.to_variable(requires_grad=True)

        del indices
        del labels
        del features
        del points
        del x
        del y
        del z
        del dim

        # Create the folder to export images if not exists
        if not os.path.exists('generated'):
            os.makedirs('generated')

    def hook_layer(self):
        def hook_function(module, grad_in, grad_out):
            # Gets the conv output of the selected filter (from selected layer)
            self.conv_output = grad_out.features[0, self.selected_filter]

        # Hook the selected layer
        self.model.sparseModel[self.selected_layer].register_forward_hook(hook_function)

    def visualise_layer_with_hooks(self):
        # Hook the selected layer
        self.hook_layer()

        # Define optimizer for the image
        # Earlier layers need higher learning rates to visualize whereas later layers need less
        optimizer = SGD([self.created_image.features], lr=5, weight_decay=1e-6)
        for i in range(1, 51):
            optimizer.zero_grad()
            # Assign create image to a variable to move forward in the model
            x = self.created_image
            for index, layer in enumerate(self.model.sparseModel):
                # Forward pass layer by layer
                # x is not used after this point because it is only needed to trigger
                # the forward hook function
                x = layer(x)
                # Only need to forward until the selected layer is reached
                if index == self.selected_layer:
                    # (forward hook function triggered)
                    break

            # Loss function is the mean of the output of the selected layer/filter
            # We try to minimize the mean of the output of that specific filter
            loss = torch.mean(self.conv_output)
            print("Iteration: {}, Loss: {:.2f}".format(i, loss.data.cpu()[0]))

            # Backward
            loss.backward()

            # Update image
            optimizer.step()

            # Save image
            if i % 5 == 0:
                path = "generated/layer_vis_l{}_f{}_iter{}".format(self.selected_layer, self.selected_filter, i)
                colors = (self.created_image.features.data > 0.7).long()
                volume = colors.view(96, 96, 96, -1).cpu().numpy()
                #np.savez(path+".npz", volume)
                fig, axes = mv.create_figure(1, size=(96,96,96), no_prediction=True)
                mv.plot_volume_matplotlib(
                    axes[0],
                    volume)
                mv.save_fig(fig, path)
                del fig
                del axes

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="visualize layers")
    infer.create_args(parser, default_shape=(96,96,96), add_data=False)
    parser.add_argument(
        "--n_features",
        default=3,
        type=int
    )
    parser.add_argument(
        "--dropout-width",
        default=True,
        action="store_true",
        help="Apply dropout after convolution operations on width"
    )
    parser.add_argument(
        "--dropout-depth",
        default=False,
        action="store_true",
        help="Apply dropout after convolution operations on depth"
    )
    parser.add_argument(
        "--dropout-p",
        default=0.7,
        type=float
    )
    parser.add_argument(
        "-l", "--layer",
        type=int,
        default=0
    )
    parser.add_argument(
        "-f", "--filter",
        type=int,
        default=0
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    model = infer.load_model(args.model, nFeatures=args.n_features, use_resnet_unet=True, dropout_depth=args.dropout_depth, dropout_width=args.dropout_width, dropout_p=args.dropout_p)
    layer_vis = Sparse3DCNNLayerVisualization(model, args.layer, args.filter)
    layer_vis.visualise_layer_with_hooks()
