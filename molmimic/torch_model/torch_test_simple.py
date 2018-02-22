import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/data/draizene/molmimic/molmimic/visualize")

import matplotlib
matplotlib.use("Agg")
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)

import argparse

import numpy as np
import torch

import molmimic.torch_model.torch_loader as loader
import mayavi_vieiwer as mv

from molmimic.torch_model import torch_infer as infer


def infer_spheres(model, shape=(96,96,96), n_samples=3, use_gpu=False):
    data = loader.SphereDataset(shape)

    fig, axes = mv.create_figure(n_samples, size=shape)

    for sample in xrange(n_samples):
        sphere = loader.sparse_collate([data[sample]])
        output, logs = infer.infer(model, sphere, input_shape=shape, use_gpu=use_gpu)

        ax = axes[sample]
        truth = np.where(sphere["truth"][0][:, 1]==1.0)
        truth_voxels = sphere["indices"][0][truth] 
        rot_z180, rot_x45 = mv.plot_volume_matplotlib(
            ax, 
            sphere["indices"][0], 
            colors=sphere["data"][0], 
            truth=truth_voxels)
        #ax.set_title("Truth", fontdict={"fontsize":20})

        ax = axes[sample+n_samples]
        _, prediction = (output>=0.7).max(dim=1)
        prediction = prediction.data.cpu().numpy()
        prediction = np.where(prediction==1)[0]
        prediction_voxels = sphere["indices"][0][prediction] 
        print prediction.shape, prediction_voxels.shape
        colors = np.tile(np.array([0.95, 0.37, 0.18]), (prediction_voxels.size,1))
        mv.plot_volume_matplotlib(
            ax, prediction_voxels, 
            colors=colors, 
            truth=truth_voxels, 
            rot_z180=rot_z180, 
            rot_x45=rot_x45)
        #ax.set_title("Prediction", fontdict={"fontsize":20})

        print logs.meter["dice_avg"].val
        print logs.meter["dice_class1"].val
        print logs.meter["mcc_avg"].val
        print
    
    mv.save_fig(fig, "sphere_inder.pdf")

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Infer spheres")
    infer.create_args(parser, default_shape=(96,96,96), add_data=False)
    parser.add_argument(
        "--n_samples",
        default=3,
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    model = infer.load_model(args.model, nFeatures=3, no_batch_norm=args.no_batch_norm, use_resnet_unet=True)
    
    infer_spheres(model, shape=args.shape, n_samples=args.n_samples, use_gpu=True)

