import sys
sys.path.append("/data/draizene/molmimic")
sys.path.append("/data/draizene/molmimic/molmimic/visualize")

import matplotlib
matplotlib.use("Agg")
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)

import argparse
from collections import Counter

import numpy as np
import torch

from Bio.PDB.Selection import unfold_entities

import molmimic.torch_model.torch_loader as loader
import mayavi_vieiwer as mv

from molmimic.torch_model import torch_infer as infer
from molmimic.biopdbtools import Structure, InvalidPDB


def infer_spheres(model, shape=(96,96,96), n_samples=3, start_index=0, index_step=1, n_features=3, combinations=False, bs_feature=None, bs_feature2=None, stripes=False, no_prediction=False, use_gpu=False, ibis_data=None, pymol=True):
    if ibis_data is None:
        data = loader.SphereDataset(shape, nFeatures=n_features, allow_feature_combos=combinations, bs_feature=bs_feature, bs_feature2=bs_feature2, stripes=stripes)
    else:
        shape = (264, 264, 264)
        if bs_feature is None:
            random_features = None
        else:
            random_features = (n_features, combinations, bs_feature, bs_feature2)
        train_test = loader.IBISDataset.get_training_and_validation(ibis_data, input_shape=shape, random_features=random_features, transform=False, use_deepsite_features=True)
        data = train_test["val"]

    if ibis_data is None or (ibis_data is not None and not pymol):
        mpl = True
        fig, axes = mv.create_figure(n_samples, size=shape, no_prediction=no_prediction)
    else:
        mpl = False

    runs = []

    for sample in range(start_index, start_index+n_samples, index_step):
        print("Running", sample, data[sample]["id"])
        sphere = loader.sparse_collate([data[sample]])

        print("center", np.mean(sphere["indices"][0], axis=0))

        if mpl:
            ax = axes[sample]
        truth = np.where(sphere["truth"][0][:, 1]==1.0)
        truth_voxels = sphere["indices"][0][truth]
        if mpl:
            rot_z180, rot_x45 = mv.plot_volume_matplotlib(
                ax,
                sphere["indices"][0],
                colors=sphere["data"][0],
                truth=truth_voxels) #sphere["indices"][0]
            #ax.set_title("Truth", fontdict={"fontsize":20})

        if no_prediction:
            if not mpl:
                view_in_pymol(data[sample]["id"], truth_voxels=truth_voxels)
            continue

        output, logs = infer.infer(model, sphere, input_shape=shape, use_gpu=use_gpu)
        if mpl:
            ax = axes[sample+n_samples]
        _, prediction = (output>=0.7).max(dim=1)
        prediction = prediction.data.cpu().numpy()
        prediction = np.where(prediction==1)[0]
        prediction_voxels = sphere["indices"][0][prediction]
        print(prediction.shape, prediction_voxels.shape)
        colors = np.tile(np.array([0.95, 0.37, 0.18]), (prediction_voxels.size,1))
        if mpl:
            mv.plot_volume_matplotlib(
                ax, prediction_voxels,
                colors=colors,
                truth=truth_voxels,
                rot_z180=rot_z180,
                rot_x45=rot_x45)
            #ax.set_title("Prediction", fontdict={"fontsize":20})
        print(data[sample]["id"])
        print(logs.meter["dice_avg"].val)
        print(logs.meter["dice_class1"].val)
        print(logs.meter["mcc_avg"].val)
        runs.append([data[sample]["id"], logs.meter["dice_avg"].val, logs.meter["dice_class1"].val, logs.meter["mcc_avg"].val])
        print()

        if not mpl:
            view_in_pymol(data[sample]["id"], prediction_voxels, truth_voxels)
    print("Saving Figure")
    print(runs)
    if mpl:
        mv.save_fig(fig, "sphere_infer.pdf")

def atom_volume(structure, atom):
    return (4/3.)*np.pi*structure.get_vdw(atom)**3

def view_in_pymol(id, predicted_voxels=None, truth_voxels=None, voxel_atom_ratio=.2):
    pdb, chain = id.split(".")
    structure = Structure.from_pdb(pdb, chain, rotate=False)

    cmd = """fetch {id}
remove hetatm
hide everything, {id}
show surface, {id}
color gray90, {id}
""".format(id=id)

    if truth_voxels is not None:
        truth_atoms = Counter()
        for v in truth_voxels:
            atoms = structure.convert_voxels(v, level="A")
            if len(atoms) > 0:
                truth_atoms[atoms[0]] += 1

        truth_atoms = [atom for atom, count in list(truth_atoms.items()) \
            if float(count)/atom_volume(structure, atom) >= voxel_atom_ratio]

        truth_residues = [str(r.get_id()[1]) for r in unfold_entities(truth_atoms, "R")]
        truth_resi = "+".join(truth_residues)

        cmd += """select true_binding_site, resi {true_resi}
color orange, true_binding_site
""".format(true_resi=truth_resi)

    if predicted_voxels is not None:
        predicted_atoms = Counter()
        for v in predicted_voxels:
            atoms = structure.convert_voxels(v, level="A")
            if len(atoms) > 0:
                predicted_atoms[atoms[0]] += 1

        predicted_atoms = [atom for atom, count in list(predicted_atoms.items()) \
            if float(count)/atom_volume(structure, atom) >= voxel_atom_ratio]

        predicted_residues = [str(r.get_id()[1]) for r in unfold_entities(predicted_atoms, "R")]
        predicted_resi = "+".join(truth_residues)

        cmd += """select predicted_binding_site, resi {predicted_resi}
color magenta, predicted_binding_site
""".format(predicted_resi=predicted_resi)

    if truth_voxels is not None and predicted_voxels is not None:
        false_postive_voxels = set(predicted_residues)-set(truth_residues)
        fp_resi = "+".join(false_postive_voxels)
        cmd += """select false_positive_binding_site, resi {fp_resi}
color blue, false_positive_binding_site
""".format(fp_resi=fp_resi)

    with open("{}_pymol.cmd".format(id), "w") as f:
        print(cmd, file=f)

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Infer spheres")
    infer.create_args(parser, default_shape=(96,96,96), add_data=False)
    parser.add_argument(
        "--n_samples",
        default=3,
        type=int
    )
    parser.add_argument(
        "--n_features",
        default=3,
        type=int
    )
    parser.add_argument(
        "--combination",
        default=False,
        action="store_true"
    )
    parser.add_argument(
        "--dropout-width",
        default=False,
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
        default=0.5,
        type=float
    )
    parser.add_argument(
        "--bs-feature",
        default=None
    )
    parser.add_argument(
        "--bs-feature2",
        default=None
    )
    parser.add_argument(
        "--stripes",
        default=False,
        action="store_true",
        help="On spherical models, apply bs-feature2 as stripes on patch"
    )
    parser.add_argument(
        "--no-prediction",
        default=False,
        action="store_true"
    )
    parser.add_argument(
        "--ibis-data",
        default=None
    )
    parser.add_argument(
        "--start-index",
        default=0,
        type=int
    )
    parser.add_argument(
        "--index-step",
        default=1,
        type=int
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    if not args.no_prediction:
        model = infer.load_model(args.model, nFeatures=args.n_features, no_batch_norm=args.no_batch_norm, use_resnet_unet=True, dropout_depth=args.dropout_depth, dropout_width=args.dropout_width, dropout_p=args.dropout_p)
    else:
        model = None

    infer_spheres(model, shape=args.shape, n_samples=args.n_samples, start_index=args.start_index, index_step=args.index_step, n_features=args.n_features, combinations=args.combination, bs_feature=args.bs_feature, bs_feature2=args.bs_feature2, stripes=args.stripes, no_prediction=args.no_prediction, use_gpu=True, ibis_data=args.ibis_data)
