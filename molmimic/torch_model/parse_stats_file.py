import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

set2 = sns.color_palette("Set2", 8)
sns.set_palette(set2)

from matplotlib.backends.backend_pdf import PdfPages

colors = [
    "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"
]

def parse_stats_file(stats_file):
    epoch = None
    stats = defaultdict(dict)
    for line in stats_file:
        if line.strip() == "": continue
        if line.startswith("epoch"):
            epoch = int(line.split(" ", 2)[1])
            line = line.split("] ", 1)[1]

        fields = line.strip().split()

        if fields[0].startswith("Acc"):
            stats[epoch][fields[0]] = float(fields[1][:-1])
            stats[epoch][fields[2]] = float(fields[3][:-1])
        elif fields[0].startswith("mAP"):
            stats[epoch]["mAP"] = float(fields[1])
        else:
            try:
                stats[epoch][fields[0]] = float(fields[1])
                stats[epoch][fields[0]+"_avg"] = float(fields[2][1:-1])
            except IndexError, ValueError:
                pass

    return stats

def plot_stats(stats_files, names=None, prefix=None, metrics=None, start_epoch=0, end_epoch=None, step_epoch=1):
    stats = [parse_stats_file(f) for f in stats_files]

    if names is not None:
        if isinstance(names, str) or (isinstance(names, (list, tuple)) and len(names) == 1):
            #Format string  with index {i}, {dir}
            name = names[0] if isinstance(names, (list, tuple)) else names
            directory = lambda f: os.path.basename(os.path.dirname(os.path.abspath(f.name)))
            names = [name.format(i=i+1, dir=directory(f)) for i, f in enumerate(stats_files)]
        elif isinstance(names, (list, tuple)) and len(names) != len(stats_files):
            raise RuntimeError("Names must be the same length as stats_files")
    else:
        names = [str(i) for i in xrange(len(stats_files))]

    print names

    if metrics is None:
        metrics = stats[0][0].keys()

    print stats[0][0].keys()
    for stat in stats:
        print stat.keys()

    num_epochs = max(stats[0].keys())+1

    x = range(len(stats_files))
    #use_raw_color = len(names)<=2**len(names[0])
    if end_epoch is None or end_epoch>num_epochs:
        end_epoch = num_epochs
    epochs = range(start_epoch, end_epoch, step_epoch)

    #All Epochs
    for metric in metrics:
        pp = PdfPages('compare_{}_all_epochs.pdf'.format(metric))
        f, ax = fig, ax = plt.subplots(figsize=(9,9))
        f.suptitle("Sparse 3D Unet Compare {}".format(metric.title()), fontsize=14)
        ax.set_xlabel("Epoch #")
        ax.set_ylabel("")
        for i, stat in enumerate(stats):
            #y = [stat[epoch][metric] for epoch in epochs]
            y=[]
            for epoch in epochs:
                print epoch
                print stat[epoch].keys()
                y.append(stat[epoch][metric])
            if False and use_raw_color:
            	ax.plot(y, label=names[i], color=[int(x) for x in names[i]])
            else:
            	ax.plot(y, label=names[i], color=colors[i])
        plt.legend(loc=2, borderaxespad=0.)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)
        plt.clf()
        plt.cla()
        del f
        del ax

    #Final Epoch
    for metric in metrics:
        pp = PdfPages('compare_{}_last_epoch.pdf'.format(metric))
        f, ax = fig, ax = plt.subplots(figsize=(8.5,11))
        f.suptitle("Sparse 3D Unet Compare Final Epoch".format(metric.title()), fontsize=14)
        ax.set_xlabel("Epoch #")
        ax.set_ylabel("")
        x = [stat[end_epoch-1][metric] for i, stat in enumerate(stats)]
        sns.kdeplot(x, shade=True)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(pp, format='pdf')
        pp.close()
        plt.close(f)
        plt.clf()
        plt.cla()
        del f
        del ax

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Parse molmimc stats file and make figures")
    parser.add_argument(
        "--prefix",
        default=None)
    parser.add_argument(
        "-m", "--metrics",
        nargs="*",
        default=None)
    parser.add_argument(
        "-n", "--names",
        nargs="*",
        default=None)
    parser.add_argument(
        "--start-epoch",
        default=0,
        type=int,
        help="Epoch to start parsing (inclusive)"
    )
    parser.add_argument(
        "--end-epoch",
        default=None,
        type=int,
        help="Epoch to stop parsing (exclusive). If None, include last epoch."
    )
    parser.add_argument(
        "--step-epoch",
        default=1,
        type=int
    )
    parser.add_argument(
        "--stats-file",
        nargs="*",
        type=argparse.FileType("r"),
        required=True)


    args = parser.parse_args()

    print args.stats_file

    return args

if __name__ == "__main__":
    args = parse_args()
    plot_stats(args.stats_file, names=args.names, prefix=args.prefix, metrics=args.metrics, start_epoch=args.start_epoch, end_epoch=args.end_epoch, step_epoch=args.step_epoch)
