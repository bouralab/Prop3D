import os, sys
sys.path.append("/data/draizene/molmimic")

import numpy as np
from molmimic.calculate_features import SwarmJob

def vary_learning_rate(name, dataset, ibis_data, start, end, step_pct, output_dir=None):
    job = SwarmJob("{}_ppi_lr".format(name), cpus=14, gpus=1, merge_output=True, threads_per_job=1)
    if output_dir is None:
        output_dir = os.getcwd()
    step = -1 if start > end else 0
    for lr in range(start, end, step):
        print(lr)
        step = 10**((-1*lr)+1)*step_pct
        print(step)
        lr_range = np.arange(0, 10**((-1*lr)+1), step)
        lr_range[0] = 10**(-1*lr)
        for _lr in lr_range:
            print(_lr)
            lr_out_dir = os.path.join(output_dir, str(_lr))
            if not os.path.exists(lr_out_dir):
                os.makedirs(lr_out_dir)
            job += "cd {}; ".format(lr_out_dir)
            job += "/data/draizene/3dcnn-torch-py2 python /data/draizene/molmimic/molmimic/torch_model/torch_train.py "
            job += "--use-resnet-unet --epochs 100 --dropout-width --dropout-p 0.7 "
            job += "--learning-rate {} ".format(_lr)
            job += "{} {} \n".format(dataset, ibis_data)
    print(job.run())

if __name__ == "__main__":
    if len(sys.argv) == 7:
        name, start, end, step_pct, dataset, ibis_data = sys.argv[1:]
        vary_learning_rate(name, dataset, ibis_data, int(start), int(end), float(step_pct))
    else:
        print("Invalid arguments: name, start, end, step, dataset, ibis_data")
