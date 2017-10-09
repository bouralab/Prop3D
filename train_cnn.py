import os
import argparse
from datetime import datetime

from tensorflow.python.client import device_lib
from keras import backend as K
K.set_image_dim_ordering('tf')

from keras.utils.io_utils import HDF5Matrix
import h5py

from unet3d.original_generator import
from unet3d.model import unet_model_3d
from unet3d.training import train_model

def get_available_gpus():
    local_device_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']

def train(h5_file, input_shape=(96,96,96,50), batch_size=20, data_split=0.8, num_gpus=None):
	model = unet_model_3d(input_shape=input_shape, num_gpus=num_gpus, batch_size=batch_size)

	with h5py.File(os.path.abspath(h5_file), "r") as f:
		size = f["data"].shape[0]

	# n_training = int(size*data_split)

	# X_train = HDF5Matrix(h5_file, 'data', start=0, end=n_training)
	# y_train = HDF5Matrix(h5_file, 'truth', start=0, end=n_training)

	# # Likewise for the test set
	# X_test = HDF5Matrix(h5_file, 'data', start=n_training, end=size)
	# y_test = HDF5Matrix(h5_file, 'truth', start=n_training, end=size)

	# train_model(model=model, 
	#             model_file=os.path.abspath("./molmimic_{}.h5".format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))), 
	#             X=X_train,
	#             y=y_train,
	#             initial_learning_rate=0.001,
	#             learning_rate_drop=0.6,
	#             learning_rate_epochs=10,
	#             n_epochs=50, 
	#             batch_size=batch_size,
	#             validation_split=data_split,
	#             validation_data=(X_test, y_test))


def parse_args():
	parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")
	parser.add_argument(
		"-s",
		"--shape",
		nargs=4,
		default=(96, 96, 96, 50))
	parser.add_argument(
		"--batch",
		type=int,
		default=20)
	parser.add_argument(
		"--split",
		type=float,
		default=0.8)

	gpus = parser.add_mutually_exclusive_group()
	gpus.add_argument(
		"--num_gpus",
		type=int,
		default=1)
	gpus.add_argument(
		"--all_gpus",
		action="store_true",
		default=False)

	parser.add_argument(
		"h5_file")

	args = parser.parse_args()

	if args.all_gpus:
		args.num_gpus = len(get_available_gpus())

	return args

if __name__ == "__main__":
	args = parse_args()
	train(args.h5_file, args.shape, args.batch, args.split, args.num_gpus)
