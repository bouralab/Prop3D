import sys
sys.path.append("/data/draizene/3DUnetCNN")
sys.path.append("/data/draizene/molmimic")
sys.path.append("/usr/share/pdb2pqr")

import os
import argparse
from datetime import datetime

from tensorflow.python.client import device_lib
from keras import backend as K
K.set_image_dim_ordering('tf')

from keras.utils.io_utils import HDF5Matrix
import h5py

#from unet3d.original_generator import
from unet3d.model import unet_model_3d
from unet3d.training import train_model_generator

from molmimic.keras_model.pdb_generator import IBISGenerator

def get_available_gpus():
    local_device_protos = device_lib.list_local_devices()
    return [x.name for x in local_device_protos if x.device_type == 'GPU']

def train(ibis_data, input_shape=(96,96,96), batch_size=2, data_split=0.8, num_gpus=None, only_aa=False):
	input_shape = input_shape+(21,) if only_aa else input_shape+(59,)
	model = unet_model_3d(input_shape=input_shape, num_gpus=num_gpus)
	if num_gpus is not None and num_gpus > 1:
		batch_size *= num_gpus
	train, validate = IBISGenerator.get_training_and_validation(ibis_data, input_shape=input_shape, batch_size=batch_size, only_aa=only_aa)

	train_model_generator(
		model=model,
        model_file=os.path.abspath("./molmimic_{}.h5".format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))),
        training_generator=train.generate(),
        validation_generator=validate.generate(),
        steps_per_epoch=train.steps_per_epoch,
        validation_steps=validate.steps_per_epoch,
        initial_learning_rate=0.001,
        learning_rate_drop=0.6,
        learning_rate_epochs=10,
        n_epochs=200
		)

	# with h5py.File(os.path.abspath(h5_file), "r") as f:
	# 	size = f["data"].shape[0]

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
		nargs=3,
		default=(96, 96, 96))
	parser.add_argument(
		"--batch",
		type=int,
		default=20)
	parser.add_argument(
		"--split",
		type=float,
		default=0.8)
	parser.add_argument(
		"--only-aa",
		default=False,
		action="store_true",
		help="Only use one feature: aa (20 features since aa is one hot encoded. Else use all 59 features.")

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
		"ibis_data")

	args = parser.parse_args()

	if args.all_gpus:
		args.num_gpus = len(get_available_gpus())

	return args

if __name__ == "__main__":
	args = parse_args()
	train(args.ibis_data, args.shape, args.batch, args.split, args.num_gpus, args.only_aa)
