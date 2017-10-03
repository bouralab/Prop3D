import os
import argparse
import datetime

from keras import backend as K
K.set_image_dim_ordering('tf')

from unet3d.model import unet_model_3d
from unet3d.training import load_old_model, train_model
from unet3d.generator import DistributedH5FileGenerator


def train(h5_files, input_shape=(144,144,144,50)):
	model = unet_model_3d(input_shape=input_shape)

	train, validation, = DistributedH5FileGenerator.get_training_and_validation(data_file, truth_files, num_classes_per_file=2, num_samples_per_class=20, batch_size=20)

	train_model(model=model, 
	            model_file=os.path.abspath("./molmomic_{}.h5".format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))), 
	            training_generator=t_generator,
	            validation_generator=v_generator, 
	            steps_per_epoch=train.num_steps,
	            validation_steps=validation.num_steps, 
	            initial_learning_rate=0.00001,
	            learning_rate_drop=0.5,
	            learning_rate_epochs=10, 
	            n_epochs=50)

def parse_args():
	parser = argparse.ArgumentParser(description="Load data and truth files to train the 3dCNN")
	parser.add_argument(
		"-s",
		"--shape",
		nargs=4,
		default=(96, 96, 96, 50))
	parser.add_argument(
		"h5_files",
		nargs="+")

	return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()
	train(args.h5_files, args.shape)