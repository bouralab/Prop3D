import os
import argparse

from keras import backend as K
K.set_image_dim_ordering('tf')

from unet3d.model import unet_model_3d
from unet3d.training import load_old_model, train_model
from unet3d.generator import FileGenerator


def train(data_files, truth_files, input_shape=(144,144,144,50))
	model = unet_model_3d(input_shape=input_shape)

	train, validation, t_gen, v_gen = FileGenerator.get_training_and_validation(input_shape, cnt=5, border=10, cnt=5, batch_size=20, n_samples=500)

	train_model(model=model, 
	            model_file=os.path.abspath("./SphereCNN.h5"), 
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
		"-d",
		"--data",
		nargs="+",
		required=True)
	parser.add_argument(
		"-t",
		"--truth",
		nargs="+",
		required=True)
	parser.add_argument(
		"-s",
		"--shape",
		nargs=4,
		default=(144,144,144,50))

	return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()
	train(args.data, args.truth, args.shape)