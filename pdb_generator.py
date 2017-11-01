import pandas as pd
import numpy as np
from biopdbtools import Structure
from unet3d.threadsafe_generator import threadsafe_generator

class IBISGenerator(object):
    'Generates data for Keras'
    def __init__(self, ibis_data, input_shape=(96, 96, 96), batch_size=2, shuffle=True, tax_glob_group="A_eukaryota", num_representatives=2, start_index=0, end_index=None):
        'Initialization'
        self.input_shape = input_shape
        self.batch_size = batch_size
        self.shuffle = shuffle

        # observations_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/unique_obs_bs-8.tab"
        # resfaces_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/unique_molresface_resi_resn_NR_all.tab"

        # observations = pd.read_table(observations_f, header=0, usecols=(0,2))
        # resfaces = pd.read_table(resfaces_f, header=0)
        # data = pd.merge(resfaces, observations, left_on="unique_obs_int", right_on="uniq_obs_int_id")
        ibis_data = pd.read_table(ibis_data, sep="\t")
        self.data = ibis_data.loc[(ibis_data["tax_glob_group"] == tax_glob_group) & (ibis_data["n"] >= num_representatives)]

        if start_index < 1:
            start_index *= self.data.shape[0] 

        if end_index is None:
            end_index = self.data.shape[0]
        elif end_index < 1:
            end_index *= self.data.shape[0]

        if start_index != 0 and end_index != self.data.shape[0]:
            self.data = self.data.iloc[start_index:end_index]

        self.indexes = np.arange(self.data.shape[0])
        self.steps_per_epoch = self.data.shape[0]/self.batch_size

    @classmethod
    def get_training_and_validation(cls, ibis_data, input_shape=(96, 96, 96), batch_size=2, shuffle=True, tax_glob_group="A_eukaryota", num_representatives=2, data_split=0.8):
        train = cls(
            ibis_data,
            input_shape=input_shape, 
            batch_size=batch_size, 
            shuffle=shuffle, 
            tax_glob_group=tax_glob_group, 
            num_representatives=num_representatives, 
            end_index=data_split)
        validate = cls(
            ibis_data,
            input_shape=input_shape, 
            batch_size=batch_size, 
            shuffle=shuffle, 
            tax_glob_group=tax_glob_group, 
            num_representatives=num_representatives, 
            end_index=1-data_split)
        return train, validate

    @threadsafe_generator
    def generate(self):
        'Generates batches of samples'
        # Infinite loop
        while 1:
            print "START"
            # Generate order of exploration of dataset
            if self.shuffle == True:
              np.random.shuffle(self.indexes)

            # Generate batches
            imax = int(len(self.indexes)/self.batch_size)
            for i in xrange(imax):
                # Find list of IDs
                indexes = self.indexes[i*self.batch_size:(i+1)*self.batch_size]
                print "    BATCH", indexes

                # Generate data
                X, y = self.__data_generation(indexes)

                yield X, y
  
    def __data_generation(self, indexes):
        'Generates data of batch_size samples' # X : (n_samples, v_size, v_size, v_size, n_channels)
        # Initialization
        X = None
        y = None
  
        # Generate data
        for i, index in enumerate(indexes):
            datum = self.data.iloc[i]
            pdb_data, pdb_truth = Structure.features_from_string(datum["pdb"], datum["chain"], datum["resi"], input_shape=self.input_shape)

            if X is None:
                X = np.zeros([self.batch_size]+list(pdb_data.shape))
                y = np.zeros([self.batch_size]+list(pdb_truth.shape))

            X[i] = pdb_data
            y[i] = pdb_truth

        return X, y
