import os

import pandas as pd
import numpy as np
import traceback

try:
    import tensorflow as tf
except ImportError:
    tf = None

from unet3d.threadsafe_generator import threadsafe_generator

# try:
from molmimic.biopdbtools import Structure, InvalidPDB
# except ImportError:
#     Structure = None
#     InvalidPDB = None
#     print "Warning: This module can be used to read in the data file, but will be unable to load any fo the features"



class IBISGenerator(object):
    'Generates data for Keras'
    def __init__(self, ibis_data, input_shape=(96, 96, 96), batch_size=2, shuffle=True, tax_glob_group="A_eukaryota", num_representatives=2, start_index=0, end_index=None, only_aa=False):
        'Initialization'
        self.input_shape = input_shape
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.only_aa = only_aa

        # observations_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/unique_obs_bs-8.tab"
        # resfaces_f = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tandemPPI/databases/unique_molresface_resi_resn_NR_all.tab"

        # observations = pd.read_table(observations_f, header=0, usecols=(0,2))
        # resfaces = pd.read_table(resfaces_f, header=0)
        # data = pd.merge(resfaces, observations, left_on="unique_obs_int", right_on="uniq_obs_int_id")
        data = pd.read_table(ibis_data, sep="\t")
        self.data = data.loc[(data["tax_glob_group"] == tax_glob_group) & (data["n"] >= num_representatives)]

        try:
            skip = "{}_skip.tab".format(os.path.splitext(ibis_data)[0])
            print skip
            with open(skip) as skip_f:
                skip_ids = [int(l.rstrip()) for l in skip_f if l]
            osize = self.data.shape[0]
            self.data = self.data.loc[~self.data["unique_obs_int"].isin(skip_ids)]
            print "Skipping {} of {}, {} remaining of {}".format(self.data.shape[0], osize, self.data.shape[0], osize)
        except IOError:
            print "No Skip ID file"
            pass

        if start_index < 1:
            start_index *= self.data.shape[0]

        if end_index is None:
            end_index = self.data.shape[0]
        elif end_index < 1:
            end_index *= self.data.shape[0]

        start_index = int(start_index)
        end_index = int(end_index)

        print "starting at", start_index
        print "ending at", end_index

        if start_index != 0 or end_index != self.data.shape[0]:
            print "Resizing"
            self.data = self.data.iloc[start_index:end_index]

        self.indexes = np.arange(self.data.shape[0])
        self.steps_per_epoch = self.data.shape[0]/self.batch_size

    @classmethod
    def get_training_and_validation(cls, ibis_data, input_shape=(96, 96, 96), batch_size=2, shuffle=True, tax_glob_group="A_eukaryota", num_representatives=2, data_split=0.8, only_aa=False):
        train = cls(
            ibis_data,
            input_shape=input_shape,
            batch_size=batch_size,
            shuffle=shuffle,
            tax_glob_group=tax_glob_group,
            num_representatives=num_representatives,
            end_index=data_split,
            only_aa=only_aa)
        validate = cls(
            ibis_data,
            input_shape=input_shape,
            batch_size=batch_size,
            shuffle=shuffle,
            tax_glob_group=tax_glob_group,
            num_representatives=num_representatives,
            start_index=data_split,
            only_aa=only_aa)
        return train, validate

    #@threadsafe_generator
    def generate(self):
        'Generates batches of samples'
        # Infinite loop
        epoch = 0
        while 1:
            epoch += 1
            print "START"
            # Generate order of exploration of dataset
            if self.shuffle:
              np.random.shuffle(self.indexes)

            # Generate batches
            imax = int(len(self.indexes)/self.batch_size)
            for i in xrange(imax):
                # Find list of IDs
                indexes = self.indexes[i*self.batch_size:(i+1)*self.batch_size]
                print "    BATCH", indexes

                # Generate data
                X, y = self.__data_generation(indexes, epoch=epoch, batch=i)

                if X.shape[0] > 0:
                    yield X, y

    def __data_generation(self, indexes, epoch=None, batch=None):
        'Generates data of batch_size samples' # X : (n_samples, v_size, v_size, v_size, n_channels)
        # Initialization
        X = np.zeros((self.batch_size, self.input_shape[0], self.input_shape[1], self.input_shape[1], 21 if self.only_aa else Structure.nFeatures)) #or Structure.Features
        y = np.zeros((self.batch_size, self.input_shape[0], self.input_shape[1], self.input_shape[1], 1))
        # Generate data
        failures = 0
        for batch_index, index in enumerate(indexes):
            datum = self.data.iloc[index]
            print "Running {} ({}.{}): {}".format(datum["unique_obs_int"], datum["pdb"], datum["chain"], ",".join(["{}{}".format(i,n) for i, n in zip(datum["resi"].split(","), datum["resn"].split(","))]))
            if self.only_aa:
                print "AA only"

            try:
                (data_idx, data), (truth_idx, truth) = Structure.features_from_string(
                    datum["pdb"], 
                    datum["chain"], 
                    datum["resi"], 
                    id=datum["unique_obs_int"], 
                    input_shape=self.input_shape, 
                    batch_index=batch_index-failures, 
                    only_aa=self.only_aa)
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                failures += 1

                X = np.resize(X, (X.shape[0]-1,)+X.shape[1:])
                y = np.resize(y, (y.shape[0]-1,)+y.shape[1:])

                trace = traceback.format_exc()
                print "Error:", trace
                with open("{}_{}_{}_{}_{}.error".format(datum["pdb"], datum["chain"], datum["unique_obs_int"], epoch, batch), "w") as ef:
                    print >> ef, trace
            	continue

            X[data_idx.T.tolist()] = data.flatten()
            y[truth_idx.T.tolist()] = truth

        return X, y
