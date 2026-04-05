import argparse

import keras
import numpy as np

from common.autoencoder import Autoencoder, load_ae_parameters

# Take ae_id as CLI input
parser = argparse.ArgumentParser(description='Grid Search for autoencoders')
parser.add_argument('ae_id')

# Parse the ae_id input
args = parser.parse_args()
ae_id = args.ae_id

ae_id, data_frame, data_type, architecture, autoencoder_seed = load_ae_parameters(ae_id)

# Set random seed, intialize RNG and generate seeds for weights initlization
keras.utils.set_random_seed(autoencoder_seed)
rng = np.random.default_rng(seed=autoencoder_seed)
layer_seeds = rng.integers(low=0, high=1000, size=len(architecture*2))

autoencoder = Autoencoder(data_frame, data_type, architecture, layer_seeds, ae_id)

autoencoder.cross_validate(plot_learning_curve=True, plot_reconstruction=True, save=True)
