import argparse
import os

import matplotlib.pyplot as plt

from common.clustering import Clustering

# Take ae_id as CLI input
parser = argparse.ArgumentParser(description='Grid Search for autoencoders')
parser.add_argument('ae_id')

# Parse the ae_id input
args = parser.parse_args()
ae_id = args.ae_id

clustering = Clustering(ae_id)
clustering.fit(save=True)
clustering.compute_metrics(save=True)

fig, axes = plt.subplots(1, 2, figsize=(8, 4))
clustering.pca_plot(axes[0])
clustering.tsne_plot(axes[1])

fig.tight_layout()
fig.savefig(os.path.join(f'results/autoencoder/{ae_id}', 'pca_tsne_plots.png'), dpi=600)
