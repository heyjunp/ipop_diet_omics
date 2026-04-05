import argparse

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from common import config
from common.clustering import Clustering


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ae_id')
    args = parser.parse_args()
    ae_id = args.ae_id

    clustering = Clustering(ae_id)
    clustering.fit()

    latent_data = clustering.get_latent_data()
    diet_clusters = clustering.df_diet_clusters['diet_cluster'].values

    tsne = TSNE(n_components=2, perplexity=60, learning_rate=100, n_iter=1000, random_state=config.random_state)
    tsne_data = tsne.fit_transform(latent_data)

    pca = PCA(n_components=2, random_state=config.random_state)
    pca_data = pca.fit_transform(latent_data)

    fig, axes = plt.subplots(2, 1, figsize=(6, 9))

    sns.scatterplot(x=tsne_data[:, 0], y=tsne_data[:, 1], hue=diet_clusters, palette='Set2', ax=axes[0])
    axes[0].set_title('t-SNE Plot of Autoencoder Features')
    axes[0].set_xlabel('t-SNE Dimension 1')
    axes[0].set_ylabel('t-SNE Dimension 2')
    axes[0].get_legend().remove()

    sns.scatterplot(x=pca_data[:, 0], y=pca_data[:, 1], hue=diet_clusters, palette='Set2', ax=axes[1])
    axes[1].set_title('PCA Plot of Autoencoder Features')
    axes[1].set_xlabel('PCA Dimension 1')
    axes[1].set_ylabel('PCA Dimension 2')

    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, ['DP 0', 'DP 1'], loc='upper right', bbox_to_anchor=(1.15, 1))
    axes[1].get_legend().remove()

    fig.tight_layout()
    fig.savefig(f'results/autoencoder/{ae_id}/tsne_pca_plot.png', bbox_inches='tight', dpi=900)


if __name__ == '__main__':
    main()
