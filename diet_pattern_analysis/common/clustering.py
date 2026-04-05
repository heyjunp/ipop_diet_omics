import os

import keras
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score

from common import config
from common.autoencoder import load_ae_parameters
from common.preprocessing import proportional_scale
from common.utils import merge_dataframes


class Clustering:

    def __init__(self, ae_id: str):
        self.ae_id = ae_id
        self.metrics = {}
        self.__latent_data = None

        _, data_frame, data_type, _, _ = load_ae_parameters(ae_id)

        self.data_frame = data_frame
        self.data_type = data_type

        self.__load_food_data()
        self.__load_encoder()

    def __load_encoder(self):
        encoder_path = os.path.join('results/autoencoder/', self.ae_id, 'encoder.keras')
        self.encoder: keras.Model = keras.models.load_model(encoder_path)

    def __load_food_data(self):
        food_path = os.path.join('data/diet', self.data_frame, f'{self.data_type}.csv')
        self.df_food = pd.read_csv(food_path)

    def get_latent_data(self):
        if self.__latent_data is None:
            food_data = self.df_food[config.food_group_cols]
            food_data = proportional_scale(food_data)

            latent_data = self.encoder.predict(food_data)
            self.__latent_data = latent_data

        return self.__latent_data

    def fit(self, n_clusters: int = 2, save: bool = False):
        latent_data = self.get_latent_data()

        k_means = KMeans(n_clusters=n_clusters, n_init=10, random_state=config.random_state)
        k_means.fit(latent_data)

        subject_ids = self.df_food['SubjectID'].to_list()
        diet_clusters = k_means.labels_

        df_diet_clusters = pd.DataFrame(
            {
                'SubjectID': subject_ids,
                'diet_cluster': diet_clusters
            }
        )

        if 'Time' in self.df_food.columns:
            df_diet_clusters['Time'] = self.df_food['Time'].to_list()
            df_diet_clusters = df_diet_clusters.reindex(
                columns=['SubjectID', 'Time', 'diet_cluster']
            )

        self.df_diet_clusters = df_diet_clusters

        if save:
            df_diet_clusters.to_csv(
                os.path.join('results/autoencoder', self.ae_id, 'diet_clusters.csv'),
                index=False
            )

        return df_diet_clusters

    def pca_plot(self, ax: plt.Axes):
        latent_data = self.get_latent_data()

        pca = PCA(n_components=2, random_state=config.random_state)
        pca_data = pca.fit_transform(latent_data)

        sns.scatterplot(
            x=pca_data[:, 0],
            y=pca_data[:, 1],
            hue=self.df_diet_clusters['diet_cluster'].to_list(),
            palette="Set2",
            ax=ax
        )

        ax.set_title('PCA Plot of Autoencoder Features')
        ax.set_xlabel('PCA Dimension 1')
        ax.set_ylabel('PCA Dimension 2')

        return

    def tsne_plot(self, ax: plt.Axes):
        latent_data = self.get_latent_data()

        tsne = TSNE(
            n_components=2,
            perplexity=60,
            learning_rate=100,
            n_iter=1000,
            random_state=config.random_state
        )
        tsne_data = tsne.fit_transform(latent_data)

        sns.scatterplot(
            x=tsne_data[:, 0],
            y=tsne_data[:, 1],
            hue=self.df_diet_clusters['diet_cluster'].to_list(),
            palette="Set2", ax=ax
        )

        ax.set_title('t-SNE Plot of Autoencoder Features')
        ax.set_xlabel('t-SNE Dimension 1')
        ax.set_ylabel('t-SNE Dimension 2')

        return

    def compute_metrics(self, save: bool = False):
        silhouette = silhouette_score(
            self.get_latent_data(), self.df_diet_clusters['diet_cluster'].to_list()
        )

        cluster_size = self.df_diet_clusters['diet_cluster'].value_counts().to_list()

        cluster_value_counts = self.df_diet_clusters['diet_cluster'].value_counts()
        cluster_size_percentage = (
            cluster_value_counts / len(self.df_diet_clusters['diet_cluster'])
        ).to_list()
        cluster_size_percentage_std = np.std(cluster_size_percentage)

        (
            cluster_wise_ir_counts,
            cluster_wise_ir_proportion,
            cluster_wise_ir_proportion_std
        ) = self.get_cluster_wise_proportion('sspg_status', 'IR')

        (
            cluster_wise_pre_T2D_counts,
            cluster_wise_pre_T2D_proportion,
            cluster_wise_pre_T2D_proportion_std
        ) = self.get_cluster_wise_proportion('baseline_a1c_status', 'preDM/T2D')

        self.metrics['silhouette'] = silhouette

        self.metrics['cluster_size'] = cluster_size
        self.metrics['cluster_size_percentage'] = cluster_size_percentage
        self.metrics['cluster_size_percentage_std'] = cluster_size_percentage_std

        self.metrics['cluster_wise_ir_counts'] = cluster_wise_ir_counts
        self.metrics['cluster_wise_ir_proportion'] = cluster_wise_ir_proportion
        self.metrics['cluster_wise_ir_proportion_std'] = cluster_wise_ir_proportion_std

        self.metrics['cluster_wise_pre_T2D_counts'] = cluster_wise_pre_T2D_counts
        self.metrics['cluster_wise_pre_T2D_proportion'] = cluster_wise_pre_T2D_proportion
        self.metrics['cluster_wise_pre_T2D_proportion_std'] = cluster_wise_pre_T2D_proportion_std

        if save:
            metrics_df = {key: [value] for key, value in self.metrics.items()}
            df_metrics = pd.DataFrame(metrics_df)
            df_metrics['AE_ID'] = [self.ae_id]

            df_metrics.to_csv(
                os.path.join('results/autoencoder', self.ae_id, 'clustering_metrics.csv'),
                index=False
            )

        return self.metrics

    def get_cluster_wise_proportion(self, demographic_col: str, value: str):
        df_demographics = pd.read_csv(
            os.path.join('data/diet', self.data_frame, 'demographics.csv')
        )
        df_diet_demographics = merge_dataframes(self.df_diet_clusters, df_demographics)

        df_value = df_diet_demographics[df_diet_demographics[demographic_col] == value]
        cluster_wise_value_counts = df_value.groupby('diet_cluster').count()['SubjectID']

        cluster_sizes = df_diet_demographics.groupby('diet_cluster').count()['SubjectID']
        cluster_wise_value_proportion = cluster_wise_value_counts / cluster_sizes

        cluster_wise_value_proportion_std = np.std(cluster_wise_value_proportion)

        return (
            cluster_wise_value_counts.to_list(),
            cluster_wise_value_proportion.to_list(),
            cluster_wise_value_proportion_std
        )
