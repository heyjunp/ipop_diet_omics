import ast
import os

import keras

import tensorflow as tf
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from typing import List
from keras.layers import Input, Dense
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold

from common import config
from common.preprocessing import proportional_scale


def load_ae_parameters(ae_id: str):
    df_ae_params = pd.read_csv('data/inputs/ae_params.csv')
    params = df_ae_params[df_ae_params['AE_ID'] == ae_id].to_numpy()[0]

    ae_id, data_frame, data_type, architecture, autoencoder_seed = params
    architecture = ast.literal_eval(architecture)

    return ae_id, data_frame, data_type, architecture, autoencoder_seed


class Autoencoder:

    def __init__(
        self,
        data_frame: str,
        data_type: str,
        architecture: List[int],
        layer_seeds: List[int],
        ae_id: str = 'AEcustom',
    ):
        self.ae_id = ae_id
        self.architecture = architecture
        self.layer_seeds = layer_seeds

        self.df_food = pd.read_csv(f'data/diet/{data_frame}/{data_type}.csv')
        self.input_shape = len(config.food_group_cols)

        df_demographics = pd.read_csv(f'data/diet/{data_frame}/demographics.csv')
        df_sspg_status = df_demographics.drop_duplicates(subset=['SubjectID'])
        df_sspg_status = df_sspg_status[['SubjectID', 'sspg_status']]
        self.df_sspg_status = df_sspg_status

        self.batch_size = 16 if data_frame == "df_a" else 8

        return

    def __create_models(self):
        layer_seed_idx = 0

        encoder_layers = []
        encoder_layers.append(Input(shape=self.input_shape))

        # Add the encoder layers
        for layer in self.architecture:
            encoder_layers.append(
                Dense(
                    layer,
                    kernel_initializer=keras.initializers.glorot_normal(
                        seed=self.layer_seeds[layer_seed_idx]
                    ),
                    activation='relu'
                )
            )

            layer_seed_idx += 1

        decoder_layers = []
        # Add the decoder layers
        for layer in reversed(self.architecture[:-1]):
            decoder_layers.append(
                Dense(
                    layer,
                    kernel_initializer=keras.initializers.glorot_normal(
                        seed=self.layer_seeds[layer_seed_idx]
                    ),
                    activation='relu'
                )
            )

            layer_seed_idx += 1

        decoder_layers.append(
            Dense(
                    self.input_shape,
                    kernel_initializer=keras.initializers.glorot_normal(
                        seed=self.layer_seeds[layer_seed_idx]
                    ),
                    activation='relu'
                )
        )

        autoencoder = keras.Sequential(encoder_layers + decoder_layers)
        autoencoder.compile(optimizer='adam', loss='mean_squared_error')

        encoder = keras.Sequential(encoder_layers)

        return autoencoder, encoder

    def __fit(self, train_data, val_data=None):
        autoencoder, encoder = self.__create_models()

        _val_data = (val_data, val_data) if val_data is not None else None

        history = autoencoder.fit(
            train_data,
            train_data,
            validation_data=_val_data,
            epochs=400,
            batch_size=self.batch_size,
            verbose=0
        )

        return autoencoder, encoder, history

    def __plot_learning_curve(self, histories: List, ax):
        train_losses_history = []
        val_losses_history = []

        for history in histories:
            train_losses_history.append(history.history['loss'])
            val_losses_history.append(history.history['val_loss'])

        train_losses_history = np.array(train_losses_history)
        val_losses_history = np.array(val_losses_history)

        mean_train_loss = np.mean(train_losses_history, axis=0)
        std_train_loss = np.std(train_losses_history, axis=0)
        mean_val_loss = np.mean(val_losses_history, axis=0)
        std_val_loss = np.std(val_losses_history, axis=0)

        epochs = range(1, len(mean_train_loss) + 1)

        ax.plot(epochs, mean_train_loss, label='Mean Training Loss', color='blue')
        ax.fill_between(
            epochs,
            mean_train_loss - std_train_loss,
            mean_train_loss + std_train_loss,
            alpha=0.1,
            color='blue'
        )

        ax.plot(epochs, mean_val_loss, label='Mean Validation Loss', color='red')
        ax.fill_between(
            epochs,
            mean_val_loss - std_val_loss,
            mean_val_loss + std_val_loss,
            alpha=0.1,
            color='red'
        )

        ax.set_title('Training and Validation Loss')
        ax.set_xlabel('Epochs')
        ax.set_ylabel('Loss')
        ax.legend()
        ax.grid(True)

    def __plot_reconstruction(
        self,
        autoencoder: keras.Model,
        train,
        val,
        ax,
    ):
        pca = PCA(n_components=2, random_state=config.random_state)

        x_values = []
        y_values = []
        original_reconstructed = []
        train_val = []

        for data, data_split in zip([train, val], ['train', 'val']):
            if data is None:
                continue

            data_pca = pca.fit_transform(data)
            data_reconstructed = autoencoder.predict(data)
            data_reconstructed_pca = pca.fit_transform(data_reconstructed)

            x_values.extend(data_pca[:, 0])
            x_values.extend(data_reconstructed_pca[:, 0])

            y_values.extend(data_pca[:, 1])
            y_values.extend(data_reconstructed_pca[:, 1])

            original_reconstructed.extend(['original'] * len(data))
            original_reconstructed.extend(['reconstructed'] * len(data_reconstructed_pca))

            train_val.extend([data_split] * len(data))
            train_val.extend([data_split] * len(data_reconstructed_pca))

        df_plot = pd.DataFrame({
            'x': x_values,
            'y': y_values,
            'original_reconstructed': original_reconstructed,
            'train_val': train_val,
        })

        sns.scatterplot(
            df_plot, x='x', y='y', hue='original_reconstructed', style='train_val', ax=ax
        )

        ax.set_xlabel('PCA Dimension 1')
        ax.set_ylabel('PCA Dimension 2')
        ax.get_legend().remove()

    def __get_results_dir(self):
        results_dir = f'results/autoencoder/{self.ae_id}'
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)

        return results_dir

    def cross_validate(
        self,
        folds: int = 5,
        plot_learning_curve: bool = False,
        plot_reconstruction: bool = False,
        save: bool = False,
    ):
        k_fold = StratifiedKFold(n_splits=folds, shuffle=True, random_state=config.random_state)

        subject_ids = self.df_sspg_status['SubjectID'].to_list()
        sspg_status = self.df_sspg_status['sspg_status'].to_list()

        losses = []
        val_losses = []
        histories = []

        fig_reconstruction, axes = plt.subplots(1, folds, figsize=(5*folds, 5))

        for idx, (train_idx, val_idx) in enumerate(k_fold.split(subject_ids, sspg_status)):
            train_subject_ids = [subject_ids[i] for i in train_idx]
            val_subject_ids = [subject_ids[i] for i in val_idx]

            train = self.df_food[self.df_food['SubjectID'].isin(train_subject_ids)]
            val = self.df_food[self.df_food['SubjectID'].isin(val_subject_ids)]

            train = train[config.food_group_cols]
            val = val[config.food_group_cols]

            train = proportional_scale(train)
            val = proportional_scale(val)

            autoencoder, _, history = self.__fit(train, val)

            if plot_reconstruction:
                ax = axes.flatten()[idx]
                self.__plot_reconstruction(autoencoder, train, val, ax)

            loss = autoencoder.evaluate(train, train)
            val_loss = autoencoder.evaluate(val, val)

            losses.append(loss)
            val_losses.append(val_loss)
            histories.append(history)

            tf.keras.backend.clear_session()

        fig_reconstruction.tight_layout()

        if plot_learning_curve:
            fig_learning_curve, ax = plt.subplots(figsize=(10, 6))
            self.__plot_learning_curve(histories, ax)

            fig_learning_curve.tight_layout()

        if save:
            results_dir = self.__get_results_dir()

            cv_result_cols = [
                'AE_ID',
                'loss_mean',
                'loss_std',
                'val_loss_mean',
                'val_loss_std'
            ]
            cv_result = [
                self.ae_id,
                np.mean(losses),
                np.std(losses),
                np.mean(val_losses),
                np.std(val_losses)
            ]

            df_result = pd.DataFrame([cv_result], columns=cv_result_cols)
            df_result.to_csv(os.path.join(results_dir, 'cv_result.csv'), index=False)

            fig_reconstruction.savefig(
                os.path.join(results_dir, 'cv_reconstruction_plot.png'),
                dpi=600
            )

            fig_learning_curve.savefig(
                os.path.join(results_dir, 'cv_learning_curve.png'),
                dpi=600
            )

        return losses, val_losses, histories

    def train_final_model(
        self,
        plot_reconstruction: bool = False,
        save: bool = False
    ):
        train = self.df_food
        train = train[config.food_group_cols]
        train = proportional_scale(train)

        autoencoder, encoder, history = self.__fit(train)

        if plot_reconstruction:
            fig, ax = plt.subplots(figsize=(5, 5))
            self.__plot_reconstruction(autoencoder, train, None, ax)
            fig.tight_layout()

        if save:
            results_dir = self.__get_results_dir()

            final_result_cols = ['AE_ID', 'loss']
            final_result = [self.ae_id, autoencoder.evaluate(train, train)]

            df_result = pd.DataFrame([final_result], columns=final_result_cols)
            df_result.to_csv(os.path.join(results_dir, 'final_result.csv'), index=False)

            fig.savefig(
                os.path.join(results_dir, 'final_reconstruction_plot.png'),
                dpi=600
            )

            autoencoder.save(os.path.join(results_dir, 'autoencoder.keras'))
            encoder.save(os.path.join(results_dir, 'encoder.keras'))

        return autoencoder, encoder, history
