import os

import pandas as pd

results_dir = 'results/autoencoder'
result_folder_names = [f for f in os.listdir(results_dir) if not f.startswith('.')]

df_autoencoder_cv_results = pd.DataFrame()

for folder_name in result_folder_names:
    folder_path = os.path.join(results_dir, str(folder_name))
    cv_result_path = os.path.join(folder_path, 'cv_result.csv')

    if os.path.exists(cv_result_path):
        df_cv_result = pd.read_csv(cv_result_path)
        df_autoencoder_cv_results = pd.concat(
            [df_autoencoder_cv_results, df_cv_result], ignore_index=True
        )

df_autoencoder_cv_results.to_csv('results/autoencoder/autoencoder_cv_results.csv', index=False)
