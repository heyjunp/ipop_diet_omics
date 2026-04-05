import os

import pandas as pd

results_dir = 'results/autoencoder'
result_folder_names = [f for f in os.listdir(results_dir) if not f.startswith('.')]

df_clustering_results = pd.DataFrame()
df_omic_summary_results = pd.DataFrame()

for folder_name in result_folder_names:
    folder_path = os.path.join(results_dir, str(folder_name))
    clustering_result_path = os.path.join(folder_path, 'clustering_metrics.csv')

    if os.path.exists(clustering_result_path):
        df_clustering_result = pd.read_csv(clustering_result_path)
        df_clustering_results = pd.concat(
            [df_clustering_results, df_clustering_result], ignore_index=True
        )

df_clustering_results.to_csv('results/autoencoder/clustering_results.csv', index=False)
