import pandas as pd

df = pd.read_csv('results/autoencoder/autoencoder_cv_results.csv')
df['val_loss_cv'] = (df['val_loss_std'] / df['val_loss_mean']) * 100

df = df[df['val_loss_mean'] < 0.05]
df = df[df['val_loss_cv'] < 20]

df['AE_ID'].to_csv('results/autoencoder/filtered_AE_IDs.txt', index=False, header=False)
