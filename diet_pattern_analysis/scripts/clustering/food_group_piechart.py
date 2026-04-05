import argparse

import pandas as pd
import matplotlib.pyplot as plt

from common.config import food_group_cols
from common.autoencoder import load_ae_parameters
from common.preprocessing import proportional_scale
from common.utils import merge_dataframes

COLUMN_LABELS = {
    'Bread_Products': 'Bread Products',
    'Dairy_Products': 'Dairy Products',
    'Fats and Oils': 'Fats & Oils',
    'Meat_Poultry': 'Meat & Poultry',
    'Non_Starchy_Vegetables': 'Non Starchy Vegetables',
    'Nuts_Seeds': 'Nuts & Seeds',
    'Starchy_Vegetables': 'Starchy Vegetables',
    'Sweets_Sweeteners': 'Sweets & Sweeteners',
    'Cereals_Grains': 'Cereals & Grains',
    'Sauces, Spices, and Herbs': 'Sauces & Spices & Herbs',
    'Legume_Products': 'Legume Products',
    'Mixed_Meal': 'Pizza & Sandwich & Wrap',
    'Rice_Pasta_Noodles': 'Rice & Pasta & Noodles',
}

COLOR_MAP = {
    'Bread Products': '#1f77b4',
    'Coffee': '#aec7e8',
    'Dairy Products': '#ff7f0e',
    'Fats & Oils': '#ffbb78',
    'Fruits': '#2ca02c',
    'Non Starchy Vegetables': '#98df8a',
    'Meat & Poultry': '#d62728',
    'Sweets & Sweeteners': '#ff9896',
    'Sauces & Spices & Herbs': '#9467bd',
    'Snacks': '#c5b0d5',
    'Nuts & Seeds': '#8c564b',
    'Legume Products': '#c49c94',
    'Seafood': '#e377c2',
    'Salad': '#f7b6d2',
    'Starchy Vegetables': '#7f7f7f',
    'Rice & Pasta & Noodles': '#c7c7c7',
    'Pizza & Sandwich & Wrap': '#bcbd22',
    'Eggs': '#dbdb8d',
    'Tea': '#17becf',
    'Cereals & Grains': '#9edae5',
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ae_id')
    args = parser.parse_args()
    ae_id = args.ae_id

    _, data_frame, data_type, _, _ = load_ae_parameters(ae_id)

    df_food = pd.read_csv(f'data/diet/{data_frame}/{data_type}.csv')
    df_clusters = pd.read_csv(f'results/autoencoder/{ae_id}/diet_clusters.csv')
    df = merge_dataframes(df_food, df_clusters)

    df_medians = df[food_group_cols + ['diet_cluster']].groupby('diet_cluster').median()
    df_medians = proportional_scale(df_medians)
    df_medians = df_medians.rename(columns=COLUMN_LABELS)

    fig, axes = plt.subplots(nrows=1, ncols=len(df_medians), figsize=(14, 8))

    for i, (_, row) in enumerate(df_medians.iterrows()):
        row = row.replace(0, None).dropna()
        colors = [COLOR_MAP[col] for col in row.index]
        axes[i].pie(row.to_list(), autopct='%1.1f%%', colors=colors, radius=1, startangle=90)
        axes[i].set_title(f'Diet Pattern {i}')

    legend_patches = [
        plt.Rectangle((0, 0), 1, 1, facecolor=color)
        for color in COLOR_MAP.values()
    ]
    fig.legend(legend_patches, COLOR_MAP.keys(), loc='lower center', ncol=5)
    fig.tight_layout()
    fig.savefig(
        f'results/autoencoder/{ae_id}/food_groups_piechart.png', bbox_inches='tight', dpi=900
    )


if __name__ == '__main__':
    main()
