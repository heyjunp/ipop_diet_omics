import argparse

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

from common.config import non_omic_cols
from common.utils import merge_dataframes

P_CUTOFF = 0.05
BH_P_CUTOFF = 0.2

OMIC_TYPES = ['clinic', 'cytokines', 'metabolites', 'stool_microbes']
OMIC_LABELS = {
    'clinic': 'Clinical Labs',
    'cytokines': 'Cytokines',
    'metabolites': 'Metabolome',
    'stool_microbes': 'Stool Microbiome',
}

COMPARISONS = [
    dict(group_col='sspg_status', label='IS vs IR'),
    dict(group_col='diet_cluster', label='DP 0 vs DP 1'),
    dict(
        group_col='diet_cluster',
        filter_col='sspg_status',
        filter_value='IS',
        label='IS - DP 0 vs DP 1',
    ),
    dict(
        group_col='diet_cluster',
        filter_col='sspg_status',
        filter_value='IR',
        label='IR - DP 0 vs DP 1',
    ),
]


def load_omic(omic_type):
    df = pd.read_csv(f'data/omics/{omic_type}.csv')
    omic_cols = [c for c in df.columns if c not in non_omic_cols]

    df = df[df['CL4'] == 'Healthy']
    df = df.dropna(subset=['Time'])

    return df, omic_cols


def kw_test(df, omic_cols, group_col, filter_col=None, filter_value=None):
    if filter_col is not None:
        df = df[df[filter_col] == filter_value]

    p_values = []
    for col in omic_cols:
        df_test = df[[col, group_col]].dropna()
        groups = df_test.groupby(group_col)[col].apply(list)
        try:
            _, p = kruskal(*groups)
        except ValueError:
            continue
        if pd.isna(p):
            continue
        p_values.append([col, p])

    df_results = pd.DataFrame(p_values, columns=['Omic', 'p-value'])
    _, bh_p_values, _, _ = multipletests(df_results['p-value'], method='fdr_bh')
    df_results['BH p-value'] = bh_p_values

    significant = int((df_results['BH p-value'] < BH_P_CUTOFF).sum())
    return significant, df_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ae_id')
    args = parser.parse_args()
    ae_id = args.ae_id

    df_diet_clusters = pd.read_csv(f'results/autoencoder/{ae_id}/diet_clusters.csv')

    if 'Time' in df_diet_clusters.columns:
        df_demographics = pd.read_csv('data/diet/df_a/demographics.csv')
        df_sspg = df_demographics[['SubjectID', 'Time', 'sspg_status']]
    else:
        df_demographics = pd.read_csv('data/diet/df_b/demographics.csv')
        df_sspg = df_demographics[['SubjectID', 'sspg_status']]

    rows = []

    for omic_type in OMIC_TYPES:
        df_omic, omic_cols = load_omic(omic_type)
        df_omic = merge_dataframes(df_omic, df_diet_clusters)
        df_omic = merge_dataframes(df_omic, df_sspg)

        for comparison in COMPARISONS:
            significant, df_results = kw_test(
                df_omic,
                omic_cols,
                comparison['group_col'],
                comparison.get('filter_col'),
                comparison.get('filter_value'),
            )

            rows.append([
                OMIC_LABELS[omic_type],
                comparison['label'],
                significant,
            ])

    df_summary = pd.DataFrame(rows, columns=['Omic', 'Comparison Type', 'Significant Markers'])

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(
        df_summary,
        y='Significant Markers',
        x='Omic',
        hue='Comparison Type',
        palette=['yellow', 'green', 'blue', 'red'],
        ax=ax,
    )
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Comparison Group')
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_ylabel('Number of significant biomarker differences')
    fig.tight_layout()
    fig.savefig(f'results/autoencoder/{ae_id}/omic_differences.png', dpi=600)

    df_summary.to_csv(f'results/autoencoder/{ae_id}/omic_differences.csv', index=False)


if __name__ == '__main__':
    main()
