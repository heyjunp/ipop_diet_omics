import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

NUTRIENT_LIST = [
    'MonSac', 'BetaCaro', 'SolFib', 'Fib', 'Biot', 'MonoFat', 'Retinol',
    'W3', 'Folate', 'Fe',
    'VitA_IU', 'VitB1', 'VitB2', 'VitB3', 'VitB6', 'VitB12',
    'VitC', 'VitD_IU', 'VitE_a_Toco', 'VitK',
]

P_CUTOFF = 0.25

COLOR_MAP = {
    'IS': '#298C8C',
    'IR': '#F1A226',
    'Both': 'red',
    'None': '#bec4c0',
}

COLUMN_RENAMES = {
    'data_set1': 'X',
    'data_set2': 'Y',
    'cor': 'Correlation',
    'p': 'p-value',
}


def get_significance(p_is, p_ir):
    sig_is = p_is < P_CUTOFF
    sig_ir = p_ir < P_CUTOFF
    if sig_is and sig_ir:
        return 'Both'
    elif sig_is:
        return 'IS'
    elif sig_ir:
        return 'IR'
    return 'None'


def main():
    df_is = pd.read_csv('data/correlation/nutrients_stool_microbes/cor_data_IS.csv')
    df_ir = pd.read_csv('data/correlation/nutrients_stool_microbes/cor_data_IR.csv')

    df_is = df_is.rename(columns=COLUMN_RENAMES)
    df_ir = df_ir.rename(columns=COLUMN_RENAMES)

    df = pd.merge(df_is, df_ir, on=['X', 'Y'], suffixes=('_IS', '_IR'))

    genus_cols = [c for c in df['Y'].unique() if 'genus' in c]
    df = df[df['X'].isin(NUTRIENT_LIST) & df['Y'].isin(genus_cols)]

    df['Significance'] = df.apply(
        lambda row: get_significance(row['p_adjust_IS'], row['p_adjust_IR']), axis=1
    )

    df_plot = df.pivot(index='X', columns='Y', values='Significance')

    value_map = {'None': 0, 'IS': 1, 'IR': 2, 'Both': 3}
    df_numeric = df_plot.map(lambda x: value_map[x])

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(
        df_numeric,
        linewidths=0.005,
        cmap=list(COLOR_MAP.values()),
        cbar=False,
        ax=ax,
    )
    ax.set_ylabel('Vitamins')
    ax.set_xlabel('Microbes')

    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='black')
        for c in COLOR_MAP.values()
    ]
    ax.legend(legend_elements, COLOR_MAP.keys(), loc='center left', bbox_to_anchor=(1, 0.5))

    fig.tight_layout()
    fig.savefig(
        'results/nutrient_microbe_significance.png', bbox_inches='tight', dpi=600
    )

    df_plot.reset_index().to_csv('results/nutrient_microbe_significance.csv', index=False)


if __name__ == '__main__':
    main()
