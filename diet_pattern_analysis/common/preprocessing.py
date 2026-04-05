import pandas as pd


def proportional_scale(data: pd.DataFrame):
    data_row_sums = data.sum(axis=1)
    data_scaled = data.div(data_row_sums, axis=0)

    return data_scaled