import pandas as pd


def merge_dataframes(left: pd.DataFrame, right: pd.DataFrame):
    if ("Time" in left.columns) and ("Time" in right.columns):
        df_merged = pd.merge(left, right, on=["SubjectID", "Time"])
    else:
        df_merged = pd.merge(left, right, on=["SubjectID"])

    return df_merged
