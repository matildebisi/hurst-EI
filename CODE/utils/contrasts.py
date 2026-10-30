"""Helpers for task-pair statistical contrast maps."""

import re
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def parse_pairwise_key(key, task_names):
    """
    Find the two task names in the key string and return them 
    as a tuple (task1, task2).
    """
    matches = []
    for task in task_names:
        for m in re.finditer(re.escape(task), key):
            matches.append((m.start(), task))

    matches = sorted(matches, key=lambda x: x[0])
    ordered = []
    for _, task in matches:
        if task not in ordered:
            ordered.append(task)

    if len(ordered) < 2:
        raise ValueError(f"Unable to identify two tasks in key='{key}'")

    return ordered[0], ordered[1]


def make_task_diverging_cmap(task1, task2, task_colors, n=256):
    """
    Create a personalized diverging colormap for two tasks, where 
    one task is associated with negative values and the other with positive values.
    """
    color_neg = task_colors[task2]   
    color_pos = task_colors[task1]   

    cmap = LinearSegmentedColormap.from_list(
        f"{task1}_vs_{task2}_cmap",
        [color_neg, "#FFFFFF", color_pos],
        N=n
    )
    return cmap


def get_pairwise_keys_from_df(df):
    """
    Find columns in the DataFrame that represent pairwise t-test results
    """
    keys = []
    for col in df.columns:
        if col.startswith("t_"):
            key = col[2:]
            if f"sig_{key}_fdr" in df.columns:
                keys.append(key)
    return keys


def build_plot_matrix(pair_df, plot_rows, plot_cols, session_order, only_significant=True):

    mat = np.full((len(plot_rows), len(plot_cols)), np.nan, dtype=float)
    annot = np.full((len(plot_rows), len(plot_cols)), "", dtype=object)

    for i, r in enumerate(plot_rows):
        for j, c in enumerate(plot_cols):
            if session_order.index(r) >= session_order.index(c):
                continue

            row = pair_df[(pair_df["s1"] == r) & (pair_df["s2"] == c)]
            if row.empty:
                continue

            row = row.iloc[0]

            if only_significant and not bool(row["sig"]):
                continue

            mat[i, j] = float(row["effect"])
            annot[i, j] = f"{row['effect']:+.2f}"

    return mat, annot
