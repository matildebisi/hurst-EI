"""Parcelwise summaries and nonparametric statistical comparisons."""

import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from scipy.stats import friedmanchisquare
from itertools import combinations
from statsmodels.stats.multitest import multipletests


def mean_parcelwise(data_df, field, task=None):
    sub_data_df = data_df.copy()
    if task is not None:
        sub_data_df = sub_data_df[sub_data_df['Task'] == task]
    mean_hurst = sub_data_df.groupby('Parcel')[field].mean().values
    return mean_hurst

def std_parcelwise(data_df, task, field):
    sub_data_df = data_df.copy()
    sub_data_df = sub_data_df[sub_data_df['Task'] == task]
    std_hurst = sub_data_df.groupby('Parcel')[field].std().values
    return std_hurst

def ci95_parcelwise(data_df, field):
    sub_data_df = data_df.copy()
    ci95_low = sub_data_df.groupby('Parcel')[field].quantile(0.025).values
    ci95_high = sub_data_df.groupby('Parcel')[field].quantile(0.975).values
    return (ci95_high - ci95_low) / 2

def zeta_score(mean_hurst):
    mean = np.mean(mean_hurst)
    std_dev = np.std(mean_hurst)
    zeta_scores = (mean_hurst - mean) / std_dev
    return zeta_scores


def preprocess_df(df, task_col, tasks):
    df = df.copy()
    df = df[df[task_col].isin(tasks)].copy()
    df[task_col] = pd.Categorical(df[task_col], categories=tasks, ordered=True)
    return df


def safe_wilcoxon(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = ~(np.isnan(x) | np.isnan(y))
    x = x[mask]
    y = y[mask]

    n = len(x)
    if n < 2:
        return np.nan, np.nan, n

    diff = x - y
    if np.allclose(diff, 0):
        return 0.0, 1.0, n

    try:
        res = wilcoxon(x, y, zero_method="wilcox", alternative="two-sided", mode="auto")
        return float(res.statistic), float(res.pvalue), n
    
    except ValueError:
        return np.nan, np.nan, n



def friedman_test(df_sub, 
                  selection_values, 
                  other_values, 
                  index_col,
                  columns_col,
                  value_col,
                  selection_col,
                  p_val_sig=0.05, 
                  correction="fdr_bh"):   
     
    overall_rows = []

    for val in selection_values:
        tmp = df_sub[df_sub[selection_col] == val].copy()

        wide = tmp.pivot_table(
            observed=False,
            index=index_col,
            columns=columns_col,
            values=value_col,
            aggfunc="mean"
        ).reindex(columns=other_values)

        complete = wide.dropna()
        n_subj = int(complete.shape[0])

        stat = np.nan
        pval = np.nan

        if n_subj >= 2:
            try:
                arrays = [complete[s].to_numpy() for s in other_values]
                stat, pval = friedmanchisquare(*arrays)
                stat = float(stat)
                pval = float(pval)
            except ValueError:
                pass

        overall_rows.append({
            selection_col: val,
            "n_subjects": n_subj,
            "friedman_chi2": stat,
            "p": pval
        })

    overall_df = pd.DataFrame(overall_rows)
    overall_df["q"] = np.nan
    overall_df["sig"] = False

    valid = overall_df["p"].notna()
    if valid.any():
        reject, qvals, _, _ = multipletests(
            overall_df.loc[valid, "p"].to_numpy(),
            alpha=p_val_sig,
            method=correction
        )
        overall_df.loc[valid, "q"] = qvals
        overall_df.loc[valid, "sig"] = reject

    sig_values = overall_df.loc[overall_df["sig"], selection_col].tolist()

    return overall_df, sig_values


def pairwise_wilcoxon(df_sub, 
                    selection_values, 
                    other_values, 
                    index_col,
                    columns_col,
                    value_col,
                    selection_col,
                    p_val_sig=0.05, 
                    correction="fdr_bh"):  

    pairwise_rows = [] 
    
    pairwise_pairs = list(combinations(other_values, 2))
    for net in selection_values:
        tmp = df_sub[df_sub[selection_col] == net].copy()

        wide = tmp.pivot_table(
            observed=False,
            index=index_col,
            columns=columns_col,
            values=value_col,
            aggfunc="mean"
        ).reindex(columns=other_values)

        for s1, s2 in pairwise_pairs:
            pair = wide[[s1, s2]].dropna()

            if pair.shape[0] >= 2:
                stat, pval, nobs = safe_wilcoxon(pair[s1].values, pair[s2].values)
                effect = np.median(pair[s1].values - pair[s2].values)
            else:
                stat, pval, nobs, effect = np.nan, np.nan, 0, np.nan

            pairwise_rows.append({
                selection_col: net,
                "Comparison": f"{s1} vs {s2}",
                "s1": s1,
                "s2": s2,
                "n_subjects": int(nobs),
                "wilcoxon_W": stat,
                "effect": effect,
                "p": pval
            })

    pairwise_df = pd.DataFrame(pairwise_rows)
    pairwise_df["q"] = np.nan
    pairwise_df["sig"] = False

    valid = pairwise_df["p"].notna()
    if valid.any():
        reject, qvals, _, _ = multipletests(
            pairwise_df.loc[valid, "p"].to_numpy(),
            alpha=p_val_sig,
            method=correction
        )
        pairwise_df.loc[valid, "q"] = qvals
        pairwise_df.loc[valid, "sig"] = reject

    return pairwise_df




def build_summary_table(df, network_col, task_col, value_col, tasks, networks=None):
    summary = (
        df.groupby([network_col, task_col], observed=False)[value_col]
        .agg(mean="mean", std="std", n="count")
        .reset_index()
    )

    if networks is None:
        networks = list(df[network_col].dropna().unique())

    mean_mat = (
        summary.pivot(index=network_col, columns=task_col, values="mean")
        .reindex(index=networks, columns=tasks)
    )

    std_mat = (
        summary.pivot(index=network_col, columns=task_col, values="std")
        .reindex(index=networks, columns=tasks)
    )

    n_mat = (
        summary.pivot(index=network_col, columns=task_col, values="n")
        .reindex(index=networks, columns=tasks)
    )

    return summary, mean_mat, std_mat, n_mat, networks
