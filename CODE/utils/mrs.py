"""Reusable MRS measurements and paired-comparison helpers.

Numerical conventions match the original research notebooks; in particular the
Wilcoxon effect-size function rounds differences before excluding zero ranks.
"""

import numpy as np
import pandas as pd
from scipy import stats


def parse_real_numeric_array(x):
    out = []
    for val in x:
        if pd.isna(val):
            out.append(np.nan)
            continue

        if isinstance(val, (int, float, np.integer, np.floating)):
            out.append(float(val))
            continue

        s = str(val).strip()

        if s == "" or s.lower() in ["nan", "na", "none"]:
            out.append(np.nan)
            continue

        try:
            z = complex(s.replace("i", "j"))
            if abs(z.imag) > 1e-8:
                print(f"[WARNING] complex value with non-zero imaginary part: {s}")
            out.append(float(z.real))
        except Exception:
            out.append(pd.to_numeric(s, errors="coerce"))

    return np.asarray(out, dtype=float)


def normality_test(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    stat, p = stats.shapiro(x)
    return stat, p


def cohens_dz_paired(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    diff = x - y
    diff = diff[np.isfinite(diff)]

    if len(diff) < 2:
        return np.nan

    sd = np.std(diff, ddof=1)
    if sd == 0:
        return np.nan
    return np.mean(diff) / sd


def rank_biserial_from_wilcoxon(x, y):
    """
    Approximate matched-pairs rank-biserial correlation from Wilcoxon signed-rank test,
    using pingouin if available.
    """
    try:
        import pingouin as pg

        tmp = pg.wilcoxon(x, y)
        # RBC column exists in pingouin
        if "RBC" in tmp.columns:
            return float(tmp["RBC"].iloc[0])
    except Exception:
        pass
    return np.nan


def wilcoxon_effect_size(x, y, round_decimals=12):
    """
    Wilcoxon signed-rank test with effect size:

        r = |z| / sqrt(N)

    where N is the number of non-zero paired differences.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Remove invalid pairs
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    # Compute paired differences directly
    differences = np.round(x - y, decimals=round_decimals)

    # Consistent with zero_method="wilcox"
    n_effective = np.sum(differences != 0)

    if n_effective == 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    result = stats.wilcoxon(
        differences,
        zero_method="wilcox",
        alternative="two-sided",
        method="approx",
        correction=False
    )

    statistic = float(result.statistic)
    p_value = float(result.pvalue)
    z_value = float(result.zstatistic)

    effect_size_r = abs(z_value) / np.sqrt(n_effective)

    return statistic, p_value, z_value, effect_size_r, n_effective


def effect_size_magnitude(r):
    if np.isnan(r):
        return np.nan
    elif r < 0.10:
        return "negligible"
    elif r < 0.30:
        return "small"
    elif r < 0.50:
        return "medium"
    else:
        return "large"
