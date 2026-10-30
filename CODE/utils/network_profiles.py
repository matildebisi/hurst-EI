"""Parcel-to-network mapping and component profile aggregations.

Keep Schaefer and Glasser labeling rules explicit; the profile functions do
not orient SVD component signs unless requested by the caller.
"""

import numpy as np
import pandas as pd

NETWORK_ORDER = ["Limbic", "SomMot", "DorsAttn", "SalVentAttn", "Vis", "Cont", "Default"]


def load_schaefer_networks(csv_path):
    df = pd.read_csv(csv_path, usecols=[0, 1], engine="python")
    roi_name_col = df.columns[1]

    def extract_network(roi_name):
        parts = str(roi_name).split("_")
        if len(parts) >= 3:
            return parts[2]
        return None

    return df[roi_name_col].apply(extract_network).tolist()


def orient_component_sign(vec, method="largest_abs_positive"):
    vec = np.asarray(vec, dtype=float).copy()

    if method == "largest_abs_positive":
        idx = np.nanargmax(np.abs(vec))
        if vec[idx] < 0:
            vec *= -1

    elif method == "sum_positive":
        if np.nansum(vec) < 0:
            vec *= -1

    return vec


def compute_network_mean_signed_profiles(loadings_vec, parcel_networks, network_order=NETWORK_ORDER):
    """
    For each network:
      1) separate parcels with positive and negative loadings
      2) average positive and absolute negative loadings, respectively
    """
    df = pd.DataFrame({
        "loading": np.asarray(loadings_vec, dtype=float),
        "network": parcel_networks
    })

    pos_profile = []
    neg_profile = []

    for net in network_order:
        vals = df.loc[df["network"] == net, "loading"].to_numpy(dtype=float)

        pos_vals = vals[vals > 0]
        neg_vals = vals[vals < 0]

        pos_profile.append(np.mean(pos_vals) if len(pos_vals) > 0 else 0.0)
        neg_profile.append(np.abs(np.mean(neg_vals)) if len(neg_vals) > 0 else 0.0)

    return np.array(pos_profile), np.array(neg_profile)


def compute_network_sign_fraction_profiles(loadings_vec, parcel_networks, network_order=NETWORK_ORDER):
    df = pd.DataFrame({
        "loading": np.asarray(loadings_vec, dtype=float),
        "network": parcel_networks
    })

    pos_frac = []
    neg_frac = []

    for net in network_order:
        vals = df.loc[df["network"] == net, "loading"].to_numpy(dtype=float)

        if len(vals) == 0:
            pos_frac.append(0.0)
            neg_frac.append(0.0)
            continue

        pos_frac.append(np.mean(vals > 0))
        neg_frac.append(np.mean(vals < 0))

    return np.array(pos_frac), np.array(neg_frac)


def compute_network_mean_profile(loadings_vec, parcel_networks, network_order=NETWORK_ORDER):
    df = pd.DataFrame({
        "loading": np.asarray(loadings_vec, dtype=float),
        "network": parcel_networks
    })

    profile = []

    for net in network_order:
        vals = df.loc[df["network"] == net, "loading"].to_numpy(dtype=float)
        profile.append(np.nanmean(vals) if len(vals) > 0 else np.nan)

    return np.array(profile)


def normalize_network_name(net):
    """
    Normalize possible network naming conventions to the
    names used in NETWORK_ORDER.
    """

    net = str(net).strip()

    mapping = {
        # abbreviations
        "LIM": "Limbic",
        "SMN": "SomMot",
        "DAN": "DorsAttn",
        "VAN": "SalVentAttn",
        "VIS": "Vis",
        "CONT": "Cont",
        "DMN": "Default",

        # possible full names
        "Limbic": "Limbic",
        "SomMot": "SomMot",
        "Somatomotor": "SomMot",
        "DorsAttn": "DorsAttn",
        "DorsalAttention": "DorsAttn",
        "Dorsal Attention": "DorsAttn",
        "SalVentAttn": "SalVentAttn",
        "VentralAttention": "SalVentAttn",
        "Ventral Attention": "SalVentAttn",
        "Vis": "Vis",
        "Visual": "Vis",
        "Cont": "Cont",
        "Control": "Cont",
        "Frontoparietal": "Cont",
        "Default": "Default",
        "DefaultMode": "Default",
        "Default Mode": "Default",
    }

    return mapping.get(net, net)


def load_glasser_networks(network_info_df):
    """
    Returns a list of network assignments ordered by Glasser
    parcel ID.

    Expected network_info_df columns:
        Parcel
        Network

    Expected parcels:
        1 ... 360
    """

    df = network_info_df[["Parcel", "Network"]].copy()

    df["Parcel"] = pd.to_numeric(
        df["Parcel"],
        errors="coerce"
    )

    df = df.dropna(
        subset=["Parcel", "Network"]
    )

    df["Parcel"] = df["Parcel"].astype(int)

    # One network assignment per parcel
    if df["Parcel"].duplicated().any():
        duplicated = df.loc[
            df["Parcel"].duplicated(keep=False),
            "Parcel"
        ].unique()

        raise ValueError(
            f"Duplicated parcel IDs in network_info_df: {duplicated}"
        )

    df = df.sort_values("Parcel")

    # Check Glasser360
    expected = np.arange(1, 361)
    observed = df["Parcel"].to_numpy()

    missing = np.setdiff1d(expected, observed)

    if len(missing) > 0:
        raise ValueError(
            f"Missing {len(missing)} Glasser parcels in network_info_df:\n"
            f"{missing}"
        )

    # Keep exactly parcels 1...360
    df = (
        df
        .set_index("Parcel")
        .reindex(expected)
    )

    parcel_networks = (
        df["Network"]
        .apply(normalize_network_name)
        .tolist()
    )

    print("\nGlasser network counts:")
    print(
        pd.Series(parcel_networks)
        .value_counts()
        .reindex(NETWORK_ORDER, fill_value=0)
    )

    return parcel_networks
