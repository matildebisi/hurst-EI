"""Shared outlier filtering and parcellation helpers. Server-side parcellation is unchanged."""

import numpy as np
import pandas as pd
import nibabel as nib
import scipy.io
import os


def remove_outliers_IQR_from_df(df, threshold=1.5, value_col="Hurst", group_cols=None):
    """Remove IQR-based outliers within each group (default: Parcel x Task, or Parcel only)."""
    if group_cols is None:
        group_cols = ["Parcel", "Task"] if "Task" in df.columns else ["Parcel"]

    grouped = df.groupby(group_cols)[value_col]
    q1 = grouped.transform("quantile", q=0.25)
    q3 = grouped.transform("quantile", q=0.75)
    n = grouped.transform("count")
    iqr = q3 - q1

    within_bounds = (df[value_col] >= q1 - threshold * iqr) & (df[value_col] <= q3 + threshold * iqr)
    keep = within_bounds | (n <= 1) | (iqr == 0) | df[value_col].isna()

    return df.loc[keep].reset_index(drop=True)


def parcellate(
    atlasfile,
    datafile,
    fname2save,
    mean_center=True,
    expected_nreg=None
):
    """
    Extract mean time series from a parcellation atlas.

    Parcel labels are read directly from the atlas, so they do NOT
    need to be consecutive.

    Parameters
    ----------
    atlasfile : str
        Path to atlas NIfTI file.
    datafile : str
        Path to 4D fMRI/other data NIfTI file.
    fname2save : str
        Path to output MAT file.
    mean_center : bool
        If True, mean-center each parcel's time series.
    expected_nreg : int or None
        Optional expected number of parcels, used as a sanity check.

    Returns
    -------
    result : numpy.ndarray
        Array of shape (n_volumes, n_regions).
    """

    # -------------------------------------------------------------
    # Load atlas
    # -------------------------------------------------------------

    atlas_img = nib.load(atlasfile)
    atlas_float = atlas_img.get_fdata()

    # Check that atlas really contains integer labels
    if not np.allclose(atlas_float, np.round(atlas_float)):
        raise ValueError(
            "Atlas contains non-integer values. "
            "It may have been resampled using linear interpolation "
            "instead of nearest-neighbour interpolation."
        )

    atlas = np.round(atlas_float).astype(int)

    # -------------------------------------------------------------
    # Get REAL parcel labels from atlas
    # -------------------------------------------------------------

    parcel_labels = np.unique(atlas)

    # Remove background
    parcel_labels = parcel_labels[parcel_labels != 0]

    # Sorted automatically:
    # 1 ... 180, 1001 ... 1180
    parcel_labels = np.sort(parcel_labels)

    print(f"Atlas contains {len(parcel_labels)} parcels.")
    print("First labels:", parcel_labels[:10])
    print("Hemisphere transition:", parcel_labels[175:185])
    print("Last labels:", parcel_labels[-10:])

    if expected_nreg is not None and len(parcel_labels) != expected_nreg:
        raise ValueError(
            f"Expected {expected_nreg} parcels, "
            f"but atlas contains {len(parcel_labels)} non-zero labels."
        )

    # -------------------------------------------------------------
    # Load fMRI data
    # -------------------------------------------------------------

    data_img = nib.load(datafile)
    data = data_img.get_fdata()

    if data.ndim != 4:
        raise ValueError(
            f"Data must be 4D, got shape {data.shape}"
        )

    if atlas.shape != data.shape[:3]:
        raise ValueError(
            f"Atlas shape {atlas.shape} does not match "
            f"data shape {data.shape[:3]}"
        )

    # -------------------------------------------------------------
    # Parcellate
    # -------------------------------------------------------------

    nvols = data.shape[3]
    nparcels = len(parcel_labels)

    result = np.full(
        (nvols, nparcels),
        np.nan,
        dtype=float
    )

    for i, parcel_label in enumerate(parcel_labels):

        mask = atlas == parcel_label

        n_voxels = np.sum(mask)

        if n_voxels == 0:
            print(
                f"[WARNING] Parcel {parcel_label} has no voxels."
            )
            continue

        # data[mask] -> shape:
        # n_voxels x n_volumes
        parcel_data = data[mask]

        result[:, i] = np.mean(
            parcel_data,
            axis=0
        )

    # -------------------------------------------------------------
    # Mean center each parcel
    # -------------------------------------------------------------

    if mean_center:
        result = result - np.nanmean(
            result,
            axis=0,
            keepdims=True
        )

    # -------------------------------------------------------------
    # Sanity checks
    # -------------------------------------------------------------

    all_nan = np.all(np.isnan(result), axis=0)

    print(
        f"Valid parcels: {np.sum(~all_nan)} / {nparcels}"
    )

    if np.any(all_nan):
        print(
            "Parcels containing only NaN:",
            parcel_labels[all_nan]
        )

    # Special Glasser check
    if nparcels == 360:

        lh_valid = np.sum(
            ~np.all(np.isnan(result[:, :180]), axis=0)
        )

        rh_valid = np.sum(
            ~np.all(np.isnan(result[:, 180:]), axis=0)
        )

        print(f"LH valid parcels: {lh_valid}/180")
        print(f"RH valid parcels: {rh_valid}/180")

    # -------------------------------------------------------------
    # Save
    # -------------------------------------------------------------

    output_folder = os.path.dirname(fname2save)

    if output_folder:
        os.makedirs(
            output_folder,
            exist_ok=True
        )

    scipy.io.savemat(
        fname2save,
        {
            "avg_time_series": result,

            # tells you which atlas label corresponds
            # to each column of avg_time_series
            "parcel_labels": parcel_labels
        }
    )

    return result