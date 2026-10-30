"""Surface mesh handling and parcelwise cortical visualizations."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import nibabel as nib
import scipy.sparse as sp

from scipy.sparse.csgraph import connected_components
from matplotlib.colors import LinearSegmentedColormap
from nilearn import datasets, plotting, surface
from nilearn.surface import load_surf_mesh
from nibabel.freesurfer.io import read_annot


# Published cortex map order: lateral L, medial L, medial R, lateral R.
# This tuple is used by both Schaefer parcellated and voxel-to-surface maps.
SURFACE_VIEWS = (
    ("left", "lateral"),
    ("left", "medial"),
    ("right", "medial"),
    ("right", "lateral"),
)


def plot_cortical_map(parcel_values, *, left_annotation, right_annotation,
                      output_path, **style):
    """Plot a parcellated cortex in the four canonical surface views.

    Values and atlas ordering are not transformed; input index ordering must
    match the original Schaefer/Glasser annotation used to estimate H.
    """
    return plot_parcelwise_with_roi_surface_direct(
        parcel_values, left_annotation, right_annotation,
        img_name=str(output_path), **style
    )

def make_partially_inflated_mesh(pial_mesh, inflated_mesh, inflate_fraction=0.20):
    pial_coords, pial_faces = load_surf_mesh(pial_mesh)
    infl_coords, infl_faces = load_surf_mesh(inflated_mesh)

    if pial_coords.shape != infl_coords.shape:
        raise ValueError("Pial and inflated surfaces have different coordinate shapes.")

    if not np.array_equal(pial_faces, infl_faces):
        raise ValueError("Pial and inflated surfaces have different face topology.")

    coords = (1.0 - inflate_fraction) * pial_coords + inflate_fraction * infl_coords
    return coords, pial_faces


def get_mesh_and_bg(fsaverage, hemi, inflate_fraction=0.20):
    if hemi == "left":
        pial = fsaverage.pial_left
        inflated = fsaverage.infl_left
        bg = fsaverage.sulc_left
    elif hemi == "right":
        pial = fsaverage.pial_right
        inflated = fsaverage.infl_right
        bg = fsaverage.sulc_right
    else:
        raise ValueError("hemi must be 'left' or 'right'")

    mesh = make_partially_inflated_mesh(
        pial,
        inflated,
        inflate_fraction=inflate_fraction
    )

    return mesh, bg


def load_surface_annotation(annot_path):
    labels, ctab, names = read_annot(annot_path, orig_ids=False)
    names = [n.decode("utf-8") if isinstance(n, bytes) else str(n) for n in names]
    return labels.astype(int), names


def is_real_parcel_name(name):
    name_low = name.lower()
    bad = ["unknown", "medial", "wall", "background", "???"]
    return not any(b in name_low for b in bad)


def build_surface_map_from_matrix(labels, names, matrix, offset=0):
    surf_map = np.full(labels.shape, np.nan, dtype=float)

    parcel_label_ids = [
        idx for idx, name in enumerate(names)
        if is_real_parcel_name(name)
    ]

    for local_idx, label_id in enumerate(parcel_label_ids):
        matrix_idx = offset + local_idx

        if matrix_idx >= len(matrix):
            continue

        surf_map[labels == label_id] = matrix[matrix_idx]

    return surf_map, len(parcel_label_ids)


def make_roi_noholes_from_surface_mask(mask, surf_mesh):
    coords, faces = load_surf_mesh(surf_mesh)
    n = coords.shape[0]

    i = np.hstack([
        faces[:, 0], faces[:, 0],
        faces[:, 1], faces[:, 1],
        faces[:, 2], faces[:, 2]
    ])
    j = np.hstack([
        faces[:, 1], faces[:, 2],
        faces[:, 0], faces[:, 2],
        faces[:, 0], faces[:, 1]
    ])

    A = sp.coo_matrix((np.ones_like(i), (i, j)), shape=(n, n)).tocsr()
    A = A + A.T
    A.setdiag(0)
    A.eliminate_zeros()

    idx = np.where(mask)[0]

    if idx.size > 0:
        ncc, labels_cc = connected_components(A[idx][:, idx], directed=False)
        keep = np.bincount(labels_cc).argmax()
        keep_idx = idx[labels_cc == keep]

        mask_big = np.zeros(n, dtype=bool)
        mask_big[keep_idx] = True
    else:
        mask_big = np.zeros(n, dtype=bool)

    bg = ~mask_big
    idx_bg = np.where(bg)[0]

    if idx_bg.size > 0:
        ncc_bg, labels_bg = connected_components(A[idx_bg][:, idx_bg], directed=False)
        keep_bg = np.bincount(labels_bg).argmax()
        outside_idx = idx_bg[labels_bg == keep_bg]

        filled_mask = np.ones(n, dtype=bool)
        filled_mask[outside_idx] = False
    else:
        filled_mask = mask_big

    return filled_mask.astype(int)


def plot_parcelwise_with_roi_surface_direct(
    matrix,
    lh_annot_path,
    rh_annot_path,
    roi_img_path=None,
    img_name="parcelwise",
    show=False,
    max_color=None,
    colormap=None,
    mesh="fsaverage6",
    inflate_fraction=0.20,
    vmin=None,
    vmax=None,
    label=None,
    label_fs=18,
    tick_fs=18,
    dpi=150
):

    matrix = np.asarray(matrix, dtype=float)

    if vmin is None:
        vmin = np.nanmin(matrix)

    if vmax is None:
        vmax = np.nanmax(matrix)

    fsaverage = datasets.fetch_surf_fsaverage(mesh=mesh)

    labels_lh, names_lh = load_surface_annotation(lh_annot_path)
    labels_rh, names_rh = load_surface_annotation(rh_annot_path)

    mesh_lh, bg_lh = get_mesh_and_bg(
        fsaverage,
        "left",
        inflate_fraction=inflate_fraction
    )

    mesh_rh, bg_rh = get_mesh_and_bg(
        fsaverage,
        "right",
        inflate_fraction=inflate_fraction
    )

    surf_map_lh, n_lh = build_surface_map_from_matrix(
        labels_lh,
        names_lh,
        matrix,
        offset=0
    )

    surf_map_rh, n_rh = build_surface_map_from_matrix(
        labels_rh,
        names_rh,
        matrix,
        offset=n_lh
    )

    colors = ["#FFFFFF", max_color if max_color is not None else "#7F0000"]
    cmap_custom = LinearSegmentedColormap.from_list(
        "custom",
        colors,
        N=256
    )
    if colormap is not None:
        cmap_custom = colormap

    roi_noholes_lh = None

    if roi_img_path is not None:
        roi_img = nib.load(roi_img_path)

        roi_surf_lh = surface.vol_to_surf(
            roi_img,
            fsaverage.pial_left,
            interpolation="linear",
            radius=3
        )

        eps = 1
        roi_mask_lh = roi_surf_lh > eps

        roi_noholes_lh = make_roi_noholes_from_surface_mask(
            roi_mask_lh,
            mesh_lh
        )

    views = SURFACE_VIEWS

    fig, axes = plt.subplots(
        1,
        len(views),
        figsize=(10 * len(views), 10),
        subplot_kw={"projection": "3d"}
    )

    axes = np.array(axes).flatten()

    fig.subplots_adjust(
        left=0.02,
        right=0.90,
        wspace=0.00
    )

    for j, (hemi, view) in enumerate(views):
        ax = axes[j]

        if hemi == "left":
            surf_mesh = mesh_lh
            surf_map = surf_map_lh
            bg_map = bg_lh
        else:
            surf_mesh = mesh_rh
            surf_map = surf_map_rh
            bg_map = bg_rh

        plotting.plot_surf_stat_map(
            surf_mesh,
            surf_map,
            hemi=hemi,
            view=view,
            cmap=cmap_custom,
            colorbar=False,
            bg_on_data=True,
            axes=ax,
            figure=fig,
            vmin=vmin,
            vmax=vmax,
            bg_map=bg_map,
           # darkness=0.8,
            avg_method="mean",
        )

        if hemi == "left" and roi_noholes_lh is not None:
            plotting.plot_surf_contours(
                mesh_lh,
                roi_noholes_lh,
                hemi=hemi,
                levels=[1],
                colors=["white"],
                linewidths=3,
                axes=ax,
                figure=fig,
            )

        ax.axis("off")

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap_custom, norm=norm)
    sm.set_array([])

    fig.colorbar(
        sm,
        ax=axes,
        shrink=0.3,
        aspect=10,
        location="right",
        label=label,
        ticks=[vmin, vmax]
    )

    cbar_ax = fig.axes[-1]
    cbar_ax.tick_params(labelsize=tick_fs)
    cbar_ax.set_ylabel(
        label,
        fontsize=label_fs,
        fontweight="bold"
    )

    fig.savefig(
        img_name,
        format="png",
        dpi=dpi,
        bbox_inches="tight"
    )

    if show:
        plt.show()
    else:
        plt.close(fig)




