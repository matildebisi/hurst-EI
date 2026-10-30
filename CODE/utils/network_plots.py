"""Network-level radar charts and task comparisons. Functions preserve historical plotting defaults."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms



def plot_network_task_radarplot_from_values(values,
                                            std_values,
                                            networks,
                                            tasks,
                                            task_colors,
                                            figsize=(4, 4),
                                            ymin=0.9,
                                            ymax=1.5,
                                            yticks=None,
                                            std_alpha=0.15,
                                            LABELSIZE = 20,
                                            TICKSIZE = 10,
                                            save_path=None,
                                            show=True):

    if len(values) != len(tasks):
        raise ValueError("values and tasks must have the same length.")

    if len(std_values) != len(tasks):
        raise ValueError("std_values and tasks must have the same length.")

    if len(task_colors) != len(tasks):
        raise ValueError("task_colors and tasks must have the same length.")

    if any(len(v) != len(networks) for v in values):
        raise ValueError(
            f"Each values entry must contain {len(networks)} network values."
        )

    if any(len(v) != len(networks) for v in std_values):
        raise ValueError(
            f"Each std_values entry must contain {len(networks)} network values."
        )

    N = len(networks)

    # Angles for each network
    angles = np.linspace(
        0,
        2 * np.pi,
        N,
        endpoint=False
    ).tolist()

    # Close the polygon
    angles_closed = angles + [angles[0]]

    if yticks is None:
        yticks = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

    theta_label = np.deg2rad(25)
    rotation = -25

    fig, ax = plt.subplots(
        figsize=figsize,
        subplot_kw=dict(polar=True)
    )

    # ---------------------------------------------------------
    # Polar settings
    # ---------------------------------------------------------

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    ax.set_ylim(ymin, ymax)

    # ---------------------------------------------------------
    # Grid
    # ---------------------------------------------------------

    ax.set_yticks(yticks)
    ax.set_yticklabels([""] * len(yticks))

    ax.yaxis.grid(
        True,
        color="#bdbdbd",
        linewidth=0.9,
        alpha=0.9,
        zorder=0
    )

    ax.xaxis.grid(
        True,
        color="#cfcfcf",
        linewidth=0.9,
        alpha=0.9,
        zorder=0
    )

    ax.spines["polar"].set_visible(False)

    # ---------------------------------------------------------
    # Network labels
    # ---------------------------------------------------------

    ax.set_xticks(angles)

    ax.set_xticklabels(
        networks,
        fontsize=LABELSIZE,
        color="black",
        fontweight="bold"
    )

    ax.tick_params(
        axis="x",
        pad=0
    )

    # Manual adjustment for CONT
    # if "CONT" in networks:

    #     label = ax.get_xticklabels()[
    #         networks.index("CONT")
    #     ]

    #     offset = mtransforms.ScaledTranslation(
    #         -0.4,
    #         0,
    #         ax.figure.dpi_scale_trans
    #     )

    #     label.set_transform(
    #         label.get_transform() + offset
    #     )

    # # Manual adjustment for DAN
    # if "DAN" in networks:

    #     label = ax.get_xticklabels()[
    #         networks.index("DAN")
    #     ]

    #     offset = mtransforms.ScaledTranslation(
    #         +0.3,
    #         0,
    #         ax.figure.dpi_scale_trans
    #     )

    #     label.set_transform(
    #         label.get_transform() + offset
    #     )

    # ---------------------------------------------------------
    # WHITE BACKGROUND OF RADIAL LABELS
    #
    # Important:
    # zorder = 0.5
    #
    # therefore:
    # grid < white background < STD bands < mean lines
    # ---------------------------------------------------------

    for r in yticks:

        txt = f"{r:.1f}"

        ax.text(
            theta_label,
            r,
            txt,

            # invisible text here:
            # we only want the white bbox
            color="none",

            fontsize=TICKSIZE,

            rotation=rotation,
            rotation_mode="anchor",

            ha="center",
            va="center",

            zorder=0.5,

            bbox=dict(
                facecolor="white",
                edgecolor="none",
                boxstyle="round,pad=0.12"
            )
        )

    # ---------------------------------------------------------
    # Plot tasks
    # ---------------------------------------------------------

    for task, task_values, task_std, color in zip(
        tasks,
        values,
        std_values,
        task_colors
    ):

        means = np.asarray(
            task_values,
            dtype=float
        )

        stds = np.asarray(
            task_std,
            dtype=float
        )

        lower = means - stds
        upper = means + stds

        # Close polygons
        means_closed = np.r_[
            means,
            means[0]
        ]

        lower_closed = np.r_[
            lower,
            lower[0]
        ]

        upper_closed = np.r_[
            upper,
            upper[0]
        ]

        # -----------------------------------------------------
        # STD band
        # -----------------------------------------------------

        ax.fill_between(
            angles_closed,
            lower_closed,
            upper_closed,

            color=color,
            alpha=std_alpha,

            linewidth=0,

            zorder=1
        )

        # -----------------------------------------------------
        # Mean line
        # -----------------------------------------------------

        ax.plot(
            angles_closed,
            means_closed,

            color=color,

            linewidth=2.5,

            marker="o",
            markersize=6,

            zorder=3,

            label=task
        )

    # ---------------------------------------------------------
    # RADIAL NUMBERS ON TOP
    #
    # Here there is NO bbox.
    # Only the gray numbers are redrawn.
    # ---------------------------------------------------------

    for r in yticks:

        txt = f"{r:.1f}"

        ax.text(
            theta_label,
            r,
            txt,

            fontsize=TICKSIZE,
            color="gray",

            rotation=rotation,
            rotation_mode="anchor",

            ha="center",
            va="center",

            zorder=20
        )

    # ---------------------------------------------------------
    # Layout
    # ---------------------------------------------------------

    plt.tight_layout()

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    if save_path is not None:

        plt.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
            facecolor="white"
        )

    # ---------------------------------------------------------
    # Show
    # ---------------------------------------------------------

    if show:
        plt.show()

    plt.close()



def plot_network_radarplot(mean_mat, std_mat, n_mat, networks, tasks,
                                         figsize=(4, 4),
                                         task_colors=None,
                                         ci_alpha=0.2,
                                         ymin=0.9,
                                         ymax=1.5,
                                         labelsize=12,
                                         ticksize=10,
                                         titlesize=14,
                                         save_path=None,
                                         show=True):

    if task_colors is None:
        task_colors = ['#ff5e69', '#86ab00', '#03b6ce', '#8491ff']

    N = len(networks)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles_closed = angles + [angles[0]]

    yticks = np.linspace(ymin, ymax, 6)[1:]
    theta_label = np.deg2rad(25)
    rotation = -25

    for task, color in zip(tasks, task_colors):
        means = mean_mat[task].to_numpy(dtype=float).tolist()
        stds  = std_mat[task].to_numpy(dtype=float).tolist()
        ns    = n_mat[task].to_numpy(dtype=float).tolist()

        ci = [1.96 * s / np.sqrt(n) for s, n in zip(stds, ns)]
        means_c  = means + means[:1]
        ci_upper = [max(m + c, ymin) for m, c in zip(means, ci)] + [max(means[0] + ci[0], ymin)]
        ci_lower = [max(m - c, ymin) for m, c in zip(means, ci)] + [max(means[0] - ci[0], ymin)]

        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_ylim(ymin, ymax)

        # grid
        ax.set_yticks(yticks)
        ax.set_yticklabels([""] * len(yticks))
        ax.yaxis.grid(True, color="#bdbdbd", linewidth=0.9, alpha=0.9, zorder=0)
        ax.xaxis.grid(True, color="#cfcfcf", linewidth=0.9, alpha=0.9, zorder=0)
        ax.spines["polar"].set_visible(False)

        ax.set_xticks(angles)
        ax.set_xticklabels(["LIM", "SMN", "DAN", "VAN", "VIS", "CONT", "DMN"],
                   fontsize=labelsize, color="black", fontweight="bold")
        ax.tick_params(axis='x', pad=0)

        # ticks
        for r in yticks:
            txt = f"{r:.1f}"
            ax.text(theta_label, r, txt,
                    fontsize=ticksize, color="gray",
                    rotation=rotation, rotation_mode="anchor",
                    ha="center", va="center", zorder=2,
                    bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.12"))
            
        ax.plot(angles_closed, means_c, color=color, linewidth=2.5,
                marker='o', markersize=6, zorder=3)
        ax.fill(angles_closed, ci_upper, color=color, alpha=ci_alpha, linewidth=0, zorder=2)
        for r in yticks:
            txt = f"{r:.1f}"
            ax.text(theta_label, r, txt,
                    fontsize=ticksize, color="gray",
                    rotation=rotation, rotation_mode="anchor",
                    ha="center", va="center", zorder=4)

        ax.set_title(f"{task}", fontsize=titlesize, fontweight='bold', pad=20)

        plt.tight_layout()

        if save_path is not None:
            path = save_path if len(tasks) == 1 else f"{save_path}/{task}.png"
            plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")

        if show:
            plt.show()
        plt.close()


def plot_network_all_task_radarplot(mean_mat,
                                std_mat,
                                n_mat,
                                networks,
                                tasks,
                                figsize=(4, 4),
                                labelsize=12,
                                ticksize=10,
                                task_colors=None,
                                ci_alpha=0.15,
                                ymin=0.9,
                                ymax=1.5,
                                save_path=None,
                                show=True):

    if task_colors is None:
        task_colors = [
            '#ff5e69',
            '#86ab00',
            '#03b6ce',
            '#8491ff'
        ]

    N = len(networks)

    angles = np.linspace(
        0,
        2 * np.pi,
        N,
        endpoint=False
    ).tolist()

    angles_closed = angles + [angles[0]]

    yticks = [1.0, 1.1, 1.2, 1.3, 1.4]#, 1.5]

    theta_label = np.deg2rad(25)
    rotation = -25

    fig, ax = plt.subplots(
        figsize=figsize,
        subplot_kw=dict(polar=True)
    )

    # ---------------------------------------------------------
    # Polar settings
    # ---------------------------------------------------------

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    ax.set_ylim(ymin, ymax)

    # ---------------------------------------------------------
    # Grid
    # ---------------------------------------------------------

    ax.set_yticks(yticks)
    ax.set_yticklabels([""] * len(yticks))

    ax.yaxis.grid(
        True,
        color="#bdbdbd",
        linewidth=0.9,
        alpha=0.9,
        zorder=0
    )

    ax.xaxis.grid(
        True,
        color="#cfcfcf",
        linewidth=0.9,
        alpha=0.9,
        zorder=0
    )

    ax.spines["polar"].set_visible(False)

    # ---------------------------------------------------------
    # Network labels
    # ---------------------------------------------------------

    labels = [
        "LIM",
        "SMN",
        "DAN",
        "VAN",
        "VIS",
        "CONT",
        "DMN"
    ]

    ax.set_xticks(angles)

    ax.set_xticklabels(
        labels,
        fontsize=labelsize,
        color="black",
        fontweight="bold"
    )

    ax.tick_params(
        axis="x",
        pad=0
    )

    # Move CONT slightly
    cont_label = ax.get_xticklabels()[
        labels.index("CONT")
    ]

    offset = mtransforms.ScaledTranslation(
        -0.4,
        0,
        ax.figure.dpi_scale_trans
    )

    cont_label.set_transform(
        cont_label.get_transform() + offset
    )

    # Move DAN slightly
    dan_label = ax.get_xticklabels()[
        labels.index("DAN")
    ]

    offset = mtransforms.ScaledTranslation(
        +0.3,
        0,
        ax.figure.dpi_scale_trans
    )

    dan_label.set_transform(
        dan_label.get_transform() + offset
    )

    # ---------------------------------------------------------
    # White background behind radial labels
    #
    # It sits BELOW the confidence bands.
    # ---------------------------------------------------------

    for r in yticks:

        txt = f"{r:.1f}"

        ax.text(
            theta_label,
            r,
            txt,
            fontsize=ticksize,

            # invisible text; only bbox is needed
            color="none",

            rotation=rotation,
            rotation_mode="anchor",

            ha="center",
            va="center",

            zorder=0.5,

            bbox=dict(
                facecolor="white",
                edgecolor="none",
                boxstyle="round,pad=0.12"
            )
        )

    # ---------------------------------------------------------
    # Plot tasks
    # ---------------------------------------------------------

    for task, color in zip(tasks, task_colors):

        means = mean_mat[task].to_numpy(dtype=float)
        stds = std_mat[task].to_numpy(dtype=float)
        ns = n_mat[task].to_numpy(dtype=float)

        # -----------------------------------------------------
        # 95% confidence interval
        #
        # CI = mean ± 1.96 * SD / sqrt(n)
        # -----------------------------------------------------

        sem = stds / np.sqrt(ns)
        ci = 1.96 * sem

        ci_lower = means - ci
        ci_upper = means + ci

        # Optional clipping to plotting limits
        ci_lower = np.maximum(ci_lower, ymin)
        ci_upper = np.minimum(ci_upper, ymax)

        # Close radar
        means_closed = np.r_[
            means,
            means[0]
        ]

        ci_lower_closed = np.r_[
            ci_lower,
            ci_lower[0]
        ]

        ci_upper_closed = np.r_[
            ci_upper,
            ci_upper[0]
        ]

        # -----------------------------------------------------
        # 95% CI continuous band
        # -----------------------------------------------------

        ax.fill_between(
            angles_closed,
            ci_lower_closed,
            ci_upper_closed,
            color=color,
            alpha=ci_alpha,
            linewidth=0,
            zorder=1
        )

        # -----------------------------------------------------
        # Mean
        # -----------------------------------------------------

        ax.plot(
            angles_closed,
            means_closed,
            color=color,
            linewidth=2.5,
            marker="o",
            markersize=6,
            zorder=3,
            label=task
        )

    # ---------------------------------------------------------
    # Draw radial numbers ABOVE bands and lines
    # ---------------------------------------------------------

    for r in yticks:

        txt = f"{r:.1f}"

        ax.text(
            theta_label,
            r,
            txt,
            fontsize=ticksize,
            color="gray",

            rotation=rotation,
            rotation_mode="anchor",

            ha="center",
            va="center",

            zorder=20
        )

    # ---------------------------------------------------------
    # Layout
    # ---------------------------------------------------------

    plt.tight_layout()

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    if save_path is not None:

        plt.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
            facecolor="white"
        )

    # ---------------------------------------------------------
    # Show
    # ---------------------------------------------------------

    if show:
        plt.show()

    plt.close()