import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def dotplot(
    df_plot,
    x_order=None,
    y_order=None,
    x_col='x',
    y_col='y',
    xticklabels=None,
    yticklabels=None,
    pval_col='p',
    diff_col='logFC',
    figsize=(18, 4.5),
    min_dot_size=20,
    max_dot_size=300,
    cmap="RdBu_r",
    cbar_title='',
    vmin=None,
    vmax=None,
    sig_edgecolor="black",
    sig_edgewidth=1.5,
    nonsig_edgecolor="0.7",
    nonsig_edgewidth=0.5,
    bg_color="white",
    title_fontsize=12,
    axis_label_fontsize=10,
    tick_fontsize=9,
    celltype_rotation=45,
    legend_fontsize=10,
    fdr_thresh=0.05,
    save_path=None,
    save_dpi=300,
    bbox_to_anchor_dot = (1.05, 1),
    ax=None
    ):
    """
    Plot a bubble chart for cell type differences across y_order.

    Parameters:
    - df_plot: DataFrame containing the data to plot.
    - y_order: List of y axis category names.
    - x_order: List of x axis category names.
    - figsize: Figure size.
    - min_dot_size, max_dot_size: Minimum and maximum bubble sizes.
    - cmap: Colormap for the bubbles.
    - vmin, vmax: Value range for the colormap.
    - sig_edgecolor, nonsig_edgecolor: Edge colors for significant and non-significant bubbles.
    - sig_edgewidth, nonsig_edgewidth: Edge widths for significant and non-significant bubbles.
    - bg_color: Background color.
    - title_fontsize, axis_label_fontsize, tick_fontsize, celltype_fontsize: Font sizes.
    - celltype_rotation: Rotation angle for cell type labels.
    - legend_fontsize: Font size for the legend.
    - fdr_thresh: FDR threshold for significance.
    - save_path: Path to save the figure.
    - save_dpi: DPI for saving the figure.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor(bg_color)

    # 1. Setup Ordering
    if x_order is None: x_order = sorted(df_plot[x_col].unique())
    if y_order is None: y_order = sorted(df_plot[y_col].unique())

    x_map = {name: i for i, name in enumerate(x_order)}
    y_map = {name: i for i, name in enumerate(y_order)}

    # 2. Prepare Data (Vectorized)
    df = df_plot.copy()
    df["neglog10fdr"] = -np.log10(df[pval_col].clip(1e-300, 1))

    # Scaling logic
    size_min, size_max = df["neglog10fdr"].min(), df["neglog10fdr"].max()
    if size_max == size_min: size_max = size_min + 1

    # Calculate sizes in points squared (for scatter)
    df["size"] = min_dot_size + (
        (df["neglog10fdr"] - size_min) / (size_max - size_min)
    ) * (max_dot_size - min_dot_size)

    # Significance logic: p-value < thresh  =>  -log10(p) > -log10(thresh)
    is_sig = df["neglog10fdr"] > -np.log10(fdr_thresh)

    # Map coordinates
    x_coords = df[x_col].map(x_map)
    y_coords = df[y_col].map(y_map)

    # Auto color range
    if vmin is None or vmax is None:
        abs_max = df[diff_col].abs().quantile(0.95)
        vmin, vmax = -abs_max, abs_max

    # 3. Plotting
    sc = ax.scatter(
        x_coords, y_coords,
        s=df["size"],
        c=df[diff_col], cmap=cmap, vmin=vmin, vmax=vmax,
        edgecolors=np.where(is_sig, sig_edgecolor, nonsig_edgecolor),
        linewidths=np.where(is_sig, sig_edgewidth, nonsig_edgewidth),
        zorder=3
    )

    if xticklabels is False:
        ax.set_xticks([])
    else:
        labels = x_order if xticklabels is None else xticklabels
        ticks = range(len(labels))
        ax.set_xticks(ticks)
        ax.set_xticklabels(
            labels,
            rotation=celltype_rotation,
            ha="right",
            fontsize=tick_fontsize
        )

    # Y axis
    if yticklabels is False:
        ax.set_yticks([])
    else:
        labels = y_order if yticklabels is None else yticklabels
        ticks = range(len(labels))
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, fontsize=tick_fontsize)

    ax.set_xlim(-0.5, len(x_order) - 0.5)
    ax.set_ylim(-0.5, len(y_order) - 0.5)

    # Colorbar
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label(cbar_title, fontsize=axis_label_fontsize)

    # 5. Legend
    legend_vals = np.linspace(size_min, size_max, 3)
    # Convert back from points^2 area to diameter for markersize
    size_handles = [
        Line2D([0], [0], marker="o", linestyle="None", color="gray",
               markersize=np.sqrt(min_dot_size + ((v - size_min) / (size_max - size_min)) * (max_dot_size - min_dot_size)) * 0.8,
               label=f"{v:.1f}")
        for v in legend_vals
    ]

    sig_handle = Line2D([0], [0], marker="o", linestyle="None", color="white",
                        markeredgecolor="black", markeredgewidth=1.5,
                        markersize=10, label=f"FDR < {fdr_thresh}")

    ax.legend(handles=size_handles + [sig_handle], title="-log₁₀(FDR)",
              title_fontsize=legend_fontsize, fontsize=legend_fontsize,
              loc="upper left", bbox_to_anchor=bbox_to_anchor_dot, frameon=False)

    plt.tight_layout()
    if save_path: fig.savefig(save_path, dpi=save_dpi, bbox_inches="tight")

    return ax
