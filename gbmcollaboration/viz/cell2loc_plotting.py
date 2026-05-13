# + taken from https://cell2location.readthedocs.io/en/latest/_modules/cell2location/plt/plot_heatmap.html#clustermap
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance
from scipy.cluster.hierarchy import linkage, leaves_list

def heatmap(
    array,
    ticks=False,
    log=False,
    figsize=None,
    equal=False,
    row_labels=None,
    col_labels=None,
    cbar=True,
    cmap="RdPu",
    title="",
    vmin=None,
    vmax=None,
):
    r"""Plot heatmap with row and column labels using plt.imshow

    :param array: np.ndarray to be visualised, or an object that can be coerced to np.ndarray
    :param ticks: boolean, show x and y axis ticks?
    :param log: boolean, color on logscale?
    :param figsize: figure size as a tuple (x, y)
    :param equal: boolean, each tile should be square (equal aspect)
    :param row_labels: names of rows (pd.Series, pd.Index or list)
    :param col_labels: names of columns (pd.Series, pd.Index or list)
    :param cbar: boolean, show colorbar?
    :param cmap: valid matplotlib colormap name
    :param title: title of the plot
    """
    if figsize is not None:
        plt.figure(figsize=figsize)

    array = np.array(array)
    if log:
        plt.imshow(array, interpolation="nearest", cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
    else:
        plt.imshow(array, interpolation="nearest", cmap=cmap)

    if cbar is True:
        plt.colorbar()

    if ticks:
        plt.xticks(range(array.shape[1]))
        plt.yticks(range(array.shape[0]))
        plt.xlim(-0.5, array.shape[1] - 0.5)
        plt.ylim(array.shape[0] - 0.5, -0.5)

    if equal:
        plt.gca().set_aspect("equal", adjustable="box")

    if row_labels is not None:
        plt.yticks(range(array.shape[0]), row_labels)
        plt.ylim(array.shape[0] - 0.5, -0.5)

    if col_labels is not None:
        plt.xticks(range(array.shape[1]), col_labels, rotation=-90)
        plt.xlim(array.shape[1] - 0.5, -0.5)

    plt.title(title)
    plt.tight_layout()



def dotplot(
    array_color,
    array_size=None,
    ticks=False,
    log=False,
    figsize=None,
    equal=False,
    row_labels=None,
    col_labels=None,
    cbar=True,
    cmap="RdPu",
    title="",
):
    r"""Plot dotplot with row and column labels

    :param array_color: np.ndarray to be visualised as dot color
    :param array_size: np.ndarray to be visualised as dor size
    :param ticks: boolean, show x and y axis ticks?
    :param log: boolean, color on logscale?
    :param figsize: figure size as a tuple (x, y)
    :param equal: boolean, each tile should be square (equal aspect)
    :param row_labels: names of rows (pd.Series, pd.Index or list)
    :param col_labels: names of columns (pd.Series, pd.Index or list)
    :param cbar: boolean, show colorbar?
    :param cmap: valid matplotlib colormap name
    :param title: title of the plot
    """
    if figsize is not None:
        plt.figure(figsize=figsize)

    array_color = np.array(array_color)
    y, x = np.indices(array_color.shape)
    if array_size is None:
        array_size = array_color

    array_size = (array_size / array_size.max() * 15) ** 2

    if log:
        plt.scatter(
            x.flatten(),
            y.flatten(),
            c=array_color.flatten(),
            s=array_size.flatten(),
            cmap=cmap,
            norm=matplotlib.colors.LogNorm(),
            edgecolor="none",
        )
    else:
        plt.scatter(
            x.flatten(),
            y.flatten(),
            c=array_color.flatten(),
            s=array_size.flatten(),
            cmap=cmap,
            norm=None,
            edgecolor="none",
        )

    if cbar is True:
        plt.colorbar()

    if ticks:
        plt.xticks(range(array_color.shape[1]))
        plt.yticks(range(array_color.shape[0]))
        plt.xlim(-0.5, array_color.shape[1] - 0.5)
        plt.ylim(array_color.shape[0] - 0.5, -0.5)

    plt.grid(False)
    if equal:
        plt.gca().set_aspect("equal", adjustable="box")

    if row_labels is not None:
        plt.yticks(range(array_color.shape[0]), row_labels)
        plt.ylim(array_color.shape[0] - 0.5, -0.5)

    if col_labels is not None:
        plt.xticks(range(array_color.shape[1]), col_labels, rotation=-90)
        plt.xlim(array_color.shape[1] - 0.5, -0.5)

    plt.title(title)
    plt.tight_layout()



def clustermap(
    df,
    cluster_rows=True,
    cluster_cols=True,
    figure_size=(5, 5),
    cmap="RdPu",
    log=False,
    return_linkage=False,
    equal=True,
    title="",
    fun_type="heatmap",
    array_size=None,
    vmin=None,
    vmax=None,
):
    r"""Plot heatmap with hierarchically clustered rows and columns using `cell2location.plt.plot_heatmap.heatmap()`
    and `cell2location.plt.plot_heatmap.dotplot()`.

    :param df: pandas.DataFrame to be visualised using heatmap and dotplot
    :param cluster_rows: cluster rows or keep the same order as df?
    :param cluster_cols: cluster columns or keep the same order as df?
    :param figure_size: tuple specifying figure dimensions, passed to .heatmap
    :param cmap: pyplot colormap, passed to .heatmap
    :param log: boolean, color on logscale?
    :param return_linkage: return the plot or the plot + linkage for rows and columns? If true returns a dictionary with 'plot', 'row_linkage' and 'col_linkage' elements.
    :param equal: boolean, each tile should be square (equal aspect)
    :param title: clustermap title
    :param fun_type: 'heatmap' or 'dotplot'
    :param array_size: pandas.DataFrame to be visualised as dotplot dot size
       - must have the same dimensions as `df` (Default None)
    """

    if cluster_rows:
        # hierarchically cluster rows
        cor_f1 = np.corrcoef(df)
        row_linkage = hierarchy.linkage(distance.pdist(cor_f1), method="average")
        row_ord = hierarchy.leaves_list(row_linkage)
    else:
        row_linkage = None
        row_ord = np.arange(df.shape[0])

    if cluster_cols:
        # hierarchically cluster columns
        cor_f1 = np.corrcoef(df.T)
        col_linkage = hierarchy.linkage(distance.pdist(cor_f1), method="average")
        col_ord = hierarchy.leaves_list(col_linkage)
    else:
        col_linkage = None
        col_ord = np.arange(df.shape[1])

    df = df.iloc[row_ord, col_ord]

    if fun_type == "heatmap":
        # plot heatmap
        heatmap(
            df.values,
            ticks=False,
            log=log,
            figsize=figure_size,
            equal=equal,
            cmap=cmap,
            row_labels=df.index,
            col_labels=df.columns,
            cbar=True,
            title=title,
            vmin=vmin,
            vmax=vmax,
        )
    elif fun_type == "dotplot":
        # plot dotplot
        dotplot(
            df.values,
            array_size=array_size,
            ticks=False,
            log=log,
            figsize=figure_size,
            equal=equal,
            cmap=cmap,
            row_labels=df.index,
            col_labels=df.columns,
            cbar=True,
            title=title,
        )

    if return_linkage:
        return {"row_linkage": row_linkage, "col_linkage": col_linkage, "row_ord": row_ord, "col_ord": col_ord}





def plot_neighborhood_marker_dotplot_from_agg(
    abund_per_nh,
    included_niches=None,
    top_n_per_niche=5,
    show_excluded_rows=False,
    size_scale=250.0,
    size_power=1.6,
    size_floor=1.0,
    cmap="RdPu",
    marker_score_mode='prod',
    center_zero=False,          # NEW
    clip_percentile=None,       # NEW (e.g. 99)
    save_path=None,
    dpi=300
):
    """
    Plots a cell-type marker dotplot from an aggregated table of counts per neighborhood.

    Parameters
    ----------
    abund_per_nh : pd.DataFrame
        Aggregated matrix of counts (Rows = Neighborhoods, Columns = Cell types).
    included_niches : list
        List of niche IDs (as strings or ints) whose top-N markers contribute to the plot.
    top_n_per_niche : int, default 5
        Number of top markers to pick from each included niche.
    show_excluded_rows : bool, default False
        If True, plots all niches but dims the excluded ones. If False, plots only included.
    size_scale, size_power, size_floor : float
        Parameters controlling dot sizing based on cell-type specificity.
    cmap : str, default "RdPu"
        Matplotlib colormap for the marker score.
    save_path : str, optional
        If provided, the figure is saved to this path.
    dpi : int, default 300
        Resolution for the saved figure.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes objects
    metrics : dict
        A dictionary containing the calculated 'composition', 'specificity', and 'marker_score'.
    """

    # Ensure indices match formatting (cast to string)
    abund_per_nh = abund_per_nh.copy()

    if included_niches is None:
        included_niches = abund_per_nh.index

    # 1. Core matrices calculations
    # Composition = fraction of a neighborhood made up of cell type X
    composition  = abund_per_nh.div(abund_per_nh.sum(axis=1), axis=0).fillna(0)

    # Specificity = fraction of total cell type X found in a given neighborhood
    specificity  = abund_per_nh.div(abund_per_nh.sum(axis=0), axis=1).fillna(0)

    # Global frequency of cell type X
    ct_global    = abund_per_nh.sum(axis=0) / abund_per_nh.values.sum()

    eps          = 1e-6
    log2fc       = np.log2((composition + eps).div(ct_global + eps, axis=1))

    if marker_score_mode == 'prod':
        marker_score = composition * log2fc
    elif marker_score_mode == 'log2fc':
        marker_score = log2fc
    elif marker_score_mode == 'props':
        marker_score = composition

    if not center_zero:
        marker_score = composition * log2fc.clip(lower=0)

    # 2. Validate included niches exist
    missing = [n for n in included_niches if n not in marker_score.index]
    if missing:
        raise ValueError(f"included_niches contains values not in data: {missing}")

    # 3. Pick top N markers per included niche, take the union
    top_celltypes = set()
    for nh_id in included_niches:
        top_n = marker_score.loc[nh_id].sort_values(ascending=False).head(top_n_per_niche)
        top_celltypes.update(top_n.index)

    print(f"Cell types in plot: {len(top_celltypes)} "
          f"(top {top_n_per_niche} from {len(included_niches)} niches)")

    # 4. Niche row order — cluster on composition
    all_niches = sorted(composition.index, key=lambda x: int(x) if x.isdigit() else x)
    plot_niches = all_niches if show_excluded_rows else included_niches

    if len(all_niches) > 2:
        Z = linkage(np.sqrt(composition.loc[all_niches].values),
                    method="average", metric="euclidean")
        full_order = [all_niches[i] for i in leaves_list(Z)]
    else:
        full_order = all_niches

    nhood_order =[n for n in full_order if n in plot_niches]

    # 5. Cell type order — by best-niche assignment, restricted to included niches
    included_idx_map = {nh: i for i, nh in enumerate(included_niches)}
    ms_included = marker_score.loc[included_niches, list(top_celltypes)]
    best_niche_per_ct = ms_included.idxmax(axis=0)
    best_score_per_ct = ms_included.max(axis=0)

    ct_table = pd.DataFrame({
        "best_niche_idx": best_niche_per_ct.map(included_idx_map),
        "best_score":     best_score_per_ct,
    })
    ct_order = ct_table.sort_values(
        by=["best_niche_idx", "best_score"],
        ascending=[True, False],
    ).index.tolist()

    # 6. Extract matrices for plotting
    color_M = marker_score.loc[nhood_order, ct_order].values
    size_M  = specificity.loc[nhood_order, ct_order].values

    # Dot sizing math
    spec_max  = specificity.loc[nhood_order, ct_order].values.max()
    dot_areas = (size_M / max(spec_max, 1e-9)) ** size_power * size_scale
    dot_areas = np.maximum(dot_areas, size_floor)

    # 7. Setup plot layout
    n_niches = len(nhood_order)
    n_cts    = len(ct_order)
    fig_w    = max(10.0, 0.30 * n_cts + 4.0)
    fig_h    = max(3.0,  0.32 * n_niches + 2.5)

    fig, ax  = plt.subplots(figsize=(fig_w, fig_h))
    xx, yy = np.meshgrid(np.arange(n_cts), np.arange(n_niches))

    # Faint grids
    for x in range(n_cts):
        ax.axvline(x, color="#ededed", linewidth=0.4, zorder=1)
    for y in range(n_niches):
        ax.axhline(y, color="#ededed", linewidth=0.4, zorder=1)

    # Dim rows that are excluded (visual cue that they didn't contribute markers)
    if show_excluded_rows:
        for i, nh_id in enumerate(nhood_order):
            if nh_id not in included_niches:
                ax.axhspan(i - 0.5, i + 0.5, color="#fafafa", zorder=0)


    # Color scaling
    color_vals = color_M.ravel()

    if clip_percentile is not None:
        if center_zero:
            vmax = np.percentile(np.abs(color_vals), clip_percentile)
            vmin = -vmax
        else:
            vmax = np.percentile(color_vals, clip_percentile)
            vmin = 0
    else:
        if center_zero:
            vmax = np.abs(color_vals).max()
            vmin = -vmax
        else:
            vmax = color_vals.max()
            vmin = 0

    # Plot dots
    sc = ax.scatter(
        xx.ravel(), yy.ravel(),
        s=dot_areas.ravel(),
        c=color_M.ravel(),
        cmap=cmap, vmin=vmin, vmax=vmax,
        edgecolors="#333333", linewidths=0.3,
        zorder=3,
    )

    # Axes Formatting
    ax.set_xticks(range(n_cts))
    ax.set_xticklabels(ct_order, fontsize=9, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlabel("Cell type", fontsize=11, labelpad=8)

    ax.set_yticks(range(n_niches))
    yticklabels =[
        f"{nh_id}" + ("" if nh_id in included_niches else " (excl.)")
        for nh_id in nhood_order
    ]
    ax.set_yticklabels(yticklabels, fontsize=10, fontweight="semibold")
    ax.set_ylabel("Neighborhood", fontsize=10, labelpad=8)

    ax.set_xlim(-0.5, n_cts - 0.5)
    ax.set_ylim(-0.5, n_niches - 0.5)
    ax.invert_yaxis()
    ax.tick_params(length=0)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)

    # title = (
    #     f"Cell-type markers per neighborhood "
    #     f"(top {top_n_per_niche} from N{', N'.join(included_niches)})"
    # )
    # ax.set_title(title, fontsize=12, fontweight="semibold", loc="left", pad=8)

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, fraction=0.012, pad=0.01, shrink=0.85, anchor=(0.0, 0.5))
    cbar.set_label("Marker score\n(composition × log₂ FC⁺)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Specificity size legend
    candidates =[0.05, 0.1, 0.2, 0.4, 0.7]
    legend_vals = [v for v in candidates if v <= spec_max * 1.05]
    if not legend_vals:
        legend_vals =[round(spec_max * f, 2) for f in [0.25, 0.5, 1.0]]

    legend_handles = [
        plt.scatter([],[],
                    s=(v / max(spec_max, 1e-9)) ** size_power * size_scale,
                    c="#7a7a7a", edgecolors="#333333", linewidths=0.3,
                    label=f"{v:.2f}")
        for v in legend_vals
    ]

    sleg = ax.legend(
        handles=legend_handles,
        loc="upper left", bbox_to_anchor=(1.2, .8), # Keep it relatively close
        title="Specificity\n(fraction of cell type)",
        title_fontsize=8, fontsize=8,
        frameon=False,
        # labelspacing=1.5,   # <--- FIX 1: Adds vertical space so big dots don't overlap each other
        # borderpad=10.5,      # <--- FIX 2: Expands the invisible bounding box so it wraps the big dots
        # borderaxespad=10.0,
        handletextpad=1.5,  # <--- FIX 3: Adds space between the big dots and their text labels
    )
    sleg.get_title().set_multialignment("left")
    # ax.add_artist(sleg)
    sleg.set_clip_on(False)

    # Keep a normal rect size so the main plot isn't squished!
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    if save_path:
        # bbox_extra_artists combined with the new 'borderpad' ensures zero clipping
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight", bbox_extra_artists=(sleg,))
        print(f"Saved: {save_path}")

    metrics = {
        "composition": composition,
        "specificity": specificity,
        "marker_score": marker_score
    }

    return fig, ax, metrics
