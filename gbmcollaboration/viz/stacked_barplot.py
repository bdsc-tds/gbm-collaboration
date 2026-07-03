import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

def clustered_stacked_barplot(
    df,
    figsize=(12, 6),

    # clustering
    cluster_rows=True,
    cluster_cols=False,
    method="ward",
    metric="euclidean",

    # stacked bar colors
    cmap="tab20",
    custom_colors=None,

    # optional annotation (categorical)
    col_annotation=None,
    col_annotation_cmap="Set2",
    col_annotation_colors=None,
    col_annotation_label="",

    # annotation layout
    annotation_height=0.06,
    annotation_pad=0.02,

    # plotting
    normalize=False,
    linewidth=0.5,
    edgecolor="white",
    bar_width=0.9,

    # labels
    xlabel="",
    ylabel="Proportion",
    title="",

    # legend
    legend=True,
    legend_kwargs=None,

    # ticks
    xtick_rotation=90,

    # axes
    ax=None,
):
    """
    Clustered stacked barplot with optional categorical x-axis annotation.

    Parameters
    ----------
    df : pd.DataFrame
        Rows = samples/groups (X-axis)
        Columns = stacked categories (Bars)
        Values = proportions or abundances

    col_annotation : pd.Series or dict-like
        Categorical annotation for the items on the x-axis.
        Must align with df.index.

    col_annotation_colors : dict, optional
        Mapping category -> color

    Returns
    -------
    fig, ax, ordered_df
    """

    plot_df = df.copy()

    # ---------------------------------------------------
    # normalize rows if requested
    # ---------------------------------------------------
    if normalize:
        plot_df = plot_df.div(plot_df.sum(axis=1), axis=0)

    # ---------------------------------------------------
    # cluster rows
    # ---------------------------------------------------
    if cluster_rows and len(plot_df) > 1:
        row_linkage = linkage(
            pdist(plot_df.values, metric=metric),
            method=method
        )
        row_order = leaves_list(row_linkage)
        plot_df = plot_df.iloc[row_order]

    # ---------------------------------------------------
    # cluster columns
    # ---------------------------------------------------
    if cluster_cols and len(plot_df.columns) > 1:
        col_linkage = linkage(
            pdist(plot_df.T.values, metric=metric),
            method=method
        )
        col_order = leaves_list(col_linkage)
        plot_df = plot_df.iloc[:, col_order]

    # ---------------------------------------------------
    # reorder annotation
    # ---------------------------------------------------
    if col_annotation is not None:
        if isinstance(col_annotation, pd.Series):
            annot = col_annotation.loc[plot_df.index]
        else:
            annot = pd.Series(col_annotation, index=df.index)
            annot = annot.loc[plot_df.index]

    # ---------------------------------------------------
    # stacked colors
    # ---------------------------------------------------
    n_cols = plot_df.shape[1]

    if custom_colors is not None:
        if isinstance(custom_colors, dict):
            colors = [custom_colors[col] for col in plot_df.columns]
        else:
            colors = custom_colors
    else:
        cmap_obj = plt.get_cmap(cmap)
        colors = [cmap_obj(i / max(n_cols - 1, 1)) for i in range(n_cols)]

    # ---------------------------------------------------
    # figure layout
    # ---------------------------------------------------
    if ax is None:
        fig = plt.figure(figsize=figsize)

        if col_annotation is not None:
            # FIX 3: Flipped the height_ratios so the annotation bar is at the bottom
            gs = fig.add_gridspec(
                2,
                1,
                height_ratios=[1, annotation_height],
                hspace=annotation_pad
            )
            ax = fig.add_subplot(gs[0, 0])
            annot_ax = fig.add_subplot(gs[1, 0])
        else:
            ax = fig.add_subplot(111)
            annot_ax = None
    else:
        fig = ax.figure
        annot_ax = None

    # ---------------------------------------------------
    # x positions
    # ---------------------------------------------------
    x = np.arange(len(plot_df.index))

    # ---------------------------------------------------
    # stacked bars
    # ---------------------------------------------------
    bottom = np.zeros(len(plot_df))

    for col, color in zip(plot_df.columns, colors):
        ax.bar(
            x,
            plot_df[col],
            bottom=bottom,
            label=col,
            color=color,
            linewidth=linewidth,
            edgecolor=edgecolor,
            width=bar_width
        )
        bottom += plot_df[col].values

    # ---------------------------------------------------
    # categorical annotation strip
    # ---------------------------------------------------
    if col_annotation is not None and annot_ax is not None:
        categories = pd.unique(annot)

        # category colors
        if col_annotation_colors is not None:
            color_dict = col_annotation_colors
        else:
            cmap_obj = plt.get_cmap(col_annotation_cmap)
            color_dict = {
                cat: cmap_obj(i / max(len(categories) - 1, 1))
                for i, cat in enumerate(categories)
            }

        # annotation blocks
        for i, cat in enumerate(annot):
            annot_ax.add_patch(
                Rectangle(
                    (i - 0.5, 0), # FIX 2: Shifted by -0.5 to align seamlessly with the bars above it
                    1,
                    1,
                    color=color_dict[cat]
                )
            )

        # FIX 1: set_xlim matching the index (rows) size, not the stack sizes (columns)
        annot_ax.set_xlim(-0.5, len(plot_df.index) - 0.5)
        annot_ax.set_ylim(0, 1)

        # Drop frames on annotation axis
        for spine in annot_ax.spines.values():
            spine.set_visible(False)

        annot_ax.set_yticks([])

        # Move the X-ticks to the bottom (the annotation axis)
        annot_ax.set_xticks([])
        # annot_ax.set_xticks(x)
        # annot_ax.set_xticklabels(
        #     plot_df.index.astype(str),
        #     rotation=xtick_rotation
        # )
        annot_ax.set_xlabel(xlabel)

        # Hide the X-ticks on the main plotting axis
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_xlabel("")

        # annotation legend
        handles = [
            Rectangle((0, 0), 1, 1, color=color_dict[cat])
            for cat in categories
        ]

        annot_ax.legend(
            handles,
            categories,
            title=col_annotation_label,
            loc="center left",
            bbox_to_anchor=(1.02, 1),
            frameon=False,
            **legend_kwargs

        )
    else:
        # Standard x-axis labels if no annotation bar exists
        ax.set_xticks(x)
        ax.set_xticklabels(
            plot_df.index.astype(str),
            rotation=xtick_rotation
        )
        ax.set_xlabel(xlabel)

    # ---------------------------------------------------
    # formatting
    # ---------------------------------------------------
    ax.tick_params(axis='y', labelsize=20)
    ax.set_ylabel(ylabel,fontsize=20)
    ax.set_title(title)

    # remove excess whitespace
    ax.set_xlim(-0.5, len(plot_df.index) - 0.5)
    ax.set_ylim(0,1)
    ax.margins(x=0)

    # cleaner aesthetics
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


    # ---------------------------------------------------
    # legend
    # ---------------------------------------------------
    if legend:
        if legend_kwargs is None:
            legend_kwargs = {}

        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            **legend_kwargs
        )

    plt.tight_layout()

    return fig, ax, plot_df
