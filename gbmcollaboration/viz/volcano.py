
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def volcano(
    df,
    logfc_col="logfoldchanges",
    pval_col="pvals_adj",
    gene_col=None,
    color_col=None,
    fc_thresh=2,
    pval_thresh=0.05,
    top_n_labels=20,
    figsize=(8, 6),
    colors=("grey", "blue", "red"),
    ax=None,
):
    """
    Volcano plot from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with DE results
    logfc_col : str
        Column name for log fold change
    pval_col : str
        Column name for adjusted p-value
    gene_col : str
        Column name for gene names; if None, use df.index
    fc_thresh : float
        Fold-change threshold for highlighting
    pval_thresh : float
        P-value threshold for significance
    top_n_labels : int
        Number of top up/down genes to label
    figsize : tuple
        Figure size
    colors : tuple
        Colors for (non-sig, down, up)
    """

    from adjustText import adjust_text

    # Ensure gene names
    if gene_col is None:
        df_plot = df.copy()
        df_plot["gene"] = df_plot.index
    else:
        df_plot = df.copy()
        df_plot["gene"] = df_plot[gene_col]

    # Calculate -log10(pval)
    df_plot["minus_log10_p"] = -np.log10(df_plot[pval_col] + 1e-300)
    is_sig = df_plot[pval_col] <= pval_thresh

    # Base plot
    if ax is None:
        f, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(
        data=df_plot[~is_sig],
        x=logfc_col,
        y="minus_log10_p",
        color=colors[0],
        s=10,
        alpha=0.7,
        label="Not significant",
        ax=ax,
    )

    # Highlight up- and down-regulated
    down = df_plot[(df_plot[logfc_col] <= -fc_thresh) & is_sig]
    up = df_plot[(df_plot[logfc_col] >= fc_thresh) & is_sig]

    if color_col is None:
        sns.scatterplot(data=down, x=logfc_col, y="minus_log10_p", color=colors[1], s=20, label="Down-regulated", ax=ax)
        sns.scatterplot(data=up, x=logfc_col, y="minus_log10_p", color=colors[2], s=20, label="Up-regulated", ax=ax)
    else:
        sns.scatterplot(data=df_plot[is_sig], x=logfc_col, y="minus_log10_p", hue=color_col, s=20, ax=ax)

    # Label top N genes by p-value
    texts = []
    if top_n_labels > 0:
        top_up = up.nsmallest(top_n_labels, pval_col)
        top_down = down.nsmallest(top_n_labels, pval_col)
        for _, r in pd.concat([top_up, top_down]).iterrows():
            texts.append(ax.text(r[logfc_col], r["minus_log10_p"], r["gene"], fontsize=8))
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5))

    # Threshold lines
    ax.axvline(-fc_thresh, color="grey", linestyle="--")
    ax.axvline(fc_thresh, color="grey", linestyle="--")
    ax.axhline(-np.log10(pval_thresh), color="grey", linestyle="--")

    ax.set_xlabel(logfc_col, fontsize=20)
    ax.set_ylabel("-log10(p-value)", fontsize=20)
    ax.legend()
    sns.despine(ax=ax)
