"""Supp Fig 3e-f: single-cell CNV clones and WES CNV clusters.

Produces Supp Fig 3e-f (`WES_scCNV_integration_dot_facets`): (e) fraction of each
single-cell inferCNV clone per WES CNV cluster, coloured by the log2 enrichment over the
mean of the other WES clusters; (f) per-patient fraction of the most frequent clones in WES
cluster 2 vs the others. Also the Mann-Whitney test of clone 4 (the chr19 clone), cluster 2
vs the others ("a chr19-enriched clone disproportionately represented in cluster 2 tumors").

Origin: analysis/infercnv/analyses_wes_sc_infercnv.py. Inputs (frozen): the WES cluster
labels (results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv, from
fig3b_suppfig3ac_infercnv_wes_analyses.py) and the obs of the malignant-cell inferCNV object
(pipelines/infercnv_sc_analyses.py). Only the panels of the paper are kept.

Comparison with the published panel: the malignant-cell object has one sample without a WES
profile (ER_G_743a, 42 cells, wes_cluster NaN). It is a column of none of the dot-plot rows,
but `wes_cluster != 2` counts it in "WES clusters 0, 1, 3" in panel f and in the clone 4 test
(the printed n of the other group includes it), as in the published figure.
"""

# %% -- Import packages
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

from gbmcollaboration.paper import output

OUT = output.setup("suppfig3ef_analyses_wes_sc_infercnv")
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results")
SC_ADATA_PATH = RESULTS_DIR / "scRNAseq/infercnv/infercnv_adata_malignant_only.h5ad"


# %% -- Functions
def frac_to_area(f):
    if fmax == fmin:
        return (size_min + size_max) / 2
    return size_min + (f - fmin) * (size_max - size_min) / (fmax - fmin)


# %% -- Load cluster label WES
df_wes_leiden_cluster = pd.read_csv(
    RESULTS_DIR / "WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
)
df_wes_leiden_cluster.rename(columns={"cnv_leiden": "wes_cluster"}, inplace=True)

# merge with the adata (obs only)
adata_backed = ad.read_h5ad(SC_ADATA_PATH, backed="r")
adata_tumor = ad.AnnData(obs=adata_backed.obs.copy())
adata_backed.file.close()
adata_tumor.obs = adata_tumor.obs.join(
    df_wes_leiden_cluster.set_index("orig.ident")[["wes_cluster"]],
    on="orig.ident",
    how="left"
)

# check
print(f"{adata_tumor.n_obs:,} malignant cells from {adata_tumor.obs['orig.ident'].nunique()} samples")
print(adata_tumor.obs["wes_cluster"].value_counts(dropna=False))
no_wes = adata_tumor.obs.loc[adata_tumor.obs["wes_cluster"].isna(), "orig.ident"].value_counts()
print("Samples without a WES cluster (counted in 'WES clusters 0, 1, 3' below):", no_wes.to_dict())

# -----------------------------
# 1) Dotplot: CNV clone composition across WES clusters
# -----------------------------

# %% Cross-tab (rows: WES clusters, cols: CNV clones) -> fractions per WES cluster
ct_wes = pd.crosstab(
    adata_tumor.obs["wes_cluster"],
    adata_tumor.obs["cnv_leiden"],
    normalize="index"
).sort_index()

# Mean fraction of each clone in "other" WES clusters (leave-one-out mean)
# shape: same as ct_wes (rows x cols)
mean_others = ct_wes.apply(lambda col: (col.sum() - col) / (len(col) - 1), axis=0)

# Relative enrichment and log2 enrichment (vectorized)
eps = 1e-6
rel_enrich = ct_wes / (mean_others + eps)
log2_enrich = np.log2(rel_enrich)

# Long format for plotting
df_dot = (
    ct_wes.stack()
    .rename("fraction")
    .to_frame()
    .join(log2_enrich.stack().rename("log2_enrichment"))
    .reset_index()
    .rename(columns={"level_0": "wes_cluster", "level_1": "cnv_leiden"})
)

# -----------------------------
# 2) Patient-level test: is clone 4 enriched in WES cluster 2?
# -----------------------------

# %% Ensure cnv_leiden is treated consistently (int if possible)
cnv_series = adata_tumor.obs["cnv_leiden"]
cnv_is_int = pd.api.types.is_integer_dtype(cnv_series) or pd.api.types.is_numeric_dtype(cnv_series)

clone4_mask = (cnv_series == 4) if cnv_is_int else (cnv_series.astype(str) == "4")

# Fraction of clone 4 per patient (patient is the unit)
frac_clone4 = (
    adata_tumor.obs.assign(is_clone4=clone4_mask)
    .groupby("orig.ident", observed=True)["is_clone4"]
    .mean()
    .rename("frac_clone4")
)

# Unique WES cluster label per patient
patient_wes = (
    adata_tumor.obs[["orig.ident", "wes_cluster"]]
    .drop_duplicates(subset=["orig.ident"])
    .set_index("orig.ident")
)

df_test = frac_clone4.to_frame().join(patient_wes, how="inner")

# Split groups
g2 = df_test.loc[df_test["wes_cluster"] == 2, "frac_clone4"]
gO = df_test.loc[df_test["wes_cluster"] != 2, "frac_clone4"]

# Mann–Whitney U test
u_stat, pval = mannwhitneyu(g2, gO, alternative="two-sided")

# Effect size: median difference (cluster2 - others)
median_diff = g2.median() - gO.median()

print(f"N patients (WES cluster 2): {g2.shape[0]}")
print(f"N patients (Other clusters): {gO.shape[0]}")
print(f"Mann–Whitney U = {u_stat:.1f}, p = {pval:.3e}")
print(f"Median difference (cluster2 - others) = {median_diff:.3f}")

# -----------------------------
# 3) Final figure with the stripplot (Supp Fig 3e-f)
# -----------------------------
# %%
# ============================================================
# Figure: (A) Dotplot WES clusters x single-cell CNV clones
#         (B) Faceted boxen+strip for selected clones (patient-level),
#             with dashed median lines per group
# ============================================================

# -----------------------------
# Params
# -----------------------------
palette_groups = {"WES cluster 2": "#E3B23C", "WES clusters 0, 1, 3": "#3B6FB6"}
vmax = 3
dot_sizes = (40, 300)
clones_keep = ["0", "1", "2", "4", "5", "7"]  # selected clones for panel B (strings are safest)

# -----------------------------
# Build patient-level fractions (Option A):
# fraction = (#cells in clone) / (total malignant cells in the patient)
# -----------------------------
obs = adata_tumor.obs.copy()
obs["cnv_leiden"] = obs["cnv_leiden"].astype(str)
obs["WES_group"] = np.where(obs["wes_cluster"] == 2, "WES cluster 2", "WES clusters 0, 1, 3")

# counts per patient x clone
counts = (
    obs.groupby(["orig.ident", "WES_group", "cnv_leiden"], observed=True)
       .size()
       .rename("n_cells")
       .reset_index()
)

# total malignant cells per patient (all clones)
totals = (
    obs.groupby(["orig.ident", "WES_group"], observed=True)
       .size()
       .rename("total_cells")
       .reset_index()
)

df_long = (
    counts.merge(totals, on=["orig.ident", "WES_group"], how="left")
          .assign(fraction=lambda x: x["n_cells"] / x["total_cells"])
)

# keep only selected clones for panel B
df_long_sel = df_long[df_long["cnv_leiden"].isin(clones_keep)].copy()
df_long_sel["cnv_leiden"] = pd.Categorical(df_long_sel["cnv_leiden"],
                                           categories=clones_keep, ordered=True)

# ============================================================
# Build the full figure canvas
# ============================================================
fig = plt.figure(figsize=(16, 6.2), dpi=300)
outer = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.8)

# ============================================================
# Panel A — Dot plot (left)
# ============================================================
ax1 = fig.add_subplot(outer[0])

norm = mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

# Use ONE cmap for both scatter + colorbar
cmap = sns.color_palette("icefire", as_cmap=True)  # or "icefire"

sns.scatterplot(
    data=df_dot,
    x="cnv_leiden",
    y="wes_cluster",
    size="fraction",
    hue="log2_enrichment",
    palette=cmap,          # <-- changed
    hue_norm=norm,
    sizes=dot_sizes,
    edgecolor="black",
    linewidth=0.3,
    legend=False,
    ax=ax1
)

ax1.set_yticks([0, 1, 2, 3])
ax1.set_yticklabels([0, 1, 2, 3])
ax1.set_xlabel("Single-cell–derived CNV clones")
ax1.set_ylabel("WES-derived CNV clusters")
ax1.set_title("Global CNV clonal composition", fontsize=12, pad=10)
sns.despine(ax=ax1)

# Colorbar (same cmap + same norm)
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax1, fraction=0.046, pad=0.02)
cbar.set_label("log2 enrichment")

# Size legend for fraction (panel A)
size_vals = [0.05, 0.15, 0.25]
size_min, size_max = dot_sizes
fmin, fmax = df_dot["fraction"].min(), df_dot["fraction"].max()

handles = [
    plt.Line2D([], [], linestyle="none", marker="o",
               markersize=np.sqrt(frac_to_area(v)),
               markerfacecolor="white", markeredgecolor="black",
               markeredgewidth=0.6)
    for v in size_vals
]
labels = [f"{v:.2f}" for v in size_vals]

ax1.legend(
    handles, labels,
    title="fraction",
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(1.18, 1),
    borderaxespad=0.0
)

# ============================================================
# Panel B — FacetGrid (right) inside an inset subgridspec
# ============================================================
# Create an axis placeholder and then delete it (FacetGrid makes its own axes)
ax_placeholder = fig.add_subplot(outer[1])
ax_placeholder.remove()

# Subgrid for 5 rows, 1 column in the right panel
right = outer[1].subgridspec(len(clones_keep), 1, hspace=0.35)

axes = [fig.add_subplot(right[i, 0]) for i in range(len(clones_keep))]

for i, clone in enumerate(clones_keep):
    ax = axes[i]
    sub = df_long_sel[df_long_sel["cnv_leiden"].astype(str) == str(clone)].copy()

    sns.boxenplot(
        data=sub,
        x="fraction",
        y="WES_group",
        orient="h",
        palette=palette_groups,
        linewidth=1.0,
        ax=ax
    )

    # strip
    sns.stripplot(
        data=sub,
        x="fraction",
        y="WES_group",
        orient="h",
        hue="WES_group",
        palette=palette_groups,
        dodge=False,
        jitter=0.20,
        size=2.0,
        alpha=0.35,
        linewidth=0,
        ax=ax
    )

    # dashed medians (one per group, per facet)
    for grp, subg in sub.groupby("WES_group"):
        med = subg["fraction"].median()
        ax.axvline(
            med,
            color=palette_groups.get(grp, "black"),
            linestyle="--",
            linewidth=1.2,
            alpha=0.9,
            zorder=6
        )

    # clean duplicated legends
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    ax.set_xlim(0, 1.0)
    ax.set_ylabel("")
    ax.set_title(f"Single-cell-derived CNV clone {clone}", loc="left", fontsize=10, pad=2)
    sns.despine(ax=ax)

    # only bottom axis gets x label / ticks
    if i < len(clones_keep) - 1:
        ax.set_xlabel("")
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=0)
    else:
        ax.set_xlabel("Fraction per patient")

# Title for the whole right panel
axes[0].text(
    0.0, 1.35,
    "Patient-level distribution of most frequent clones",
    transform=axes[0].transAxes,
    fontsize=12,
    ha="left",
    va="bottom"
)

# Global margins
fig.subplots_adjust(left=0.06, right=0.80, top=0.90, bottom=0.12, wspace=0.25)

output.save(OUT, "WES_scCNV_integration_dot_facets", fig, df_dot, panel="Supp Fig 3e-f")
plt.show()
plt.close(fig)
print(log2_enrich.round(2).to_string())

# %%
print(output.versions(["numpy", "pandas", "scipy", "seaborn", "matplotlib"]).to_string(index=False))
