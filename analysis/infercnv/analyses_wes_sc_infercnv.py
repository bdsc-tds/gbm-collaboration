"""Script to perform additional analyses on WES and single cell CNV clustering."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scanpy as sc
import scipy
import scipy.sparse as sp
import infercnvpy as cnv
from pathlib import Path
import anndata as ad
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# %% -- Load cluster label WES
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results")

df_wes_leiden_cluster = pd.read_csv(
    RESULTS_DIR / "WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
)
df_wes_leiden_cluster.rename(columns={"cnv_leiden": "wes_cluster"}, inplace=True)

# merge with the adata
SC_ADATA_PATH = RESULTS_DIR / "scRNAseq/infercnv/infercnv_adata_malignant_only.h5ad"
adata_tumor = ad.read_h5ad(SC_ADATA_PATH)
adata_tumor.obs = adata_tumor.obs.merge(
    df_wes_leiden_cluster[["orig.ident", "wes_cluster"]],
    on="orig.ident",
    how="left"
)

# check
adata_tumor.obs["wes_cluster"].value_counts()

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

# Dot plot
fig, ax = plt.subplots(figsize=(10, 3.5), dpi=300)
sns.scatterplot(
    data=df_dot,
    x="cnv_leiden",
    y="wes_cluster",
    size="fraction",
    hue="log2_enrichment",
    palette="RdBu_r",
    sizes=(20, 300),
    edgecolor="black",
    linewidth=0.3,
    ax=ax
)
ax.set_xlabel("Single-cell-derived CNV clones")
ax.set_ylabel("WES-derived CNV cluster")

# Clean legend: keep seaborn's combined legend but move it outside
ax.legend(
    title="",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()
plt.show()

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
    .groupby("orig.ident")["is_clone4"]
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
# 3) Final figure with the stripplot
# -----------------------------
# %%
# ---------- choisir une représentation robuste ----------
adata_tumor.obs["cnv_leiden_plot"] = adata_tumor.obs["cnv_leiden"].astype(str)
clones_zoom = ["5", "7", "1", "2", "4"]

# ---------- construire df_zoom_frac (patient-level fractions) ----------
df_zoom = (
    adata_tumor.obs
    .loc[adata_tumor.obs["cnv_leiden_plot"].isin(clones_zoom)]
    .assign(
        WES_group=lambda x: np.where(x["wes_cluster"] == 2, "WES cluster 2", "Other")
    )
)

# Si df_zoom est vide, on stoppe avec un message explicite
if df_zoom.shape[0] == 0:
    raise ValueError("df_zoom est vide. Vérifie le type de cnv_leiden et clones_zoom.")

counts = (
    df_zoom
    .groupby(["orig.ident", "WES_group", "cnv_leiden_plot"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)

totals = (
    counts
    .groupby("orig.ident", observed=True)["n_cells"]
    .sum()
    .rename("total_cells")
    .reset_index()
)

df_zoom_frac = (
    counts.merge(totals, on="orig.ident")
    .assign(fraction=lambda x: x["n_cells"] / x["total_cells"])
)

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

# -----------------------------
# Helper: add dashed medians in FacetGrid panels
# -----------------------------
def add_group_medians(data, **kws):
    ax = plt.gca()
    for grp, sub in data.groupby("WES_group"):
        med = sub["fraction"].median()
        ax.axvline(
            med,
            color=palette_groups.get(grp, "black"),
            linestyle="--",
            linewidth=1.2,
            alpha=0.9,
            zorder=6
        )

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

def frac_to_area(f):
    if fmax == fmin:
        return (size_min + size_max) / 2
    return size_min + (f - fmin) * (size_max - size_min) / (fmax - fmin)

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

plt.show()

# Optional save:
fig.savefig(RESULTS_DIR / "WES/infercnv/WES_scCNV_integration_dot_facets.pdf", bbox_inches="tight")

# %% -- Heatmap showing the log2FC of malignant transcriptional states
# between cluster 4 and the other clusters

clone_of_interest = "4"

# distributions des états dans clone 4 et dans le reste
s4 = adata_tumor.obs.loc[adata_tumor.obs["cnv_leiden"] == clone_of_interest, "celltype_level4_scanvi"]
sR = adata_tumor.obs.loc[adata_tumor.obs["cnv_leiden"] != clone_of_interest, "celltype_level4_scanvi"]

p4 = (s4.value_counts(normalize=True)).rename("p_clone4")
pR = (sR.value_counts(normalize=True)).rename("p_rest")

df = pd.concat([p4, pR], axis=1).fillna(0)

# pseudo-count pour éviter division par 0
eps = 1e-4
df["log2FC_clone4_vs_rest"] = np.log2((df["p_clone4"] + eps) / (df["p_rest"] + eps))

# tri par déplétion (valeurs les plus négatives)
df_sorted = df.sort_values("log2FC_clone4_vs_rest")
df_sorted.head(10)

# garder tous les états, mais ordonner par log2FC
heat = df_sorted[["log2FC_clone4_vs_rest"]].T

plt.figure(figsize=(10, 1.2), dpi=300)
sns.heatmap(
    heat,
    cmap="coolwarm",
    center=0,
    cbar_kws={"label": "log2FC (clone 4 vs rest)"}
)
plt.yticks([0.5], ["Clone 4 vs others"], rotation=0)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()

# %%
# ============================================================
# Differential cell-type composition: WES cluster 2 vs others
# Patient-level stats + publication-ready figure
# ============================================================

RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/scRNAseq/infercnv")
SC_ADATA_PATH = RESULTS_DIR / "infercnv_adata.h5ad"

WES_CLUSTER_PATH = Path(
    "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
)

# -----------------------------
# Load + basic filtering
# -----------------------------
adata = ad.read_h5ad(SC_ADATA_PATH)

celltype_col = "celltype_level4_scanvi"  # <-- use level4 since you plot level4
drop_level1 = ["Doublet", "Doublet_TNK_nos"]
drop_level4 = ["Other_LowCounts"]

# If you want to keep your original "celltype_level1_curated" filtering too:
if "celltype_level1_curated" in adata.obs.columns:
    adata = adata[~adata.obs["celltype_level1_curated"].isin(drop_level1)].copy()

if celltype_col in adata.obs.columns:
    adata = adata[~adata.obs[celltype_col].isin(drop_level4)].copy()
else:
    raise ValueError(f"{celltype_col} not found in adata.obs")

# -----------------------------
# Add WES clusters (patient-level label)
# -----------------------------
df_wes = pd.read_csv(WES_CLUSTER_PATH).rename(columns={"cnv_leiden": "wes_cluster"})
df_wes = df_wes[["orig.ident", "wes_cluster"]].drop_duplicates()

# map (safer than merge into .obs)
adata.obs = adata.obs.join(df_wes.set_index("orig.ident"), on="orig.ident")

# keep only patients with wes_cluster defined
adata = adata[adata.obs["wes_cluster"].notna()].copy()
adata.obs["wes_cluster"] = adata.obs["wes_cluster"].astype(int)

# define group of interest
cluster_of_interest = 2
adata.obs["WES_group"] = np.where(
    adata.obs["wes_cluster"].eq(cluster_of_interest),
    "WES cluster 2",
    "Other"
)

# -----------------------------
# Patient-level fractions per cell type
# (THIS avoids pseudo-replication)
# -----------------------------
# %% counts per patient x celltype
df_counts = (
    adata.obs
    .groupby(["orig.ident", "WES_group", celltype_col], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)

# total malignant/tumor cells per patient (within the filtered adata)
df_counts["total_cells_patient"] = df_counts.groupby(["orig.ident"], observed=True)["n_cells"].transform("sum")
df_counts["fraction"] = df_counts["n_cells"] / df_counts["total_cells_patient"]

# to ensure celltypes absent in a patient appear as 0, build a full grid
patients = (
    adata.obs[["orig.ident", "WES_group"]]
    .drop_duplicates()
    .set_index("orig.ident")
)

celltypes = pd.Index(sorted(adata.obs[celltype_col].dropna().unique()), name=celltype_col)

full_index = pd.MultiIndex.from_product(
    [patients.index, ["WES cluster 2", "Other"], celltypes],
    names=["orig.ident", "WES_group", celltype_col]
)

df_full = (
    df_counts
    .set_index(["orig.ident", "WES_group", celltype_col])[["fraction", "n_cells"]]
    .reindex(full_index)
    .fillna({"fraction": 0.0, "n_cells": 0})
    .reset_index()
)

# However the product above duplicates patients in the wrong group. Fix by merging true group:
df_full = df_full.drop(columns=["WES_group"]).merge(
    patients.reset_index(),
    on="orig.ident",
    how="left"
)

# Recompute fraction is already correct (0 for missing celltypes)
# Also keep per-patient total cells for filtering if needed:
df_full = df_full.merge(
    df_counts.groupby("orig.ident", observed=True)["total_cells_patient"].first().rename("total_cells_patient"),
    on="orig.ident",
    how="left"
)

# -----------------------------
# Filter cell types (avoid very rare types)
# -----------------------------
min_cells_overall = 50
celltype_totals = adata.obs[celltype_col].value_counts()
keep_celltypes = celltype_totals[celltype_totals >= min_cells_overall].index

df_full = df_full[df_full[celltype_col].isin(keep_celltypes)].copy()

# optional: keep only celltypes seen in >= N patients in at least one group
min_patients = 5
pt_counts = (
    df_full.assign(present=df_full["fraction"] > 0)
    .groupby(["WES_group", celltype_col], observed=True)["present"]
    .sum()
    .unstack("WES_group")
    .fillna(0)
)
keep2 = pt_counts[(pt_counts["WES cluster 2"] >= min_patients) | (pt_counts["Other"] >= min_patients)].index
df_full = df_full[df_full[celltype_col].isin(keep2)].copy()

# -----------------------------
# Stats per cell type: Mann–Whitney + FDR
# Effect size = median difference + log2FC of medians (with epsilon)
# -----------------------------
eps = 1e-6
rows = []

for ct in sorted(df_full[celltype_col].unique()):
    sub = df_full[df_full[celltype_col] == ct]
    g2 = sub.loc[sub["WES_group"] == "WES cluster 2", "fraction"].values
    go = sub.loc[sub["WES_group"] == "Other", "fraction"].values

    # require some minimal sample size per group
    if (len(g2) < 5) or (len(go) < 5):
        continue

    stat, p = mannwhitneyu(g2, go, alternative="two-sided")

    med2 = np.median(g2)
    medO = np.median(go)

    rows.append({
        "celltype": ct,
        "median_cluster2": med2,
        "median_other": medO,
        "median_diff": med2 - medO,
        "log2FC_median": np.log2((med2 + eps) / (medO + eps)),
        "pval": p
    })

df_stats = pd.DataFrame(rows)
df_stats["padj"] = multipletests(df_stats["pval"], method="fdr_bh")[1]

# sort by effect size (log2FC of medians)
df_stats = df_stats.sort_values("log2FC_median")

# -----------------------------
# Figure: barplot of log2FC + FDR marks
# -----------------------------
plt.figure(figsize=(16, 3.2), dpi=300)

ax = sns.barplot(
    data=df_stats,
    x="celltype",
    y="log2FC_median",
    palette=sns.color_palette("vlag", as_cmap=False),
    edgecolor="none"
)

ax.axhline(0, color="black", linewidth=0.8)
ax.set_ylabel("log2 fold-change of median fraction\n(WES cluster 2 vs others)")
ax.set_xlabel("")
ax.set_title("Differential cell-type composition in WES CNV cluster 2", pad=10)

# annotate significance (FDR)
for i, row in enumerate(df_stats.itertuples(index=False)):
    if row.padj < 0.05:
        y = row.log2FC_median
        ax.text(
            i,
            y + (0.03 if y >= 0 else -0.03),
            "*" if row.padj >= 1e-3 else "**" if row.padj >= 1e-6 else "***",
            ha="center",
            va="bottom" if y >= 0 else "top",
            fontsize=9
        )

plt.xticks(rotation=45, ha="right")
sns.despine(ax=ax)
plt.tight_layout()
plt.show()

# -----------------------------
# Optional: if still want the 1-row heatmap style
# now with the patient-level log2FC (median) and stats table
# -----------------------------
# heat = df_stats.set_index("celltype")[["log2FC_median"]].T
# fig, ax = plt.subplots(figsize=(16, 2.6), dpi=300)
# sns.heatmap(heat, ax=ax, cmap="vlag", center=0, cbar_kws={"label": "log2FC (median fraction)"},
#             linewidths=0.3, linecolor="white")
# ax.set_yticks([0.5])
# ax.set_yticklabels(["WES cluster 2 vs others"], rotation=0)
# ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
# plt.tight_layout()
# plt.show()

# %%
df_stats.head(3)

# %%
df_stats.tail(3)

## don't see significant differences in the cell type composition
# between cluster 2 and the others.
