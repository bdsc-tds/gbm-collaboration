"""Pipeline: inferCNV on the single-nucleus data, and the malignant-cell object.

Outputs (frozen on EFS in results/scRNAseq/infercnv/; this file now writes to its output
folder, `output.setup`, not over them):
- Part 1: infercnv_adata.h5ad (all 1,350,338 nuclei of the 266 tumours, loader exclusions,
  inferCNV matrix obsm["X_cnv"], CNV Leiden clusters and CNV score; X emptied),
  genomic_windows.csv (the 1,346 bins of 200 genes, reused as WES bins by
  fig3b_0_cnv_wes_matrix.py), ith_scores.csv, tumor_clone_annotation.csv, cnv heatmaps.
- Part 2: infercnv_adata_malignant_only.h5ad (782,542 tumour-lineage nuclei: level-1
  Malignant or Other, doublets and Other_LowCounts removed, re-clustered on their CNV
  profiles with Leiden 0.4: the 15 single-cell CNV clones of Supp Fig 3e-f),
  tumor_clone_annotation_clean.csv, ith_scores_clean and the exploratory figures in between
  (iCNV heatmaps, KDE, Sankey, UMAP composites, clustermap, heterogeneity violin; not in the
  paper).
Read by: fig3b_suppfig3ac_infercnv_wes_analyses.py (sample IDs),
suppfig3ef_analyses_wes_sc_infercnv.py (obs of the malignant-cell object),
fig3b_0_cnv_wes_matrix.py (bins).

Origin: analysis/infercnv/infercnv_sc_analyses.py.

Inputs: AllDataCenters.load_singlecell (loader exclusions), the EFS annotation key
resources/SingleCell_Annotation_Key_v3.csv.

Runtime and memory: hours (inferCNV on 1.35 M nuclei, then ITH scores), well over 60 GB of
RAM (run under `flock /tmp/claude-1000/paper_heavy.lock`). Not re-run for the paper.

Comparison with the frozen outputs:
- The frozen malignant-cell object was written by an earlier version of Part 2 with other
  column names (`tumor_lineage`, `cnv_status`, `cnv_burden`, `tumor_cnv_group`). On the
  frozen obs, level-1 Malignant or Other after the filters gives exactly the 782,542 frozen
  cells, the subset made here (lineage_status == "tumor-lineage", stored as `tumor_lineage`).
- `ithgex` is computed on the Part 1 object whose X was emptied before saving, so it is not
  informative (not used in the paper).
- Leiden on 1.16 M / 0.78 M cells is not bit-for-bit reproducible across igraph versions;
  the figures use the frozen clone labels.
"""

# %% -- Import packages and prepare functions
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scanpy as sc
import scipy
import infercnvpy as cnv
from pathlib import Path
import anndata as ad
import plotly.graph_objects as go

from gbmcollaboration import AllDataCenters
from gbmcollaboration.utils.infercnv_utils import InferCNVAnalyzer
from gbmcollaboration.paper import output

############################################################
############### PART 1 ##############################
############################################################
# %% Params
RESULTS_DIR = output.setup("fig3b_0_infercnv_sc")
SIMPLIFIED_ANNOTATION_MAPPING = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/SingleCell_Annotation_Key_v3.csv")
REFERENCE_COLUMN = "celltype_level4_scanvi"
BATCH_KEY = "orig.ident"
REF_CELL_TYPES = [
    'Astrocyte_homeostatic',
    'B_nos',
    'Plasma_nos',
    'Endothelial_cap','Endothelial_tip',
    'Neutrophils_nos',
    'Microglia',
    'Monocyte_nos',
    'Pericyte',
    'Neuron','Neuron_exc_CDH24','Neuron_exc_nos','Neuron_neurogliaform',
    'Oligodendrocyte_nos',
    'NK_KLRC1','NK_nos',
    'T_CD4_nos','T_CD8_nos',
]
THRESHOLD_N_GENES=400
THRESHOLD_N_COUNTS=10000

# %% Create results dir
RESULTS_DIR.mkdir(parents=True,exist_ok=True)

# %%-- Load data
adata = AllDataCenters.load_singlecell(
    exclude_idh1_mut=True,
    exclude_patient_ids = True,
    exclude_relapse = True,
    )

# ensure raw counts are in .X. clear log1p status if present
if 'log1p' in adata.uns:
    del adata.uns['log1p']
adata.X = adata.layers['unnormalized_counts']

# create simplified annotations from lvl1 to lvl4
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"Mapping {lvl}")
    ct_mapping = (
        pd.read_csv(SIMPLIFIED_ANNOTATION_MAPPING)
        .set_index(REFERENCE_COLUMN)[f"celltype_{lvl}_curated"])
    adata.obs[f"celltype_{lvl}_curated"] = adata.obs[REFERENCE_COLUMN].astype(str).replace(ct_mapping)

# %% -- Initialize the InferCNVAnalyzer class
infercnv_obj = InferCNVAnalyzer(adata=adata)

# add a "reference" column to add "low_genes" and "high_counts" for downstream filtering
infercnv_obj.add_reference_column(
    cell_type_annotation_key=REFERENCE_COLUMN,
    threshold_n_genes=THRESHOLD_N_GENES,
    threshold_n_counts=THRESHOLD_N_COUNTS
)

# %% -- Run infercnv inference
start_time = time.time()
infercnv_obj.run_infercnv(list_reference_cat=REF_CELL_TYPES)

infercnv_obj.run_leiden_clustering(resolution=0.4)

end_time = time.time()
print(f"Total time running: {round((end_time-start_time)/60)} minutes")

# %% -- Assign each cnv cluster a CNV score
infercnv_obj.run_cnv_score()

# %% -- Run ITH scoring
start_time = time.time()
infercnv_obj.run_ith_scores(groupby=BATCH_KEY, label_key='celltype_level1_curated', label_to_score='Malignant', save_path=RESULTS_DIR / "ith_scores.csv")
end_time = time.time()
print(f"Total time running: {round((end_time-start_time)/60)} minutes")
# %% -- Perform visualizations
infercnv_obj.composite_plots(save_path=RESULTS_DIR / "cnv_plots.png")

# %% -- Create the annotation file
infercnv_obj.tumor_clone_annotation(save_path=RESULTS_DIR / "tumor_clone_annotation.csv")

# %% -- Save the summary heatmap
cnv.pl.chromosome_heatmap_summary(
    infercnv_obj.adata[~infercnv_obj.adata.obs[REFERENCE_COLUMN].isin(REF_CELL_TYPES+['low_genes','high_counts'])],
    groupby="cnv_leiden",
    dendrogram=True,
    show=False,
    vmax=0.38

)
fig = plt.gcf()
fig.savefig(RESULTS_DIR / "cnv_heatmap.png", dpi=300, bbox_inches='tight')
plt.close(fig)

cnv.pl.chromosome_heatmap_summary(
    infercnv_obj.adata,
    groupby="cnv_leiden",
    dendrogram=True,
    show=False,
    vmax=0.38
)
fig = plt.gcf()
fig.savefig(RESULTS_DIR / "cnv_heatmap_malignant.png", dpi=300, bbox_inches='tight')
plt.close(fig)

# %% -- Create gene expression UMAPs colored by cell type or by cluster
# infercnv_obj.plot_gene_exp_umaps(
#     cell_type_lvl="celltype_level2_scanvi",
#     save_path=RESULTS_DIR / "umaps_lvl2_CNVcluster.png",
#     reset_umap=False # reset to False the second time
# )

# %% -- Recompute genomic window and derive start-end position
infercnv_obj.get_bin_positions(
    save_path=RESULTS_DIR / "genomic_windows.csv"
)

# %% Save this temporary object
infercnv_obj.adata.layers.clear()
infercnv_obj.adata.X = scipy.sparse.csr_matrix(infercnv_obj.adata.X.shape)
infercnv_obj.adata.write(RESULTS_DIR / "infercnv_adata.h5ad")

########################################################
############## PART 2 POST-HOC ANALYSES ##############
########################################################
# %% -- Reload the adata object
SC_ADATA_PATH = RESULTS_DIR / "infercnv_adata.h5ad"
single_cell_obj = ad.read_h5ad(SC_ADATA_PATH)
celltype_col = "celltype_level1_curated"

# -- Drop the "doublet" compartment before clustering
single_cell_obj_filt = single_cell_obj[~single_cell_obj.obs[celltype_col].isin(['Doublet', 'Doublet_TNK_nos'])]
# drop also other_low_counts in level 4
single_cell_obj_filt = single_cell_obj_filt[~single_cell_obj_filt.obs["celltype_level4_scanvi"].isin(["Other_LowCounts"])]
single_cell_obj_filt.obs[celltype_col].value_counts()

# %% -- Rerun the functions to create the new neighbooring and clusters
cnv.tl.pca(single_cell_obj_filt)
cnv.pp.neighbors(single_cell_obj_filt)
cnv.tl.leiden(single_cell_obj_filt, resolution=0.4,flavor="igraph", n_iterations=2)

# %% -- Need to do some processing to reconstruct the dendogram
# 1) Ensure cnv_leiden is a clean categorical with no unused categories
g = "cnv_leiden"
single_cell_obj_filt.obs[g] = single_cell_obj_filt.obs[g].astype("category")
single_cell_obj_filt.obs[g] = single_cell_obj_filt.obs[g].cat.remove_unused_categories()

# 2) Remove stale cached dendrogram + palette for this groupby
dendro_key = f"dendrogram_{g}"
color_key = f"{g}_colors"

for k in [dendro_key, color_key]:
    if k in single_cell_obj_filt.uns:
        del single_cell_obj_filt.uns[k]

# 3) Recompute dendrogram with previous settings
sc.tl.dendrogram(
    single_cell_obj_filt,
    groupby=g,
    use_rep="X_cnv_pca",
    cor_method="pearson",
    linkage_method="complete",
    optimal_ordering=True,
)

# 4) Rerun the figure
cnv.pl.chromosome_heatmap_summary(
    single_cell_obj_filt,
    groupby=g,
    dendrogram=True,
    show=False,
    vmax=0.38,
)
# We still see weak signals from 3 clusters
# We would like to create a new visualization after
# removal of them, but we need to correctly do that

# %% -- Use reference cells to set a threshold for malignancy
REF_CELL_TYPES = [
    'Astrocyte_homeostatic',
    'B_nos',
    'Plasma_nos',
    'Endothelial_cap','Endothelial_tip',
    'Neutrophils_nos',
    'Microglia',
    'Monocyte_nos',
    'Pericyte',
    'Neuron','Neuron_exc_CDH24','Neuron_exc_nos','Neuron_neurogliaform',
    'Oligodendrocyte_nos',
    'NK_KLRC1','NK_nos',
    'T_CD4_nos','T_CD8_nos',
]
ref_mask = single_cell_obj_filt.obs["celltype_level4_scanvi"].isin(REF_CELL_TYPES)

iCNV_threshold = np.percentile(
    single_cell_obj_filt.obs.loc[ref_mask, "cnv_score"],
    99
)

# iCNV status = inferCNV-derived CNV burden class
single_cell_obj_filt.obs["iCNV_status"] = np.where(
    single_cell_obj_filt.obs["cnv_score"] > iCNV_threshold,
    "iCNV-high",
    "iCNV-low"
)

# %% -- Plot the heatmap with inferCNV burden classes (optional quick check)
cnv.pl.chromosome_heatmap(
    single_cell_obj_filt,
    groupby="iCNV_status",
    use_rep="cnv",
    figsize=(16, 3),
    cmap="bwr"
)

# %% -- Cross-tab lineage (level 1 curated) vs inferCNV burden class
# (No claim of malignancy here: this is purely iCNV signal vs lineage.)
ctab = pd.crosstab(
    single_cell_obj_filt.obs["celltype_level1_curated"],
    single_cell_obj_filt.obs["iCNV_status"],
    normalize="index"
)
ctab
# We see that most healthy cells have more than 99% non-malignant CNV
# But about 60% of tumor-lineage cells are CNV-low.
# In GBM this is expected and biologically meaningful
# These cells are:
# # tumor lineage by transcriptome
# # but low-amplitude / focal / subclonal in CNV space
# # often infiltrative, stem-like, or cycling slowly
# So a global CNV threshold is under-calling real malignant cells and
# shouldn't be used alone

# %% -- Define transcriptomic lineage status (tumor-lineage vs non-tumor lineage)
tumor_types = ["Malignant", "Other"]  # you said "Other" contains malignant
single_cell_obj_filt.obs["lineage_status"] = np.where(
    single_cell_obj_filt.obs["celltype_level1_curated"].isin(tumor_types),
    "tumor-lineage",
    "non-tumor lineage"
)

# Store numeric score too (kept as-is)
single_cell_obj_filt.obs["iCNV_score"] = single_cell_obj_filt.obs["cnv_score"]

# Combined group (for summaries/plots)
single_cell_obj_filt.obs["lineage_iCNV_group"] = (
    single_cell_obj_filt.obs["lineage_status"].astype(str)
    + " / "
    + single_cell_obj_filt.obs["iCNV_status"].astype(str)
).astype("category")

# Optional: enforce a consistent order in plots
desired_order = [
    "non-tumor lineage / iCNV-low",
    "non-tumor lineage / iCNV-high",
    "tumor-lineage / iCNV-low",
    "tumor-lineage / iCNV-high",
]
present = [x for x in desired_order if x in single_cell_obj_filt.obs["lineage_iCNV_group"].cat.categories]
single_cell_obj_filt.obs["lineage_iCNV_group"] = single_cell_obj_filt.obs["lineage_iCNV_group"].cat.reorder_categories(
    present, ordered=True
)

# %% -- Heatmap summary by combined group
save_path = RESULTS_DIR / "heatmap_summary_lineage_iCNV_group.pdf"
# ---- Add n and % to group labels for summary heatmap ----
group_col = "lineage_iCNV_group"

grp = single_cell_obj_filt.obs[group_col].astype(str)
counts = grp.value_counts()
pct = counts / counts.sum() * 100

label_map = {
    g: f"{g}\n(n={counts[g]:,}, {pct[g]:.1f}%)"
    for g in counts.index
}

single_cell_obj_filt.obs["lineage_iCNV_group_labeled"] = (
    grp.map(label_map).astype("category")
)

# Preserve group order
orig_order = single_cell_obj_filt.obs[group_col].cat.categories
new_order = [label_map[g] for g in orig_order if g in label_map]

single_cell_obj_filt.obs["lineage_iCNV_group_labeled"] = (
    single_cell_obj_filt.obs["lineage_iCNV_group_labeled"]
    .cat.reorder_categories(new_order, ordered=True)
)

axs = cnv.pl.chromosome_heatmap_summary(
    single_cell_obj_filt,
    groupby="lineage_iCNV_group_labeled",
    figsize=(10, 5),
    show=False,
)
ax = max(axs.values(), key=lambda a: len(a.get_yticklabels()))
fig = ax.figure
ax.tick_params(axis="y", labelsize=12)
ax.tick_params(axis="x", labelsize=10)
ax.set_ylabel("")

for a in fig.axes:
    if a is not ax:
        a.tick_params(labelsize=10)

fig.savefig(save_path, bbox_inches="tight")
plt.close(fig)

# %% -- Heatmap not summarize to see the size of the clusters
cnv.pl.chromosome_heatmap(
    single_cell_obj_filt,
    groupby="lineage_iCNV_group"
)

# %% -- Distribution plot (KDE) of inferCNV CNV score by combined group
fig, ax = plt.subplots(figsize=(7, 3.6), dpi=300)
sns.kdeplot(
    data=single_cell_obj_filt.obs,
    x="iCNV_score",
    hue="lineage_iCNV_group",
    fill=True,
    common_norm=False,
    alpha=0.35,
    linewidth=2,
    ax=ax
)
ax.axvline(iCNV_threshold, color="black", linestyle="--", linewidth=1)
ax.set_xlabel("inferCNV CNV score (iCNV score)")
ax.set_ylabel("Density")
sns.despine(ax=ax)
plt.tight_layout()

save_path = RESULTS_DIR / "kde_iCNV_score.pdf"
plt.savefig(save_path, format="pdf", bbox_inches="tight")
plt.show()

# %% -- Sankey plot showing the counts of cell types in each
# group
# params
LEVEL2_COL = "celltype_level2_scanvi"
GROUP_COL  = "lineage_iCNV_group"

MIN_LEVEL2_CELLS = 30          # collapse rare level 2 categories

obs = single_cell_obj_filt.obs.copy()

# Sanity checks
for col in [LEVEL2_COL, GROUP_COL]:
    if col not in obs.columns:
        raise KeyError(f"{col} not found in obs.columns")

# Collapse rare level 2 categories for readability (robust to categoricals)
lvl2 = obs[LEVEL2_COL].astype(str)

lvl2_counts = lvl2.value_counts()
rare_lvl2 = lvl2_counts[lvl2_counts < MIN_LEVEL2_CELLS].index

obs["level2_sankey"] = lvl2.where(~lvl2.isin(rare_lvl2), other="Other_level2")

# Build contingency table
ct = (
    obs[["level2_sankey", GROUP_COL]]
    .value_counts()
    .reset_index(name="value")
)

# Nodes
sources = ct["level2_sankey"].unique().tolist()
targets = ct[GROUP_COL].astype(str).unique().tolist()
labels = sources + targets
idx = {lab: i for i, lab in enumerate(labels)}

# Fix node positions to reduce horizontal stretch
x_sources = [0.05] * len(sources)
x_targets = [0.55] * len(targets)
x_positions = x_sources + x_targets

y_positions = None  # let Plotly stack them automatically

# Sankey plot
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        label=labels,
        x=x_positions,
        pad=12,
        thickness=12,
        line=dict(width=0.5)
    ),
    link=dict(
        source=ct["level2_sankey"].map(idx).values,
        target=ct[GROUP_COL].astype(str).map(idx).values,
        value=ct["value"].values
    )
)])

fig.update_layout(
    font_size=14,
    width=800,
    height=500,
    margin=dict(l=10, r=10, t=40, b=10),
)

fig.update_traces(
    node=dict(thickness=16, pad=14)
)

# Always works: interactive HTML
html_path = RESULTS_DIR / "sankey_level2_scanvi_to_lineage_iCNV_group.html"
fig.write_html(str(html_path))

# Optional static exports (requires kaleido)
for ext in ["pdf", "svg", "png"]:
    out = RESULTS_DIR / f"sankey_level2_scanvi_to_lineage_iCNV_group.{ext}"
    try:
        fig.write_image(str(out))
    except Exception as e:
        print(f"Could not export {ext}: {e}")

fig.show()
print("Saved:", html_path)

# %% -- Now treat only the malignant cells for visualization
single_cell_obj_filt.obs["tumor_lineage"] = single_cell_obj_filt.obs["lineage_status"] == "tumor-lineage"
adata_tumor = single_cell_obj_filt[
    single_cell_obj_filt.obs["tumor_lineage"]
].copy()

# %% -- Rerun the functions to create the new neighbooring and clusters
cnv.tl.pca(adata_tumor)
cnv.pp.neighbors(adata_tumor)
cnv.tl.leiden(adata_tumor, resolution=0.4,flavor="igraph", n_iterations=2)

# %% -- Need to do some processing to reconstruct the dendogram
# 1) Ensure cnv_leiden is a clean categorical with no unused categories
g = "cnv_leiden"
adata_tumor.obs[g] = adata_tumor.obs[g].astype("category")
adata_tumor.obs[g] = adata_tumor.obs[g].cat.remove_unused_categories()

# 2) Remove stale cached dendrogram + palette for this groupby
dendro_key = f"dendrogram_{g}"
color_key = f"{g}_colors"

for k in [dendro_key, color_key]:
    if k in adata_tumor.uns:
        del adata_tumor.uns[k]

# 3) Recompute dendrogram with previous settings
sc.tl.dendrogram(
    adata_tumor,
    groupby=g,
    use_rep="X_cnv_pca",
    cor_method="pearson",
    linkage_method="complete",
    optimal_ordering=True,
)

# 4) Rerun the figure
cnv.pl.chromosome_heatmap_summary(
    adata_tumor,
    groupby=g,
    dendrogram=True,
    show=False,
    vmax=0.38,
    save = RESULTS_DIR / "cnv_heatmap_malignant.pdf"
)

# %% -- Make some statistics for the results section
# 1) Calculate how many clusters and the proportion of each in the total
# population
# cluster sizes
cluster_counts = adata_tumor.obs["cnv_leiden"].value_counts().sort_index()
cluster_props = cluster_counts / cluster_counts.sum()

cluster_stats = pd.DataFrame({
    "n_cells": cluster_counts,
    "fraction": cluster_props,
    "percent": (cluster_props * 100).round(2)
})

cluster_stats.sort_values(by="percent", ascending=False)

# %% -- Calculate CNV score on the CNV profile
# Generate UMAP
cnv.tl.umap(adata_tumor)
cnv.tl.cnv_score(adata_tumor, groupby="cnv_leiden")

# %% -- Create composite figure
save_path = RESULTS_DIR / "cnv_scores.pdf"

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(12, 5))

# Plot UMAPs into specific axes
cnv.pl.umap(
    adata_tumor,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False
)
cnv.pl.umap(adata_tumor, color="cnv_score", ax=ax2, show=False)

# Set plot params
ax1.set_title("CNV profiles UMAP colored by CNV leiden clusters")
ax2.set_title("CNV profiles UMAP colored by CNV scores")
ax1.title.set_size(16)
ax2.title.set_size(16)
ax1.set_xlabel("UMAP 1", fontsize=14)
ax1.set_ylabel("UMAP 2", fontsize=14)
ax2.set_xlabel("UMAP 1", fontsize=14)
ax2.set_ylabel("UMAP 2", fontsize=14)

# Adjust the layout
fig.tight_layout()

# Save the figure if a path is provided
fig.savefig(save_path, dpi=300)

# %% -- Save again the figures
# Composite plots
save_path = RESULTS_DIR / "cnv_plots.pdf"

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
# Plot normal umap colored by the cnv score or leiden clusters
sc.pl.umap(
    adata_tumor,
    color="cnv_score",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
sc.pl.umap(adata_tumor, color="cnv_leiden", ax=ax2, show=False)

# Plot UMAPs into specific axes
cnv.pl.umap(adata_tumor, color="cnv_score", ax=ax3, show=False)

cnv.pl.umap(
    adata_tumor,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax4,
    show=False,
)

# Set the title for each subplot
ax1.set_title("Gene expression UMAP colored by CNV score")
ax2.set_title("Gene expression UMAP colored by CNV leiden clusters")
ax3.set_title("CNV profiles UMAP colored by CNV scores")
ax4.set_title("CNV profiles UMAP colored by CNV leiden clusters")

# Reduce the title size
ax1.title.set_size(16)
ax2.title.set_size(16)
ax3.title.set_size(16)
ax4.title.set_size(16)

# modify axes labels for ax3 and ax4
ax3.set_xlabel("UMAP 1", fontsize=14)
ax3.set_ylabel("UMAP 2", fontsize=14)
ax4.set_xlabel("UMAP 1", fontsize=14)
ax4.set_ylabel("UMAP 2", fontsize=14)

# Adjust the layout
fig.tight_layout()

# Save the figure if a path is provided
fig.savefig(save_path, dpi=300)

# %% -- Calculate some statistics
# count the number of cells in each cluster per sample
counts_df = adata_tumor.obs.groupby(
    ["orig.ident", "cnv_leiden"]
).size().unstack().fillna(0)

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 14,
    "legend.title_fontsize": 14
})

# clustermap
g = sns.clustermap(
    counts_df,
    cmap="viridis",
    figsize=(17, 15),
    z_score=0,
    row_cluster=True,
    col_cluster=True,
)
g.ax_heatmap.set_xlabel("CNV-based clonal clusters")
g.ax_heatmap.set_ylabel("Sample IDs")
g.ax_heatmap.set_yticklabels([])     # masque les IDs des samples

# save
g.savefig(
    RESULTS_DIR / "clustermap_cnv_clusters.pdf",
    bbox_inches="tight"
)

# %% -- Deep dive clusters 9 to 14 very patient specific
patient_frac = (
    adata_tumor.obs
    .groupby(["cnv_leiden", "orig.ident"])
    .size()
    .groupby(level=0)
    .apply(lambda x: x / x.sum())
)
top_patient_frac = patient_frac.groupby("cnv_leiden").max()

# %% -- Heterogeneity score per sample
# calculate an "heterogeneity" score that will show how many cluster clone each sample has
heterogeneity_score = counts_df.apply(
    lambda x: len(x[x > 0]), axis=1
)

# plot the heterogeneity_score
mpl.rcParams.update(mpl.rcParamsDefault)

fig, ax = plt.subplots(figsize=(4.8, 2.5), dpi=300)

sns.violinplot(
    x=heterogeneity_score,
    orient="h",
    color="lightgray",
    inner="quartile",
    linewidth=1,
    ax=ax
)

# labels
ax.set_xlabel("Number of CNV-based clones per tumor sample", fontsize=11)
ax.set_ylabel("GBM tumor samples", fontsize=11)

# y axis minimal mais visible
ax.set_yticks([0])
ax.set_yticklabels([""])

# option : médiane explicite
median_val = np.median(heterogeneity_score)
ax.axvline(
    median_val,
    color="black",
    linestyle="--",
    linewidth=1
)

sns.despine(ax=ax, top=True, right=True)
plt.tight_layout()

plt.savefig(RESULTS_DIR / "cnv_clonal_heterogeneity_violin.pdf", bbox_inches="tight")
plt.show()

# %% -- Cross cnv scores with level 4 annotations
df_subset = adata_tumor.obs[["celltype_level4_scanvi", "cnv_score"]].copy()
df_subset["cnv_score"] = pd.to_numeric(df_subset["cnv_score"], errors="coerce")
df_subset = df_subset.dropna(subset=["celltype_level4_scanvi", "cnv_score"])

print("N cells used:", len(df_subset))
print("N states:", df_subset["celltype_level4_scanvi"].nunique())

state_cnv = (
    df_subset.groupby("celltype_level4_scanvi")["cnv_score"]
      .agg(["mean", "median", "std", "count"])
      .sort_values("mean")
)

plt.figure(figsize=(9, 4), dpi=300)
ax = sns.violinplot(
    data=df_subset,
    x="celltype_level4_scanvi",
    y="cnv_score",
    inner="quartile",
    cut=0
)
ax.set_xlabel("")
ax.set_ylabel("CNV score")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()

# %% -- Cross level 4 annotaiton of tumors with clonality
ct_label_cluster = pd.crosstab(
    adata_tumor.obs["cnv_leiden"],
    adata_tumor.obs["celltype_level4_scanvi"],
    normalize="index"
)

dominant_state = ct_label_cluster.idxmax(axis=1)
dominant_fraction = ct_label_cluster.max(axis=1)

dominant_df = (
    pd.DataFrame({
        "dominant_state": dominant_state,
        "fraction": dominant_fraction
    })
    .sort_index()
)

print(dominant_df)

plt.figure(figsize=(10, 6), dpi=300)

sns.heatmap(
    ct_label_cluster,
    cmap="viridis",
    vmin=0,
    vmax=1,
    cbar_kws={"label": "Fraction of cells"}
)

plt.xlabel("Transcriptional state (scanvi)")
plt.ylabel("CNV-based clonal cluster")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()

# %% -- Save other stuffs
print("Computing ITH scores based on CNAs")
cnv.tl.ithcna(adata_tumor, groupby="orig.ident", use_rep='X_cnv', key_added='ithcna', inplace=True)
print("Computing ITH scores based on GEX")
cnv.tl.ithcna(adata_tumor, groupby="orig.ident", key_added='ithgex', inplace=True)

save_path = RESULTS_DIR / "ith_scores_clean"
adata_tumor.obs[['ithcna','ithgex']].to_csv(save_path)

save_path = RESULTS_DIR / "tumor_clone_annotation_clean.csv"
df_tumor_clones = pd.DataFrame(
    adata_tumor.obs[["cnv_leiden", "cnv_score"]].copy()
)
df_tumor_clones.to_csv(save_path, index=True)

# %% Save this latest adata containing clean (as much as possible) malignant cells
adata_tumor.layers.clear()
adata_tumor.X = scipy.sparse.csr_matrix(adata_tumor.X.shape)
adata_tumor.write(RESULTS_DIR / "infercnv_adata_malignant_only.h5ad")
