"""Script to associate WES cluster with ECM neighborhood or cell landscape."""

# %% -- Load packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.stats.distance import permanova
from skbio.diversity import beta_diversity
from scipy.stats import kruskal
import statsmodels.stats.multitest as smm
import anndata as ad
from pathlib import Path
from gbmcollaboration import AllDataCenters

RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/figure3")

# -- Load data
cnv_wes_cluster_path = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
ecm_cluster_path = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/Analysis_Repo/Notebooks/ecm_neighborhood_annotations_v4.csv"

wes_cluster = pd.read_csv(cnv_wes_cluster_path)
ecm_cluster = pd.read_csv(ecm_cluster_path)
ecm_cluster = ecm_cluster.set_index("spot_id")

print(wes_cluster.head(2))
print(ecm_cluster.head(2))

####################################################################
## Part 1 : correlation with ECM neighborhood ######################
####################################################################
# %%
# -- spot-level dataframe
ecm_df = ecm_cluster.copy()
# Processing of sample ids
ecm_df["sample"] = ecm_df["sample"].astype(str).str.strip()
ecm_df["sample"] = ecm_df["sample"].str[:9]

ecm_df["ECM_neighborhood"] = pd.to_numeric(ecm_df["ECM_neighborhood"], errors="raise").astype(int)

# -- counts patients x ECM
all_ecm = np.sort(ecm_df["ECM_neighborhood"].unique())
counts = (
    pd.crosstab(
        index=ecm_df["sample"],
        columns=pd.Categorical(ecm_df["ECM_neighborhood"], categories=all_ecm),
        dropna=False
    )
    .sort_index()
)

# -- proportions
props = counts.div(counts.sum(axis=1), axis=0)

print("counts shape:", counts.shape)
print("props shape:", props.shape)
print("ECM columns:", props.columns.tolist())
print(props.head())

# -- Processing of wes_cluster
wes = wes_cluster.copy()
wes = wes.rename(columns={"orig.ident": "sample"})
wes["sample"] = wes["sample"].astype(str).str.strip()

wes = wes.drop_duplicates(subset="sample", keep="first")

# -- inner join
data = (
    props
    .merge(wes[["sample", "cnv_leiden"]], left_index=True, right_on="sample", how="inner")
    .set_index("sample")
)

# -- checks
assert data.index.is_unique, "Erreur: samples non uniques"
assert data["cnv_leiden"].notna().all(), "Erreur: cnv_leiden contient des NA"

print("data shape:", data.shape)
print(data.head())

# %% -- Permanova analysis
ecm_cols = [c for c in data.columns if c != "cnv_leiden"]
ecm_matrix = data[ecm_cols].copy()
grouping = data["cnv_leiden"].copy()

# -- force index in simple strings
ecm_matrix.index = ecm_matrix.index.astype(str)
grouping.index = grouping.index.astype(str)

# -- check alignement
assert ecm_matrix.index.equals(grouping.index)

ids = ecm_matrix.index.to_list()

dm = beta_diversity(
    metric="braycurtis",
    counts=ecm_matrix.to_numpy(),
    ids=ids
)

res = permanova(dm, grouping=grouping.to_numpy(), permutations=999)
print(res)

# %% -- Krukal Wallis test
results = []

for col in ecm_cols:
    groups = [ecm_matrix[grouping == g][col] for g in np.unique(grouping)]
    stat, p = kruskal(*groups)
    results.append((col, stat, p))

res_df = pd.DataFrame(results, columns=["ECM", "H", "p"])
res_df["p_adj"] = smm.multipletests(res_df["p"], method="fdr_bh")[1]
res_df.sort_values("p_adj", inplace=True)

print(res_df)

# %% -- Mean proportions
mean_props = data.groupby("cnv_leiden")[ecm_cols].mean()
print(mean_props)

# no association between WES CNA profiles and ECM neighborhood

####################################################################
## Part 2 : correlation with cellular landscape ####################
####################################################################
# %% Params
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

# %% -- Load single data and do small processing as for inferCNV analyses
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

# -- drop the "doublet" compartment before running initial analysis
single_cell_obj_filt = adata[~adata.obs["celltype_level1_curated"].isin(['Doublet', 'Doublet_TNK_nos'])]
# drop also other_low_counts in level 4
single_cell_obj_filt = single_cell_obj_filt[~single_cell_obj_filt.obs["celltype_level4_scanvi"].isin(["Other_LowCounts"])]

# %% -- extract metadata cellule
cell_df = single_cell_obj_filt.obs[["orig.ident", "celltype_level4_curated"]].copy()
cell_df.rename(columns={"orig.ident":"sample", "celltype_level4_curated":"cell_type"}, inplace=True)
cell_df["sample"] = cell_df["sample"].astype(str).str.strip()
cell_df["cell_type"] = cell_df["cell_type"].astype(str)

# -- counts patients x cell types
counts = pd.crosstab(
    cell_df["sample"],
    cell_df["cell_type"]
)

# -- proportions
props = counts.div(counts.sum(axis=1), axis=0)

print("Shape:", props.shape)
print(props.head())

# %% -- prepare WES cluster data and merge
wes = wes_cluster.rename(columns={"orig.ident": "sample"}).copy()
wes["sample"] = wes["sample"].astype(str).str.strip()
wes = wes.drop_duplicates("sample")

data_cells = (
    props
    .merge(wes[["sample", "cnv_leiden"]],
           left_index=True,
           right_on="sample",
           how="inner")
    .set_index("sample")
)

print(data_cells.shape)

# %% -- permanova analysis
cell_cols = [c for c in data_cells.columns if c != "cnv_leiden"]

X = data_cells[cell_cols].to_numpy()
grouping = data_cells["cnv_leiden"].to_numpy()

dm = beta_diversity(
    metric="braycurtis",
    counts=X,
    ids=data_cells.index.to_list()
)

res = permanova(dm, grouping=grouping, permutations=999)
print(res)

# %% -- Kruskal Wallis test
results = []

for col in cell_cols:
    groups = [data_cells[data_cells["cnv_leiden"]==g][col]
              for g in np.unique(data_cells["cnv_leiden"])]

    stat, p = kruskal(*groups)
    results.append((col, stat, p))

res_df = pd.DataFrame(results, columns=["cell_type","H","p"])
res_df["p_adj"] = smm.multipletests(res_df["p"], method="fdr_bh")[1]
res_df.sort_values("p_adj", inplace=True)

print(res_df)

# %% -- Stacked barplot
# top cell types (by FDR)
top_n = 15
top_types = res_df.sort_values("p_adj").head(top_n)["cell_type"].tolist()

# mean proportions per CNA cluster
mean_props = data_cells.groupby("cnv_leiden")[top_types].mean()
mean_props = mean_props.clip(lower=0)

# ensure cluster order 0,1,2,3 and labels
cluster_order = [0, 1, 2, 3]
mean_props = mean_props.reindex(cluster_order)

x_labels = [str(c) for c in mean_props.index]

# CNA cluster colors (as requested)
cluster_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
    3: "tab:red",
}

fig, ax = plt.subplots(figsize=(6.5, 3.2))
bottom = np.zeros(mean_props.shape[0])

for ct in mean_props.columns:
    vals = mean_props[ct].to_numpy()
    ax.bar(x_labels, vals, bottom=bottom, linewidth=0.6, edgecolor="white", label=ct)
    bottom += vals
for tick, c in zip(ax.get_xticklabels(), cluster_order):
    tick.set_color(cluster_colors[c])
    tick.set_fontweight("bold")

ax.set_ylabel("Mean proportion")
ax.set_ylim(0, 1.0)
ax.set_title(f"Mean cell-type composition per CNA cluster (top {top_n} by FDR)", pad=8)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="x", length=0)
ax.tick_params(axis="y", width=1)

ax.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False,
    title="Cell type",
    title_fontsize=9,
    fontsize=8,
)

plt.tight_layout()
save_path = RESULTS_DIR / "stacked_barplot_lvl4_top15.pdf"
plt.savefig(save_path, format="pdf", bbox_inches="tight")
plt.show()

# %% -- test if AC_OPC ratio is significantly represented
# -- only if lvl3 is used, otherwise need to regroup
data_cells["AC_OPC_ratio"] = (
    data_cells["Malignant_AC-like"] /
    (data_cells["Malignant_OPC-like"] + 1e-6)
)

groups = [data_cells[data_cells["cnv_leiden"]==g]["AC_OPC_ratio"]
          for g in sorted(data_cells["cnv_leiden"].unique())]
print(kruskal(*groups))
# yes, significant

# %% -- dot plot with significance

# -- colors matching fig2b
cluster_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
    3: "tab:red",
}

# -- exclude non-cell-type columns if present
exclude = {"cnv_leiden", "AC_OPC_ratio"}
cell_cols = [c for c in data_cells.columns if c not in exclude]

# -- mean proportions per cluster (force cluster order 0,1,2,3)
mean_by = data_cells.groupby("cnv_leiden")[cell_cols].mean().reindex([0, 1, 2, 3])

# -- select top cell types by FDR (from univariate table)
top_n = 15
top_types = res_df.sort_values("p_adj").head(top_n)["cell_type"].tolist()
top_types = [ct for ct in top_types if ct in mean_by.columns]

# -- build plot dataframe: rows=cell types, cols=clusters
plot_df = mean_by[top_types].T
# order rows for readability (e.g., by cluster 2 mean)
plot_df = plot_df.sort_values(by=2, ascending=True)

# -- map adjusted p-values onto displayed cell types
p_map = res_df.set_index("cell_type")["p_adj"]
plot_df["p_adj"] = plot_df.index.map(p_map)

# -- stars helper
def stars(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return ""

y = np.arange(plot_df.shape[0])

fig, ax = plt.subplots(figsize=(7.6, 5.6))

for i, ct in enumerate(plot_df.index):
    vals = plot_df.loc[ct, [0, 1, 2, 3]].to_numpy()
    ax.plot(vals, [i] * len(vals), color="lightgray", linewidth=2, zorder=1)
for cl in [0, 1, 2, 3]:
    ax.scatter(
        plot_df[cl].to_numpy(),
        y,
        s=55,
        color=cluster_colors[cl],
        label=f"CNA {cl}",
        zorder=2
    )

xmax = plot_df[[0, 1, 2, 3]].to_numpy().max()
text_x = xmax * 1.12
for i, (ct, row) in enumerate(plot_df.iterrows()):
    ax.text(
        text_x,
        i,
        stars(row["p_adj"]),
        va="center",
        fontsize=11,
        fontweight="bold",
        color="#333333"
    )

ax.set_yticks(y)
ax.set_yticklabels(plot_df.index)
ax.set_xlabel("Mean proportion")
ax.set_title(f"Cell-type composition across CNA clusters (top {len(top_types)} by FDR)", pad=8)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.tick_params(axis="x", width=1)

ax.set_xlim(0, xmax * 1.25)

ax.legend(
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0,
    title="CNA cluster"
)

plt.tight_layout()
save_path = RESULTS_DIR / "fig3d.pdf"
plt.savefig(save_path, format="pdf", bbox_inches="tight")
plt.show()

# %%
