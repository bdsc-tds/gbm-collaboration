"""Fig 3d and the PERMANOVA numbers: WES CNV clusters vs ECM neighborhoods and cell types.

Produces:
- The PERMANOVA numbers of the Results / Methods: ECM neighborhood composition vs CNV
  cluster (pseudo-F = 0.58, p = 0.86), cell-type composition vs CNV cluster (pseudo-F =
  4.56, p = 0.001).
- Fig 3d: mean proportion of the 15 level-4 cell types with the lowest FDR per CNV cluster,
  stars from a Kruskal-Wallis test across the four clusters, BH FDR over all level-4 types
  (* < 0.05, ** < 0.01, *** < 0.001). The original script saved it as figures/figure3/fig3e.pdf; its
  table (fig3e.csv) is the published Fig 3d's, and is compared below. The figure legend's
  "T-test of significant association vs total cohort" is not what ran.

Origin: analysis/multimodal_correlation/correlation_cnv_wes_ecm.py, Part 1 (ECM) and Part 2
(cellular landscape); only the outputs of the paper are kept. Inputs (frozen): the WES
cluster labels (results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv), the ECM
neighborhood labels per spot (Analysis_Repo/Notebooks/ecm_neighborhood_annotations_v4.csv)
and the single-cell obs (level-4 scANVI labels).

Comparison with the published values (printed):
- The Fig 3d table is identical to the published figures/figure3/fig3e.csv.
- The code runs 999 permutations (Methods: 9,999); the p values depend on the permutations,
  the pseudo-F statistics do not. The PERMANOVA is computed here with numpy/scipy (scikit-bio
  is not in the environment), same pseudo-F, p from a seeded permutation.
- Cell-type PERMANOVA: the Methods say level-4 proportions; the published pseudo-F 4.56 is
  the one of the level-3 curated proportions (27 types). The level-4 table (the one of
  Fig 3d) gives 4.02 (p = 0.001 for both). Both are computed below.
"""

# %% -- Load packages
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smm
from scipy.spatial.distance import pdist, squareform
from scipy.stats import kruskal

from gbmcollaboration import AllDataCenters
from gbmcollaboration.paper import output

RESULTS_DIR = output.setup("fig3d_correlation_cnv_wes_ecm")

PERMUTATIONS = 999
PUBLISHED_FIG3D_TABLE = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/figure3/fig3e.csv"


# %% -- Functions
def permanova(X, grouping, permutations: int = 999, seed: int = 0) -> pd.Series:
    """PERMANOVA on the Bray-Curtis dissimilarities of the rows of X, as scikit-bio's
    `permanova(beta_diversity("braycurtis", X, ids), grouping, permutations)`.

    The pseudo-F statistic is deterministic; the p value depends on the permutations
    (seeded here, scikit-bio used the global numpy state).
    """
    d2 = squareform(pdist(np.asarray(X, dtype=float), metric="braycurtis")) ** 2
    grouping = np.asarray(grouping)
    n = len(grouping)
    _, codes = np.unique(grouping, return_inverse=True)
    n_groups = codes.max() + 1
    s_t = d2[np.triu_indices(n, 1)].sum() / n

    def pseudo_f(c):
        sizes = np.bincount(c)
        s_w = sum(np.triu(d2[np.ix_(c == k, c == k)], 1).sum() / sizes[k] for k in range(n_groups))
        return ((s_t - s_w) / (n_groups - 1)) / (s_w / (n - n_groups))

    f_obs = pseudo_f(codes)
    rng = np.random.default_rng(seed)
    f_perm = np.array([pseudo_f(rng.permutation(codes)) for _ in range(permutations)])
    p = ((f_perm >= f_obs).sum() + 1) / (permutations + 1)
    return pd.Series({"method name": "PERMANOVA", "test statistic name": "pseudo-F",
                      "sample size": n, "number of groups": n_groups,
                      "test statistic": f_obs, "p-value": p, "number of permutations": permutations})


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


# %% -- Load data
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

# %% -- Permanova analysis (Methods: pseudo-F = 0.58, p = 0.86)
ecm_cols = [c for c in data.columns if c != "cnv_leiden"]
ecm_matrix = data[ecm_cols].copy()
grouping = data["cnv_leiden"].copy()

res = permanova(ecm_matrix.to_numpy(), grouping=grouping.to_numpy(), permutations=PERMUTATIONS)
print(res)
permanova_results = {"ECM neighborhoods": res}

# no association between WES CNA profiles and ECM neighborhood

####################################################################
## Part 2 : correlation with cellular landscape ####################
####################################################################
# %% Params
SIMPLIFIED_ANNOTATION_MAPPING = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/SingleCell_Annotation_Key_v3.csv")
REFERENCE_COLUMN = "celltype_level4_scanvi"
BATCH_KEY = "orig.ident"


# %% -- Load single data (obs only)
adata = ad.AnnData(obs=AllDataCenters.load_singlecell(
    exclude_idh1_mut=True,
    exclude_patient_ids = True,
    exclude_relapse = True,
    obs_only=True,
    ))

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
print(f"{single_cell_obj_filt.n_obs:,} nuclei from {single_cell_obj_filt.obs['orig.ident'].nunique()} samples")

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

# %% -- permanova analysis (Methods: pseudo-F = 4.56, p = 0.001)
cell_cols = [c for c in data_cells.columns if c != "cnv_leiden"]

X = data_cells[cell_cols].to_numpy()
grouping = data_cells["cnv_leiden"].to_numpy()

res = permanova(X, grouping=grouping, permutations=PERMUTATIONS)
print(res)
permanova_results["Level-4 cell types"] = res

# level-3 curated proportions (published pseudo-F 4.56)
counts_l3 = pd.crosstab(single_cell_obj_filt.obs["orig.ident"].astype(str),
                        single_cell_obj_filt.obs["celltype_level3_curated"].astype(str))
props_l3 = counts_l3.div(counts_l3.sum(axis=1), axis=0).loc[data_cells.index]
permanova_results["Level-3 cell types"] = permanova(props_l3.to_numpy(), grouping=grouping,
                                                    permutations=PERMUTATIONS)
permanova_results = pd.DataFrame(permanova_results).T
permanova_results["published pseudo-F"] = [0.58, 4.56, 4.56]
permanova_results["published p"] = [0.86, 0.001, 0.001]
print(permanova_results[["sample size", "test statistic", "p-value", "published pseudo-F", "published p"]])
output.save_table(RESULTS_DIR, "permanova_cnv_cluster", permanova_results)

# %% -- Kruskal Wallis test
results = []
cell_cols = [c for c in data_cells.columns if c != "cnv_leiden"]

for col in cell_cols:
    groups = [data_cells[data_cells["cnv_leiden"]==g][col]
              for g in np.unique(data_cells["cnv_leiden"])]

    stat, p = kruskal(*groups)
    results.append((col, stat, p))

res_df = pd.DataFrame(results, columns=["cell_type","H","p"])
res_df["p_adj"] = smm.multipletests(res_df["p"], method="fdr_bh")[1]
res_df = res_df.sort_values(
    ["p_adj", "p", "H"],
    ascending=[True, True, False],
    kind="mergesort"
)

print(res_df.head(20))
output.save_table(RESULTS_DIR, "kruskal_wallis_test_cellular_landscape_cluster", res_df)

published_kw = pd.read_csv(PUBLISHED_FIG3D_TABLE, index_col=0)
check = res_df.reset_index(drop=True).merge(published_kw, on="cell_type", suffixes=("", "_published"))
print(f"Kruskal–Wallis table vs the published figure3/fig3e.csv: {len(check)} of {len(published_kw)} types matched; "
      f"max |ΔH| = {(check['H'] - check['H_published']).abs().max():.2e}; same order: "
      f"{res_df['cell_type'].tolist() == published_kw['cell_type'].tolist()}")

# %% -- top cell types (by FDR)
top_n = 15
top_types = res_df.head(top_n)["cell_type"].tolist()

# %% -- dot plot with significance (Fig 3d)
# -- colors matching fig2b
cluster_colors = {
    0: "tab:blue",
    1: "tab:orange",
    2: "tab:green",
    3: "tab:red",
}

# -- exclude non-cell-type columns if present
exclude = {"cnv_leiden"}
cell_cols = [c for c in data_cells.columns if c not in exclude]

# -- mean proportions per cluster (force cluster order 0,1,2,3)
mean_by = data_cells.groupby("cnv_leiden")[cell_cols].mean().reindex([0, 1, 2, 3])

# -- select top cell types by FDR (from univariate table)
top_types = [ct for ct in top_types if ct in mean_by.columns]

# -- build plot dataframe: rows=cell types, cols=clusters
plot_df = mean_by[top_types].T
# order rows for readability (e.g., by cluster 2 mean)
plot_df = plot_df.sort_values(by=2, ascending=True)

# -- map adjusted p-values onto displayed cell types
p_map = res_df.set_index("cell_type")["p_adj"]
plot_df["p_adj"] = plot_df.index.map(p_map)

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
ax.set_title(f"Cell-type composition across CNA clusters (top {len(top_types)} by FDR and effect size)", pad=8)

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
output.save(RESULTS_DIR, "fig3d_celltype_composition_cna_clusters", fig, res_df, panel="Fig 3d")
plt.show()
plt.close(fig)
print(plot_df.assign(stars=plot_df["p_adj"].map(stars)).round(4).to_string())

# %%
print(output.versions(["numpy", "pandas", "scipy", "statsmodels", "anndata", "matplotlib"]).to_string(index=False))
