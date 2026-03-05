"""Script to associate WES cluster with ECM neighborhood or cell landscape."""

# %% -- Load packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from skbio.stats.distance import permanova
#from skbio.diversity import beta_diversity
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

print(res_df)
res_df.to_csv("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/RNAseq/wes_cluster_rna/kruskal_wallis_test_cellular_landscape_cluster.csv")

# %% -- Stacked barplot
# top cell types (by FDR)
top_n = 15
top_types = res_df.head(top_n)["cell_type"].tolist()

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
save_path = RESULTS_DIR / "fig3e.pdf"
plt.savefig(save_path, format="pdf", bbox_inches="tight")
plt.show()

####################################################################
## Part 3 : pseudobulk single-cell analyses across WES cluster #####
####################################################################

# %% --- Utils
import scipy.sparse as sp

def _sum_counts_by_groups(X, group_codes, n_groups):
    """
    Efficient sum of sparse/dense matrix X (cells x genes) into (groups x genes),
    where group_codes is int array length n_cells in [0..n_groups-1].
    """
    if sp.issparse(X):
        out = sp.csr_matrix((n_groups, X.shape[1]), dtype=np.float64)
        for g in range(n_groups):
            idx = np.where(group_codes == g)[0]
            if idx.size == 0:
                continue
            out[g, :] = X[idx, :].sum(axis=0)
        return out
    else:
        out = np.zeros((n_groups, X.shape[1]), dtype=np.float64)
        for g in range(n_groups):
            idx = np.where(group_codes == g)[0]
            if idx.size == 0:
                continue
            out[g, :] = X[idx, :].sum(axis=0)
        return out


def make_pseudobulk_counts(
    adata,
    celltype_col="celltype_level4_curated",
    sample_col="orig.ident",
    target_celltypes=None,
    layer="unnormalized_counts",
    min_cells_per_sample=30,
):
    """
    Returns:
      pb_counts: dict[celltype] -> (counts_df: samples x genes) with raw counts
      pb_meta:   dict[celltype] -> (meta_df: samples x metadata columns)
    """
    if layer is None:
        X = adata.X
    else:
        X = adata.layers[layer]

    genes = adata.var_names.astype(str)

    if target_celltypes is None:
        target_celltypes = adata.obs[celltype_col].astype(str).unique().tolist()

    pb_counts = {}
    pb_meta = {}

    for ct in target_celltypes:
        sub = adata[adata.obs[celltype_col].astype(str) == str(ct)].copy()
        if sub.n_obs == 0:
            continue

        # counts matrix
        Xct = sub.layers[layer] if layer is not None else sub.X

        # group by sample
        samples = sub.obs[sample_col].astype(str).values
        sample_u, inv = np.unique(samples, return_inverse=True)
        n_groups = len(sample_u)

        # filter low-cell samples (min cells per sample)
        cell_counts = np.bincount(inv, minlength=n_groups)
        keep = cell_counts >= min_cells_per_sample
        if keep.sum() < 3:
            # not enough samples to run DESeq2
            continue

        # reindex groups after filtering
        kept_samples = sample_u[keep]
        # map inv -> new codes
        old_to_new = -np.ones(n_groups, dtype=int)
        old_to_new[np.where(keep)[0]] = np.arange(keep.sum())
        new_codes = old_to_new[inv]
        ok_cells = new_codes >= 0

        # sum counts
        summed = _sum_counts_by_groups(Xct[ok_cells, :], new_codes[ok_cells], keep.sum())

        # to dense dataframe (DESeq2 prefers dense ints)
        if sp.issparse(summed):
            summed = summed.toarray()

        counts_df = pd.DataFrame(summed, index=kept_samples, columns=genes)
        counts_df = counts_df.round().astype(int)

        # metadata per sample
        meta_cols = [sample_col]
        for col in ["cnv_leiden"]:
            if col in sub.obs.columns:
                meta_cols.append(col)
        meta_df = sub.obs.loc[ok_cells, meta_cols].copy()
        meta_df[sample_col] = meta_df[sample_col].astype(str)
        meta_df = meta_df.drop_duplicates(subset=[sample_col]).set_index(sample_col).reindex(kept_samples)

        pb_counts[str(ct)] = counts_df
        pb_meta[str(ct)] = meta_df

    return pb_counts, pb_meta


def run_deseq2_cluster2_vs_others(
    counts_df,
    meta_df,
    cluster_col="cnv_leiden",
    cluster2_value=2,
    min_total_counts=10,
    min_samples_per_group=2,
):
    """
    counts_df: samples x genes raw counts
    meta_df:   samples x covariates (index aligned with counts_df.index)
    Returns: results_df with log2FoldChange, padj, etc.
    """
    # align
    meta_df = meta_df.loc[counts_df.index].copy()

    # define condition
    cond = np.where(meta_df[cluster_col].astype(int).values == int(cluster2_value), "CNA2", "Other")
    meta = pd.DataFrame({"condition": pd.Categorical(cond)}, index=counts_df.index)

    # group sizes
    vc = meta["condition"].value_counts()
    if (vc.get("CNA2", 0) < min_samples_per_group) or (vc.get("Other", 0) < min_samples_per_group):
        return None

    # filter genes low counts
    keep_genes = (counts_df.sum(axis=0) >= min_total_counts)
    counts_f = counts_df.loc[:, keep_genes].copy()
    if counts_f.shape[1] < 50:
        return None

    dds = DeseqDataSet(
        counts=counts_f,
        metadata=meta,
        design_factors="condition",
        refit_cooks=True,
        n_cpus=1,
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=("condition", "CNA2", "Other"), n_cpus=1)
    stat_res.summary()

    res = stat_res.results_df.copy()
    # harmonize column names
    # pydeseq2 gives: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    return res

# %% --- (D) Pseudobulk par celltype (top) puis DESeq2 cluster2 vs others
# PyDESeq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# remove the current cnv_leiden (clones) by the wes clusters
adata.obs["cnv_leiden"] = np.nan
mapping = wes_cluster.set_index("orig.ident")["cnv_leiden"]

adata.obs["cnv_leiden"] = adata.obs["orig.ident"].map(mapping)

pb_counts, pb_meta = make_pseudobulk_counts(
    adata,
    celltype_col="celltype_level4_curated",
    sample_col="orig.ident",
    target_celltypes=top_types,
    layer="unnormalized_counts",
    min_cells_per_sample=30,
)

all_deseq = {}
all_gene_lists = {}

for ct in top_types:
    if ct not in pb_counts:
        print(f"[SKIP] {ct}: pas assez d'échantillons après filtrage.")
        continue

    counts_df = pb_counts[ct]
    meta_df = pb_meta[ct]

    if "cnv_leiden" not in meta_df.columns:
        raise ValueError("cnv_leiden absent de meta_df. Il faut l'ajouter (merge WES -> adata.obs)")

    res = run_deseq2_cluster2_vs_others(
        counts_df,
        meta_df,
        cluster_col="cnv_leiden",
        cluster2_value=2,
        min_total_counts=10,
        min_samples_per_group=2,
    )

    if res is None:
        print(f"[SKIP] {ct}: DESeq2 non lancé (pas assez de gènes/samples).")
        continue

    # store
    res = res.sort_values("padj")
    all_deseq[ct] = res

    # gene lists (significant)
    sig = res[(res["padj"] < 0.05) & res["log2FoldChange"].notna()].copy()
    up = sig[sig["log2FoldChange"] > 0].index.astype(str).tolist()
    down = sig[sig["log2FoldChange"] < 0].index.astype(str).tolist()
    all_gene_lists[ct] = {"up": up, "down": down}

    print(f"[OK] {ct}: sig={sig.shape[0]} (up={len(up)}, down={len(down)})")

# %%
# --- Charger le DEA bulkRNA
bulk = pd.read_csv("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/RNAseq/wes_cluster_rna/DEG_cluster2_vs_cluster0_1_3.csv")

# nom du gène
bulk = bulk.rename(columns={"Unnamed: 0": "gene"})
bulk["gene"] = bulk["gene"].astype(str)

# filtrage significance
padj_thr = 0.05
lfc_thr = 0.25

bulk_sig = bulk[(bulk["padj"] < padj_thr) & bulk["log2FoldChange"].notna()].copy()

bulk_up = set(bulk_sig[bulk_sig["log2FoldChange"] > lfc_thr]["gene"])
bulk_down = set(bulk_sig[bulk_sig["log2FoldChange"] < -lfc_thr]["gene"])

print("Bulk up:", len(bulk_up))
print("Bulk down:", len(bulk_down))

# %%
import numpy as np
import pandas as pd

bulk_up_genes = sorted(set(bulk_sig.loc[bulk_sig["log2FoldChange"] > lfc_thr, "gene"]))
bulk_down_genes = sorted(set(bulk_sig.loc[bulk_sig["log2FoldChange"] < -lfc_thr, "gene"]))
bulk_genes = sorted(set(bulk_up_genes + bulk_down_genes))

print(f"Bulk sig genes: {len(bulk_genes)} (up={len(bulk_up_genes)}, down={len(bulk_down_genes)})")


# --- 2) Construire une table longue: (celltype, gene, sc_log2FC, sc_padj, sc_pvalue)
rows = []
for ct, res in all_deseq.items():
    # res: index = genes, colonnes: log2FoldChange, padj, pvalue, ...
    # on ne garde que les gènes bulk
    sub = res.loc[res.index.intersection(bulk_genes), ["log2FoldChange", "padj", "pvalue"]].copy()
    sub = sub.reset_index().rename(columns={"index": "gene", "log2FoldChange": "sc_log2FC", "padj": "sc_padj", "pvalue": "sc_pvalue"})
    sub["celltype"] = ct
    rows.append(sub)

sc_on_bulk_long = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
print("Rows in sc_on_bulk_long:", sc_on_bulk_long.shape)

# Sauvegarde table longue
sc_on_bulk_long.to_csv("sc_pseudobulk_stats_on_bulkgenes_long.csv", index=False)


# --- 3) Matrices celltype x gene (log2FC, padj, -log10(padj))
lfc_mat = sc_on_bulk_long.pivot(index="celltype", columns="gene", values="sc_log2FC")

padj_mat = sc_on_bulk_long.pivot(index="celltype", columns="gene", values="sc_padj")
neglog10padj_mat = -np.log10(padj_mat)

# option: clip inf si padj=0 (rare)
neglog10padj_mat = neglog10padj_mat.replace([np.inf, -np.inf], np.nan)

# Sauvegarde matrices
lfc_mat.to_csv("sc_pseudobulk_log2FC_celltype_x_bulkgenes.csv")
padj_mat.to_csv("sc_pseudobulk_padj_celltype_x_bulkgenes.csv")
neglog10padj_mat.to_csv("sc_pseudobulk_neglog10padj_celltype_x_bulkgenes.csv")


# --- 4) Bonus: score simple de "concordance bulk" par celltype
# (combien de gènes bulk up qui sont aussi up en sc et significatifs, idem down)
bulk_up = set(bulk_up_genes)
bulk_down = set(bulk_down_genes)

summary = []
for ct, res in all_deseq.items():
    # res index = genes
    sub = res.loc[res.index.intersection(bulk_genes)].copy()
    sub = sub[(sub["padj"] < padj_thr) & sub["log2FoldChange"].notna()]
    sc_up = set(sub.index[sub["log2FoldChange"] > 0])
    sc_down = set(sub.index[sub["log2FoldChange"] < 0])

    summary.append({
        "celltype": ct,
        "bulkUP_scUP_sig": len(bulk_up & sc_up),
        "bulkDOWN_scDOWN_sig": len(bulk_down & sc_down),
        "bulkUP_scDOWN_sig": len(bulk_up & sc_down),
        "bulkDOWN_scUP_sig": len(bulk_down & sc_up),
        "n_sc_sig_on_bulkgenes": sub.shape[0],
    })

summary_df = pd.DataFrame(summary).sort_values(["bulkUP_scUP_sig", "bulkDOWN_scDOWN_sig"], ascending=False)
#summary_df.to_csv("sc_on_bulkgenes_concordance_summary.csv", index=False)
print(summary_df)

# %%
sig_plot = neglog10padj_mat.copy()

# limiter pour lisibilité
sig_plot = sig_plot.loc[:, sig_plot.max(axis=0) > 2]

plt.figure(figsize=(16,6))

sns.heatmap(
    sig_plot,
    cmap="viridis",
    cbar_kws={"label": "-log10(padj)"}
)

plt.ylabel("Cell type")
plt.xlabel("Bulk significant genes")
plt.title("Significance of bulk genes in single-cell pseudobulk")

plt.xticks(rotation=90)
plt.tight_layout()

plt.savefig("heatmap_sc_pseudobulk_significance.png", dpi=300)
plt.show()


# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# PARAMETRES
# =========================
TOP_N_BULK = 100          # 30-60 en général
MAX_SIG = 10             # cap -log10(padj) pour la taille des points
LFC_CAP = 2              # cap log2FC pour la couleur

# =========================
# 1) Choisir les gènes bulk à afficher
# =========================
bulk2 = bulk.copy()
if "gene" not in bulk2.columns:
    bulk2 = bulk2.rename(columns={"Unnamed: 0": "gene"})
bulk2["gene"] = bulk2["gene"].astype(str)

bulk2 = bulk2[bulk2["padj"].notna()].sort_values("padj")
top_genes = bulk2.head(TOP_N_BULK)["gene"].tolist()

# garder seulement les gènes présents dans scRNA
genes_sc = set(sc_on_bulk_long["gene"])

top_genes = [g for g in top_genes if g in genes_sc]

# =========================
# 2) Préparer table sc pour plot
# =========================
df = sc_on_bulk_long[sc_on_bulk_long["gene"].isin(top_genes)].copy()

# enlever gènes sans aucun signal
gene_signal = (
    df.groupby("gene")["sc_log2FC"]
    .apply(lambda x: x.abs().max())
)

genes_keep = gene_signal[gene_signal > 0.1].index

df = df[df["gene"].isin(genes_keep)]

df["sig"] = -np.log10(df["sc_padj"])
df["sig"] = df["sig"].replace([np.inf, -np.inf], np.nan).fillna(0).clip(0, MAX_SIG)

df["lfc"] = df["sc_log2FC"].clip(-LFC_CAP, LFC_CAP)

# ordre celltypes (mettre tumoraux en bas -> lisible)
celltype_order = sorted(df["celltype"].unique())
# ordre gènes: dans l'ordre bulk (padj croissant)
gene_order = top_genes

# encoder positions
x_map = {g: i for i, g in enumerate(gene_order)}
y_map = {c: i for i, c in enumerate(celltype_order)}

df["x"] = df["gene"].map(x_map)
df["y"] = df["celltype"].map(y_map)

# =========================
# 3) Plot
# =========================
plt.figure(figsize=(max(10, 0.35*len(gene_order)), max(5, 0.35*len(celltype_order))))

sc = plt.scatter(
    df["x"],
    df["y"],
    s=10 + df["sig"] * 25,          # taille ~ significativité
    c=df["lfc"],                    # couleur ~ log2FC
    cmap="RdBu_r",
    vmin=-LFC_CAP,
    vmax=LFC_CAP,
    linewidths=0.2,
    edgecolors="black",
)

cb = plt.colorbar(sc)
cb.set_label("sc pseudobulk log2FC (cluster2 vs others)")

plt.yticks(range(len(celltype_order)), celltype_order)
plt.xticks(range(len(gene_order)), gene_order, rotation=90)

plt.xlabel("Top bulk DE genes (ranked by bulk padj)")
plt.ylabel("Cell types")
plt.title("Which cell types express bulk DE genes (sc pseudobulk dotplot)")

plt.tight_layout()
plt.savefig("dotplot_top_bulk_genes.png", dpi=300)
plt.show()

# %% -- specific cell type dotplot
import numpy as np
import matplotlib.pyplot as plt

TARGET_CT = "T_reg"
TOP_N_CT = 30

df_ct = sc_on_bulk_long[sc_on_bulk_long["celltype"] == TARGET_CT].copy()
df_ct["sig"] = -np.log10(df_ct["sc_padj"])
df_ct["sig"] = df_ct["sig"].replace([np.inf, -np.inf], np.nan).fillna(0)

# garder seulement les gènes qui sont bulk-significatifs (si tu veux)
bulk_sig_genes = set(bulk2[bulk2["padj"] < 0.05]["gene"])
df_ct = df_ct[df_ct["gene"].isin(bulk_sig_genes)].copy()

# top gènes dans ce celltype (par padj sc)
df_ct = df_ct.sort_values("sc_padj").head(TOP_N_CT)

genes_ct = df_ct["gene"].tolist()

# préparer df multi-celltypes sur ces gènes
df = sc_on_bulk_long[sc_on_bulk_long["gene"].isin(genes_ct)].copy()
df["sig"] = -np.log10(df["sc_padj"])
df["sig"] = df["sig"].replace([np.inf, -np.inf], np.nan).fillna(0).clip(0, 10)
df["lfc"] = df["sc_log2FC"].clip(-2, 2)

celltype_order = sorted(df["celltype"].unique())
gene_order = genes_ct

x_map = {g: i for i, g in enumerate(gene_order)}
y_map = {c: i for i, c in enumerate(celltype_order)}
df["x"] = df["gene"].map(x_map)
df["y"] = df["celltype"].map(y_map)

plt.figure(figsize=(max(10, 0.4*len(gene_order)), max(5, 0.35*len(celltype_order))))

sc = plt.scatter(
    df["x"], df["y"],
    s=10 + df["sig"] * 35,
    c=df["lfc"],
    cmap="RdBu_r",
    vmin=-2, vmax=2,
    linewidths=0.2,
    edgecolors="black",
)

cb = plt.colorbar(sc)
cb.set_label("sc pseudobulk log2FC (cluster2 vs others)")

plt.yticks(range(len(celltype_order)), celltype_order)
plt.xticks(range(len(gene_order)), gene_order, rotation=90)

plt.xlabel(f"Bulk DE genes most significant in {TARGET_CT} (by sc padj)")
plt.ylabel("Cell types")
plt.title(f"Bulk signature genes highlighted by {TARGET_CT}")

plt.tight_layout()
plt.savefig(f"dotplot_bulk_genes_focus_{TARGET_CT}.png", dpi=300)
plt.show()



# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# PARAMETRES (à ajuster)
# ============================================================
TOP_BULK = 30                 # nb de gènes bulk dans le bloc "global"
TOP_PER_CELLTYPE = 15         # nb de gènes spécifiques par celltype
SC_PADJ_THR = 0.05            # seuil sc pour sélectionner gènes spécifiques
MIN_ABS_LFC_SC = 0.25         # optionnel: effet minimal sc pour spécificité

MAX_SIG = 10                  # cap -log10(padj) pour la taille des points
LFC_CAP = 2                   # cap log2FC pour la couleur

# Cell types à utiliser pour le module "spécifique"
CELLTYPES_OF_INTEREST = ["T_cell_CD8", "T_cell_CD4", "Microglia_homeostatic"]


# ============================================================
# 1) Préparer bulk + sc tables (noms colonnes)
# ============================================================
bulk2 = bulk.copy()
if "gene" not in bulk2.columns:
    if "Unnamed: 0" in bulk2.columns:
        bulk2 = bulk2.rename(columns={"Unnamed: 0": "gene"})
    else:
        raise ValueError("bulk doit contenir une colonne 'gene' ou 'Unnamed: 0'.")

bulk2["gene"] = bulk2["gene"].astype(str)
bulk2 = bulk2[bulk2["padj"].notna()].sort_values("padj")

sc = sc_on_bulk_long.copy()
sc["gene"] = sc["gene"].astype(str)

genes_sc = set(sc["gene"])


# ============================================================
# 2) Bloc global : top gènes bulk (présents en sc)
# ============================================================
bulk_top = bulk2.head(200)["gene"].tolist()  # on prend large puis on filtre
bulk_top = [g for g in bulk_top if g in genes_sc]
bulk_global_genes = bulk_top[:TOP_BULK]


# ============================================================
# 3) Bloc spécifique : top gènes sc par cell type (sur la liste bulk)
#    Ici : on prend dans chaque celltype les gènes bulk les + significatifs en sc
# ============================================================
specific_genes = []
for ct in CELLTYPES_OF_INTEREST:
    if ct not in sc["celltype"].unique():
        continue

    sub = sc[sc["celltype"] == ct].copy()
    sub = sub[sub["sc_padj"].notna()]
    sub = sub[sub["sc_padj"] < SC_PADJ_THR]

    # option: exiger un minimum d'effet
    sub = sub[sub["sc_log2FC"].abs() >= MIN_ABS_LFC_SC]

    sub = sub.sort_values("sc_padj").head(TOP_PER_CELLTYPE)
    specific_genes += sub["gene"].tolist()

# unique en gardant l'ordre
specific_genes = list(dict.fromkeys(specific_genes))


# ============================================================
# 4) Construire la liste finale (global + overlap + specific)
#    - global_only
#    - overlap (gènes présents à la fois dans global ET spécifiques)
#    - specific_only
# ============================================================
global_set = set(bulk_global_genes)
spec_set = set(specific_genes)

overlap = [g for g in bulk_global_genes if g in spec_set]  # overlap ordonné comme global
global_only = [g for g in bulk_global_genes if g not in spec_set]
specific_only = [g for g in specific_genes if g not in global_set]

gene_order = global_only + overlap + specific_only

# si jamais aucun spécifique trouvé (rare), on garde juste global
if len(gene_order) == 0:
    raise ValueError("Aucun gène sélectionné (vérifie TOP_BULK / SC_PADJ_THR / MIN_ABS_LFC_SC).")


# ============================================================
# 5) Sous-table sc pour plot + metrics (taille/couleur)
# ============================================================
df = sc[sc["gene"].isin(gene_order)].copy()

df["sig"] = -np.log10(df["sc_padj"])
df["sig"] = df["sig"].replace([np.inf, -np.inf], np.nan).fillna(0).clip(0, MAX_SIG)

df["lfc"] = df["sc_log2FC"].clip(-LFC_CAP, LFC_CAP)


# ============================================================
# 6) Ordre celltypes (on met ceux d'intérêt en haut + le reste)
# ============================================================
all_celltypes = list(df["celltype"].unique())

# ordre: d'abord ceux d'intérêt (dans l'ordre donné), puis les autres alphabétiques
celltype_order = [ct for ct in CELLTYPES_OF_INTEREST if ct in all_celltypes]
celltype_order += sorted([ct for ct in all_celltypes if ct not in celltype_order])


# ============================================================
# 7) Encoder positions + enlever colonnes vides (sécurité)
# ============================================================
x_map = {g: i for i, g in enumerate(gene_order)}
y_map = {c: i for i, c in enumerate(celltype_order)}

df["x"] = df["gene"].map(x_map)
df["y"] = df["celltype"].map(y_map)

# enlever lignes sans mapping
df = df.dropna(subset=["x", "y"]).copy()
df["x"] = df["x"].astype(int)
df["y"] = df["y"].astype(int)

# enlever gènes qui n'ont aucun point (au cas où)
genes_with_points = df["gene"].unique().tolist()
gene_order = [g for g in gene_order if g in genes_with_points]
x_map = {g: i for i, g in enumerate(gene_order)}
df["x"] = df["gene"].map(x_map).astype(int)

# %%
# ============================================================
# 8) Plot
# ============================================================
fig_w = max(12, 0.32 * len(gene_order))
fig_h = max(5, 0.35 * len(celltype_order))

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# --- scatter
scat = ax.scatter(
    df["x"], df["y"],
    s=8 + df["sig"] * 28,
    c=df["lfc"],
    cmap="RdBu_r",
    vmin=-LFC_CAP, vmax=LFC_CAP,
    linewidths=0.2,
    edgecolors="black",
)

# --- colorbar (avec pad pour laisser la place à la size legend)
cb = plt.colorbar(scat, ax=ax, fraction=0.04, pad=0.06)
cb.set_label("sc pseudobulk log2FC (cluster2 vs others)")

# --- axes ticks/labels
ax.set_yticks(range(len(celltype_order)))
ax.set_yticklabels(celltype_order)

ax.set_xticks(range(len(gene_order)))
ax.set_xticklabels(gene_order, rotation=90)

ax.set_ylabel("")
ax.set_xlabel("")  # on enlève le long xlabel pour éviter les conflits
ax.set_title("Representation of bulk DE genes in sc pseudobulk")

# --- labels de blocs SOUS les ticks (sans chevauchement)
# réserve une marge en bas
fig.subplots_adjust(bottom=0.40)

n = len(gene_order)

def xfrac(idx):
    return idx / max(n - 1, 1)

def center_frac(i0, i1):
    return xfrac((i0 + i1) / 2)

# --- size legend (significance) placée sous la colorbar
size_values = [1, 2, 5]
size_labels = ["p = 0.1", "p = 0.01", "p = 1e-5"]

handles = []
for s, lab in zip(size_values, size_labels):
    handles.append(
        plt.scatter([], [], s=8 + s * 28,
                    color="grey", edgecolors="black", label=lab)
    )

size_legend = fig.legend(
    handles=handles,
    title="Significance\n(-log10 padj)",
    frameon=False,
    loc="center left",
    bbox_to_anchor=(0.92, 0.15)   # position dans la FIGURE
)
ax.add_artist(size_legend)

# --- layout final (évite de recouper les textes)
plt.tight_layout(rect=[0, 0.08, 1, 1])

# --- save/show
plt.savefig(
    "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/suppl_figures3/dotplot_global_plus_celltype_specific.png",
    dpi=800
)
plt.show()

print("Global genes:", len(bulk_global_genes))
print("Specific genes:", len(specific_genes))
print("Overlap:", len(overlap))
print("Final plotted genes:", len(gene_order))
# %%
