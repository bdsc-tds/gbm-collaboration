"""Script to run infercnv on WES CNV matrix."""

# To run successfully this script, you need to have pre-computed the
# CNV WES matrix of all samples across centers, and this should have been
# saved in:
# "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
# to do so, run PART 1 of the script "generate_cnv_wes_matrix.py"

# %% -- Imports
import pandas as pd
import numpy as np
from pathlib import Path
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import re
import infercnvpy as cnv

# %% -- Paths
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/")

SC_ADATA_PATH = RESULTS_DIR / "scRNAseq/infercnv/infercnv_adata.h5ad"
WES_MATRIX_PATH = RESULTS_DIR / "WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
WES_BINS_ANNOT_PATH = RESULTS_DIR / "scRNAseq/infercnv/genomic_windows.csv"

SAMPLE_COL = "orig.ident"

# %% -- Load scRNA AnnData
single_cell_obj = ad.read_h5ad(SC_ADATA_PATH)

# %% -- Create pseudobulk expression AnnData (samples x genes)
# (This is not used for WES plotting, but useful for sample IDs + metadata)
X = single_cell_obj.X.toarray() if not isinstance(single_cell_obj.X, np.ndarray) else single_cell_obj.X
df_expr = pd.DataFrame(X, index=single_cell_obj.obs[SAMPLE_COL].astype(str))

# sum across cells per sample -> samples x genes
pseudobulk_matrix = df_expr.groupby(df_expr.index).sum()

pseudobulk_adata = ad.AnnData(X=pseudobulk_matrix.values)
pseudobulk_adata.obs[SAMPLE_COL] = pseudobulk_matrix.index.astype(str)
pseudobulk_adata.obs_names = pseudobulk_matrix.index.astype(str)
pseudobulk_adata.var_names = single_cell_obj.var_names
pseudobulk_adata.var = single_cell_obj.var.copy()

# %% -- Load WES CNV matrix (samples x bins)
wes_profile = pd.read_csv(WES_MATRIX_PATH, index_col=0)
wes_profile.index = wes_profile.index.astype(str)

# %% -- Align samples between pseudobulk and WES
pb_ids = pseudobulk_adata.obs[SAMPLE_COL].astype(str)
wes_ids = wes_profile.index.astype(str)

missing = pb_ids[~pb_ids.isin(wes_ids)]
if len(missing) > 0:
    raise ValueError(f"Missing WES profiles for samples: {sorted(missing.unique().tolist())}")

# %% -- Run if WES profiles are missing
keep_mask = pb_ids.isin(wes_ids).values
print(f"Dropping {(~keep_mask).sum()} pseudobulk sample(s): {pb_ids[~keep_mask].unique().tolist()}")

pseudobulk_adata = pseudobulk_adata[keep_mask].copy()
pseudobulk_adata.obs[SAMPLE_COL] = pseudobulk_adata.obs[SAMPLE_COL].astype(str)

# reorder WES to match pseudobulk sample order
wes_sub = wes_profile.loc[pseudobulk_adata.obs[SAMPLE_COL]].copy()

# %% -- Load bin annotation (bin, chromosome, start, end)
bins_annot = pd.read_csv(WES_BINS_ANNOT_PATH, index_col=0)

# Enforce types
bins_annot["bin"] = bins_annot["bin"].astype(int)
bins_annot["chromosome"] = bins_annot["chromosome"].astype(str)
bins_annot["start"] = bins_annot["start"].astype(int)
bins_annot["end"] = bins_annot["end"].astype(int)
bins_annot = bins_annot.set_index("bin")

# %% -- Map WES columns -> bin IDs (preserves the WES column order)
def col_to_bin_id(col) -> int:
    """
    Parse a bin ID from a column name.
    Works for '123', 123, 'bin_123', 'Bin123', etc.
    If your columns are exactly integers/strings of integers, this is safe.
    """
    s = str(col)
    m = re.search(r"(\d+)", s)
    if not m:
        raise ValueError(f"Could not parse a numeric bin id from WES column name: {col}")
    return int(m.group(1))

wes_bin_ids = [col_to_bin_id(c) for c in wes_sub.columns]

missing_bins = sorted(set(wes_bin_ids) - set(bins_annot.index))
if missing_bins:
    raise ValueError(f"{len(missing_bins)} WES bins missing from annotation (example): {missing_bins[:10]}")

# %% -- Build WES-bin AnnData (samples x bins) with proper .var coordinates
wes_adata = ad.AnnData(
    X=wes_sub.values,  # samples x bins
    obs=pseudobulk_adata.obs.copy(),
    var=pd.DataFrame(index=[str(b) for b in wes_bin_ids])  # var_names are bin ids
)

annot_for_cols = bins_annot.loc[wes_bin_ids].copy()
annot_for_cols.index = annot_for_cols.index.astype(str)  # align with wes_adata.var_names
wes_adata.var = wes_adata.var.join(annot_for_cols, how="left")

if wes_adata.var[["chromosome", "start", "end"]].isna().any().any():
    bad = wes_adata.var[wes_adata.var["chromosome"].isna()].index[:10].tolist()
    raise ValueError(f"Some bins have missing annotation. Example bin ids: {bad}")

# %% -- Optional: verify genome ordering and sort only if needed
def chr_rank(c: str) -> int:
    c = str(c).replace("chr", "")
    if c.isdigit():
        return int(c)      # 1..22
    if c == "X":
        return 23
    if c == "Y":
        return 24
    return 99

rank = np.array([chr_rank(c) for c in wes_adata.var["chromosome"].values], dtype=int)
start = wes_adata.var["start"].to_numpy(dtype=int)

# monotonic check: chr_rank non-decreasing; within chromosome, start non-decreasing
ordered = np.all((rank[1:] > rank[:-1]) | ((rank[1:] == rank[:-1]) & (start[1:] >= start[:-1])))
print("WES bins already genome-ordered (chr/start)?", ordered)

if not ordered:
    # stable, explicit sort
    wes_adata.var["chr_rank"] = rank
    wes_adata.var = wes_adata.var.sort_values(["chr_rank", "start", "end"])
    wes_adata = wes_adata[:, wes_adata.var.index].copy()
    wes_adata.var.drop(columns=["chr_rank"], inplace=True)

# %% -- Make infercnvpy treat this like a CNV representation
# infercnvpy uses obsm["X_cnv"] as the CNV matrix when use_rep="cnv"
wes_adata.obsm["X_cnv"] = wes_adata.X.copy()

# Build chr_pos required by infercnvpy plotting (start index of each chromosome block)
chroms = wes_adata.var["chromosome"].astype(str).to_numpy()

chr_pos = {}
for i, c in enumerate(chroms):
    if c not in chr_pos:
        chr_pos[c] = i

# infercnvpy expects this structure
wes_adata.uns["cnv"] = {"chr_pos": chr_pos}

# %% -- Start calculation and plotting (infercnvpy, on WES bins)
cnv.tl.pca(wes_adata)
cnv.pp.neighbors(wes_adata)
cnv.tl.leiden(wes_adata)

# %% -- Dendrogram using WES CNV
sc.tl.dendrogram(
    wes_adata,
    groupby="cnv_leiden",
    use_rep="X_cnv",
    cor_method="pearson",
    linkage_method="single",      # or single or average
    optimal_ordering=True
)

# %% -- Plotting using infercnvpy
cnv.pl.chromosome_heatmap_summary(
    wes_adata,
    groupby="cnv_leiden",
    figsize=(15,5),
    dendrogram=True)

# %% -- Count how many samples found in each cluster
print(wes_adata.obs["cnv_leiden"].value_counts())
# cnv_leiden
# 0    65
# 1    48
# 2    38
# 3    36
# 4    32
# 5    28
# 6    18

# %% -- Save the Leiden cluster
leiden_cluster_baseline = pd.DataFrame(wes_adata.obs["cnv_leiden"])
leiden_cluster_baseline.to_csv(RESULTS_DIR / "WES/infercnv/leiden_cluster_wes_baseline_samples.csv", sep=",", index=True)

# %% -- Calculate CNV score on the WES profile
# Generate UMAP
cnv.tl.umap(wes_adata)
cnv.tl.cnv_score(wes_adata, groupby="cnv_leiden")

# Create composite figure
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(11, 5))
#ax4.axis("off")  # turn off the empty subplot

# Plot UMAPs into specific axes
cnv.pl.umap(
    wes_adata,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False
)
cnv.pl.umap(wes_adata, color="cnv_score", ax=ax2, show=True)
