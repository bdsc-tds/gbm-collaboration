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
import pyranges as pr

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

from gbmcollaboration.viz.figure_helpers_py import savefig, set_mpl_rcparams

set_mpl_rcparams()

# %% -- Paths
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/")

SC_ADATA_PATH = RESULTS_DIR / "scRNAseq/infercnv/infercnv_adata.h5ad"
WES_MATRIX_PATH = RESULTS_DIR / "WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
WES_BINS_ANNOT_PATH = RESULTS_DIR / "scRNAseq/infercnv/genomic_windows.csv"

SAMPLE_COL = "orig.ident"

EXCLUDE_CHROMS = {"chrX", "chrY"}

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
    print(f"Missing WES profiles for samples: {sorted(missing.unique().tolist())}")

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

# %% -- Exclude chrX and Y (and anything in EXCLUDE_CHROMS)
before = wes_adata.n_vars
mask = ~wes_adata.var["chromosome"].astype(str).isin(EXCLUDE_CHROMS)
wes_adata = wes_adata[:, mask].copy()
after = wes_adata.n_vars
print(f"Excluded chromosomes {sorted(EXCLUDE_CHROMS)}: dropped {before - after} bins")

# Optional sanity check
assert not wes_adata.var["chromosome"].astype(str).isin(EXCLUDE_CHROMS).any()

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
cnv.tl.leiden(wes_adata, resolution=0.4)

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
    dendrogram=True,
    save=True
    )

# %% -- Count how many samples found in each cluster
print(wes_adata.obs["cnv_leiden"].value_counts())

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

# %% -- Calculate CNV scores per WES CNA cluster to report in manuscript
wes_adata.obs.groupby("cnv_leiden").agg({"cnv_score": "mean"})

# %% -- Pursue exploration to highlight the difference between cluster 2 and the other clusters
# ============================================================
# FULL, SINGLE BLOCK (Tukey-only bin-level signature vs cluster 2)
# - Computes Tukey HSD per bin
# - Δ = mean(cluster2) - mean(cluster k) for k in {0,1,3}
# - Global FDR across all (bin × comparison)
# - Keeps ONLY bins significant in ALL 3 comparisons (2vs0, 2vs1, 2vs3)
# - Optional: require consistent direction across the 3 deltas
# - Optional: group consecutive bins into segments (after filtering)
# - Plot with chromosome-only ticks (centered), no bin labels
# - Exports tables for downstream gene annotation
# ============================================================

# -- parms
CL2 = "2"
COMPARISONS = ["0", "1", "3"]
COMP_ORDER = [f"2 vs {k}" for k in COMPARISONS]

Q_FDR = 0.05

REQUIRE_CONSISTENT_DIRECTION = True   # True = deltas all >0 or all <0 across 2vs0/1/3

GROUP_TO_SEGMENTS = True              # True = merge consecutive bins into segments
GAP_BP = 2_000_000                    # max gap between consecutive bins to be in same segment

# Plot typography
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# %% -- inputs
if "cluster_labels" not in globals():
    cluster_labels = wes_adata.obs["cnv_leiden"].astype(str)
cluster_labels = cluster_labels.astype(str)

cnv_df = pd.DataFrame(
    wes_adata.obsm["X_cnv"],
    index=wes_adata.obs_names,
    columns=wes_adata.var_names
)

var = wes_adata.var.copy()
var["chromosome"] = var["chromosome"].astype(str)
var["start"] = var["start"].astype(int)
var["end"] = var["end"].astype(int)

cluster_labels = cluster_labels.reindex(cnv_df.index)
assert cluster_labels.notna().all(), "cluster_labels and CNV matrix indices must align."

for c in [CL2] + COMPARISONS:
    assert c in set(cluster_labels.unique()), f"Missing cluster label '{c}'."

print("Cluster sizes:")
print(cluster_labels.value_counts().sort_index())

# orders genes by genomic position
def chr_rank(c):
    c = str(c).replace("chr", "")
    if c.isdigit(): return int(c)
    if c == "X": return 23
    if c == "Y": return 24
    return 99

var2 = var.loc[cnv_df.columns, ["chromosome","start","end"]].copy()
var2["chr_rank"] = var2["chromosome"].map(chr_rank)
bin_order = var2.sort_values(["chr_rank","start","end"]).index.astype(str)

# -- long table (bin × comparison)
# delta = mean(C2) - mean(k)
rows = []
g = cluster_labels.values

for b in bin_order:
    y = cnv_df[b].values
    tk = pairwise_tukeyhsd(endog=y, groups=g, alpha=0.05)
    tk_df = pd.DataFrame(tk.summary().data[1:], columns=tk.summary().data[0])
    # statsmodels Tukey: meandiff = mean(group2) - mean(group1)

    for k in COMPARISONS:
        comp = f"2 vs {k}"

        r_k2 = tk_df[(tk_df["group1"] == k) & (tk_df["group2"] == CL2)]  # meandiff = mean(2) - mean(k)
        r_2k = tk_df[(tk_df["group1"] == CL2) & (tk_df["group2"] == k)]  # meandiff = mean(k) - mean(2)

        if not r_k2.empty:
            delta = float(r_k2["meandiff"].iloc[0])
            p_adj = float(r_k2["p-adj"].iloc[0])
        elif not r_2k.empty:
            delta = -float(r_2k["meandiff"].iloc[0])
            p_adj = float(r_2k["p-adj"].iloc[0])
        else:
            continue

        rows.append({
            "bin": str(b),
            "comp": comp,
            "delta": delta,
            "abs_delta": abs(delta),
            "p_tukey_adj": p_adj
        })

tukey_long = pd.DataFrame(rows)
tukey_long = tukey_long.join(var2[["chromosome","start","end","chr_rank"]], on="bin")

# Global FDR across all bin×comparison points
_, qvals, _, _ = multipletests(tukey_long["p_tukey_adj"].values, method="fdr_bh")
tukey_long["q_fdr"] = qvals

# -- keep only bins significant in ALL 3 comparisons
pivot_q = tukey_long.pivot_table(index="bin", columns="comp", values="q_fdr", aggfunc="min").reindex(columns=COMP_ORDER)
mask_all = pivot_q.notna().all(axis=1) & (pivot_q < Q_FDR).all(axis=1)
bins_strict = pivot_q.index[mask_all].astype(str)

print(f"\nBins significant in ALL 3 comparisons (FDR<{Q_FDR}): {len(bins_strict)}")
if len(bins_strict) == 0:
    raise ValueError("No bins pass the 'significant in all 3 comparisons' filter. Try Q_FDR=0.1.")

# consistent direction across all three comparisons
if REQUIRE_CONSISTENT_DIRECTION:
    pivot_delta = tukey_long.pivot_table(index="bin", columns="comp", values="delta", aggfunc="mean").reindex(columns=COMP_ORDER)
    dir_ok = ((pivot_delta > 0).all(axis=1)) | ((pivot_delta < 0).all(axis=1))
    bins_strict = bins_strict[dir_ok.loc[bins_strict].values]
    print(f"Bins after requiring consistent direction: {len(bins_strict)}")
    if len(bins_strict) == 0:
        raise ValueError("No bins left after consistent-direction filter. Set REQUIRE_CONSISTENT_DIRECTION=False.")

# Keep strict points only
P = tukey_long[tukey_long["bin"].isin(bins_strict)].copy()

# Ordered bins metadata for plotting / segmentation
sig_bins = (
    P[["bin","chromosome","start","end","chr_rank"]]
    .drop_duplicates()
    .sort_values(["chr_rank","start","end"])
    .reset_index(drop=True)
)

# optional: group consecutive bins into segments (after strict filtering)
if GROUP_TO_SEGMENTS:
    bdf = sig_bins.copy()

    seg_id = []
    cur = -1
    prev_chr = None
    prev_end = None
    for _, r in bdf.iterrows():
        chr_ = r["chromosome"]
        st = int(r["start"])
        en = int(r["end"])
        if (prev_chr is None) or (chr_ != prev_chr) or (prev_end is None) or (st - prev_end > GAP_BP):
            cur += 1
        seg_id.append(cur)
        prev_chr, prev_end = chr_, en

    bdf["segment"] = seg_id
    bin_to_seg = bdf.set_index("bin")["segment"]

    # segment coords (includes chr_rank)
    seg_coords = (
        bdf.groupby("segment", as_index=False)
           .agg(
               chromosome=("chromosome","first"),
               chr_rank=("chr_rank","first"),
               start=("start","min"),
               end=("end","max"),
               n_bins=("bin","count"),
               bins=("bin", lambda x: ",".join(x.astype(str)))
           )
           .sort_values(["chr_rank","start","end"])
           .reset_index(drop=True)
    )

    # aggregate per segment×comparison
    P["segment"] = P["bin"].map(bin_to_seg)
    P = (P.groupby(["segment","comp"], as_index=False)
           .agg(
               delta=("delta","mean"),
               abs_delta=("abs_delta","max"),
               q_fdr=("q_fdr","min")
           )
         ).merge(seg_coords, on="segment", how="left")

    # assign x positions from seg_coords order
    seg_coords["x"] = np.arange(len(seg_coords))
    P = P.merge(seg_coords[["segment","x"]], on="segment", how="left")

    # chromosome array along x for separators/ticks
    chroms_x = seg_coords["chromosome"].values
    n_x = len(chroms_x)

else:
    # bin-level x positions from genome order
    sig_bins["x"] = np.arange(len(sig_bins))
    P = P.merge(sig_bins[["bin","x"]], on="bin", how="left")

    chroms_x = sig_bins["chromosome"].values
    n_x = len(chroms_x)

# Ensure y exists
y_map = {c:i for i,c in enumerate(COMP_ORDER)}
P["y"] = P["comp"].map(y_map)

# AESTHETICS: sizes & color scale
absd = P["abs_delta"].values
absd95 = max(np.quantile(absd, 0.95), 1e-6)
sizes = 18 + 290 * np.clip(absd / absd95, 0, 1)

d = P["delta"].values
d95 = max(np.quantile(np.abs(d), 0.95), 1e-6)
vmin, vmax = -d95, d95

# Chromosome boundaries and center ticks
chr_change = np.where(chroms_x[:-1] != chroms_x[1:])[0]
chr_centers = {}
for chr_ in pd.unique(chroms_x):
    idx = np.where(chroms_x == chr_)[0]
    chr_centers[chr_] = idx.mean()

# %% -- Plotting
fig_w = max(6, 0.18 * n_x)
fig_h = 2.8
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

sc = ax.scatter(
    P["x"], P["y"],
    c=P["delta"], s=sizes,
    cmap="bwr", vmin=vmin, vmax=vmax,
    edgecolors="black", linewidths=0.5, alpha=0.95,
    rasterized=True
)

for pos in chr_change:
    ax.axvline(pos + 0.5, color="black", lw=0.7, alpha=0.6)

ax.set_yticks(range(len(COMP_ORDER)))
ax.set_yticklabels(COMP_ORDER)
ax.set_ylabel("Clusters")

# Chromosome ticks at centers
chr_items = sorted(chr_centers.items(), key=lambda kv: kv[1])
chr_ticks = [v for _, v in chr_items]
chr_labels = [k for k, _ in chr_items]
ax.set_xticks(chr_ticks)
ax.set_xticklabels(chr_labels, rotation=30, ha="right")
ax.set_xlabel("")

ax.set_title(
    "Segment-level comparison between WES CNV clusters")

ax.set_xlim(-0.8, n_x - 0.2)
ax.set_ylim(-0.6, len(COMP_ORDER) - 0.4)
ax.set_axisbelow(True)
ax.grid(axis="y", linewidth=0.3, alpha=0.25)

cbar = fig.colorbar(sc, ax=ax, fraction=0.02, pad=0.02)
cbar.set_label("Δ mean log2ratio")

plt.tight_layout(pad=1)
savefig(
    df=P,
    out_file="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/figure3/fig3c.pdf"
)
plt.show()

# %% -- EXPORT TABLES (for gene overlap)
# Long-form table (one row per locus × comparison)
strict_long = P.copy()

# Summary per locus (bin or segment)
if GROUP_TO_SEGMENTS:
    strict_summary = (
        strict_long.groupby(["segment","chromosome","start","end","n_bins"], as_index=False)
                .agg(
                    chr_rank=("chr_rank","first"),     # <--- add this
                    min_q=("q_fdr","min"),
                    max_abs_delta=("abs_delta","max"),
                    comps=("comp", lambda x: ";".join(sorted(set(x))))
                )
                .sort_values(["chr_rank","start","end"])
                .drop(columns=["chr_rank"])
    )
else:
    strict_summary = (
        tukey_long[tukey_long["bin"].isin(bins_strict)]
        .groupby(["bin","chromosome","start","end"], as_index=False)
        .agg(
            min_q=("q_fdr","min"),
            max_abs_delta=("abs_delta","max"),
            comps=("comp", lambda x: ";".join(sorted(set(x))))
        )
        .merge(sig_bins[["bin","chr_rank"]], on="bin", how="left")
        .sort_values(["chr_rank","start","end"])
        .drop(columns=["chr_rank"])
    )

print("\nStrict loci kept:", strict_summary.shape[0])
print(strict_summary.head(10))

strict_long.to_csv(RESULTS_DIR / "WES/cluster2_strict_loci_tukey_long.csv", index=False)
strict_summary.to_csv(RESULTS_DIR / "WES/cluster2_strict_loci_summary.csv", index=False)

# %% -- Identify the genes on those segments and append them to the table
segments = strict_summary.copy()

# IMPORTANT: pyranges expect 0-based coordinates
# we convert the segments in BED-like: Start0 = start-1, End = end
segments_df = pd.DataFrame({
    "Chromosome": segments["chromosome"].astype(str),
    "Start": (segments["start"].astype(int) - 1).clip(lower=0),
    "End": segments["end"].astype(int),
    "Segment": segments["segment"].astype(str) if "segment" in segments.columns else segments.index.astype(str)
})

gr_segs = pr.PyRanges(segments_df)

# -- read gtf
gtf_path = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/infercnv/gencode.v32.annotation.gtf"  # mets ton chemin
gtf = pr.read_gtf(gtf_path)

# -- keep only genes (no exon/transcript)
genes = gtf[gtf.Feature == "gene"]

# -- keep useful columns
genes_df = genes.df[["Chromosome","Start","End","gene_name","gene_id","Strand"]].drop_duplicates()

gr_genes = pr.PyRanges(genes_df)

# -- intersection segments x genes
hits = gr_segs.join(gr_genes)

hits_df = hits.df.sort_values(["Segment","Chromosome","Start"])
hits_df.head()

# -- list of genes per segment
genes_by_segment = (hits_df.groupby("Segment")["gene_name"]
                    .apply(lambda x: sorted(set(x)))
                    .reset_index(name="genes"))

genes_by_segment["n_genes"] = genes_by_segment["genes"].apply(len)
genes_by_segment.head()

# %% -- Save the gene list inside the previous dataframe
genes_by_segment.rename(columns={"Segment": "segment"}, inplace=True)
genes_by_segment["segment"] = genes_by_segment["segment"].astype(int)
strict_summary["segment"] = strict_summary["segment"].astype(int)

merged_df = pd.merge(strict_summary, genes_by_segment, on="segment")

# %% -- Save the table for the supplementary results
merged_df.to_csv(RESULTS_DIR / "WES/segment_cluster_genes.csv", index=False)
