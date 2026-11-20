"""Fig 3b, Supp Fig 3a, Supp Fig 3c: WES copy-number clusters.

Produces:
- Fig 3b: mean WES log2 copy ratio per genomic bin of each CNV Leiden cluster
  (`chromosome_heatmap_summary`, clusters ordered by the dendrogram).
- Supp Fig 3a: UMAP of the WES CNV profiles coloured by cluster and by CNV score.
- The CNV scores of the Results ("0.196 vs 0.170, 0.172, 0.183").
- Supp Fig 3c: bins specific to cluster 2 (Tukey HSD vs each other cluster, BH FDR over
  bins x comparisons, merged into segments), and the genes of each segment
  (`segment_cluster_genes`).
- With RUN_CLUSTERING = True: the Leiden clusters (frozen as
  results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv, read by the other Fig 3 files
  and the data release).

Origin: analysis/infercnv/infercnv_wes_analyses.py. Inputs (frozen): the WES CNV matrix
(pipelines/generate_cnv_wes_matrix.py), the genomic bins and the sample IDs of the
single-cell inferCNV object (pipelines/infercnv_sc_analyses.py).

Comparison with the published outputs (printed by the check cells):
- The CNV matrix has no column for chr21 (bin 1274, the only chr21 bin) nor chrY and keeps
  chrX, which is dropped here: Fig 3b shows chr1-22 without chr21, as published.
- The published Supp Fig 3c table and CNV scores (0.196 for cluster 2) were computed with
  sample ER_G_673a in cluster 2, the label a re-run of the Leiden step gives
  (RUN_CLUSTERING = True), while the frozen cluster table, used by every other Fig 3 panel,
  has it in cluster 0. With the frozen table: same three segments (chr10, chr19, chr20), same
  bins and genes, delta differing in the third decimal, cluster-2 CNV score 0.197.
- With RUN_CLUSTERING = False (default) the Leiden labels are read from the frozen table
  instead of being recomputed; with True the new table is written to the output folder.
"""

# %% -- Imports
import re
from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

from gbmcollaboration.paper import output

# %% -- Paths
OUT = output.setup("fig3b_suppfig3ac_infercnv_wes_analyses")
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/")

SC_ADATA_PATH = RESULTS_DIR / "scRNAseq/infercnv/infercnv_adata.h5ad"
WES_MATRIX_PATH = RESULTS_DIR / "WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
WES_BINS_ANNOT_PATH = RESULTS_DIR / "scRNAseq/infercnv/genomic_windows.csv"
CNV_CLUSTERS_PATH = RESULTS_DIR / "WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
GTF_PATH = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/infercnv/gencode.v32.annotation.gtf"

SAMPLE_COL = "orig.ident"

EXCLUDE_CHROMS = {"chrX", "chrY"}

# Re-run the Leiden clustering (True) or read the frozen labels (False)
RUN_CLUSTERING = False
LEIDEN_RESOLUTION = 0.4

# Published values, for the checks
PUBLISHED_CNV_SCORE = {"0": 0.170, "1": 0.172, "2": 0.196, "3": 0.183}
PUBLISHED_STRICT_LOCI = RESULTS_DIR / "WES/cluster2_strict_loci_summary.csv"
PUBLISHED_STRICT_LONG = RESULTS_DIR / "WES/cluster2_strict_loci_tukey_long.csv"
PUBLISHED_SEGMENT_GENES = RESULTS_DIR / "WES/segment_cluster_genes.csv"

# %% -- Functions
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


def chr_rank(c: str) -> int:
    c = str(c).replace("chr", "")
    if c.isdigit():
        return int(c)      # 1..22
    if c == "X":
        return 23
    if c == "Y":
        return 24
    return 99


def read_gtf_genes(gtf_path) -> pd.DataFrame:
    """Gene records of a GTF as pyranges.read_gtf gives them (Start = GTF start - 1):
    Chromosome, Start, End, gene_name, gene_id, Strand."""
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                      names=["Chromosome", "Feature", "Start", "End", "Strand", "attributes"],
                      dtype={"Chromosome": str})
    genes = gtf[gtf["Feature"] == "gene"].copy()
    del gtf
    genes["Start"] = genes["Start"] - 1
    genes["gene_name"] = genes["attributes"].str.extract(r'gene_name "([^"]+)"', expand=False)
    genes["gene_id"] = genes["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
    return genes[["Chromosome", "Start", "End", "gene_name", "gene_id", "Strand"]].drop_duplicates()


def genes_overlapping_segments(segments_df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    """Join of segments (Chromosome, Start, End, ...) and genes on any overlap of the 0-based
    half-open intervals, as `PyRanges(segments_df).join(PyRanges(genes_df)).df` (gene
    coordinates suffixed `_b`)."""
    hits = segments_df.merge(genes_df, on="Chromosome", suffixes=("", "_b"))
    return hits[(hits["Start_b"] < hits["End"]) & (hits["End_b"] > hits["Start"])]


# %% -- Sample IDs of the scRNA inferCNV object (obs only)
single_cell_obj = ad.read_h5ad(SC_ADATA_PATH, backed="r")
sc_obs = single_cell_obj.obs[[SAMPLE_COL]].copy()
single_cell_obj.file.close()
pb_ids = pd.Series(np.sort(sc_obs[SAMPLE_COL].astype(str).unique()), name=SAMPLE_COL)

# %% -- Load WES CNV matrix (samples x bins)
wes_profile = pd.read_csv(WES_MATRIX_PATH, index_col=0)
wes_profile.index = wes_profile.index.astype(str)

# %% -- Align samples between pseudobulk and WES
wes_ids = wes_profile.index.astype(str)

missing = pb_ids[~pb_ids.isin(wes_ids)]
if len(missing) > 0:
    print(f"Missing WES profiles for samples: {sorted(missing.unique().tolist())}")

# %% -- Run if WES profiles are missing
keep_mask = pb_ids.isin(wes_ids).values
print(f"Dropping {(~keep_mask).sum()} pseudobulk sample(s): {pb_ids[~keep_mask].unique().tolist()}")
pb_ids = pb_ids[keep_mask]

# reorder WES to match pseudobulk sample order
wes_sub = wes_profile.loc[pb_ids].copy()

# %% -- Load bin annotation (bin, chromosome, start, end)
bins_annot = pd.read_csv(WES_BINS_ANNOT_PATH, index_col=0)

# Enforce types
bins_annot["bin"] = bins_annot["bin"].astype(int)
bins_annot["chromosome"] = bins_annot["chromosome"].astype(str)
bins_annot["start"] = bins_annot["start"].astype(int)
bins_annot["end"] = bins_annot["end"].astype(int)
bins_annot = bins_annot.set_index("bin")

absent = bins_annot.index.difference(wes_profile.columns.astype(int))
print("Bins of the genomic windows absent from the WES matrix:",
      {int(b): bins_annot.loc[b, "chromosome"] for b in absent})

# %% -- Map WES columns -> bin IDs (preserves the WES column order)
wes_bin_ids = [col_to_bin_id(c) for c in wes_sub.columns]

missing_bins = sorted(set(wes_bin_ids) - set(bins_annot.index))
if missing_bins:
    raise ValueError(f"{len(missing_bins)} WES bins missing from annotation (example): {missing_bins[:10]}")

# %% -- Build WES-bin AnnData (samples x bins) with proper .var coordinates
wes_adata = ad.AnnData(
    X=wes_sub.values,  # samples x bins
    obs=pd.DataFrame({SAMPLE_COL: pb_ids.to_numpy()}, index=pb_ids.to_numpy()),
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
print(f"{wes_adata.n_obs} samples x {wes_adata.n_vars} bins; chromosomes: {list(chr_pos)}")

# %% -- Start calculation and plotting (infercnvpy, on WES bins)
cnv.tl.pca(wes_adata)
cnv.pp.neighbors(wes_adata)
frozen_clusters = pd.read_csv(CNV_CLUSTERS_PATH).set_index(SAMPLE_COL)["cnv_leiden"].astype(str)
if RUN_CLUSTERING:
    cnv.tl.leiden(wes_adata, resolution=LEIDEN_RESOLUTION)
    new = wes_adata.obs["cnv_leiden"].astype(str)
    diff = new.index[new != frozen_clusters.reindex(new.index)]
    print("Samples with a label different from the frozen table (frozen -> new):",
          {s: f"{frozen_clusters.get(s)} -> {new[s]}" for s in diff})
else:
    assert set(frozen_clusters.index) == set(wes_adata.obs_names), "frozen table and samples differ"
    wes_adata.obs["cnv_leiden"] = pd.Categorical(frozen_clusters.loc[wes_adata.obs_names].to_numpy())

# %% -- Dendrogram using WES CNV
sc.tl.dendrogram(
    wes_adata,
    groupby="cnv_leiden",
    use_rep="X_cnv",
    cor_method="pearson",
    linkage_method="single",      # or single or average
    optimal_ordering=True
)

# %% -- Plotting using infercnvpy (Fig 3b)
cnv.pl.chromosome_heatmap_summary(
    wes_adata,
    groupby="cnv_leiden",
    figsize=(15,5),
    dendrogram=True,
    show=False
    )
fig3b_data = pd.DataFrame(wes_adata.obsm["X_cnv"], columns=wes_adata.var_names).groupby(
    wes_adata.obs["cnv_leiden"].to_numpy()).mean()
output.save(OUT, "heatmap_summary_wes_cnv_leiden", plt.gcf(), fig3b_data, panel="Fig 3b")
plt.show()
plt.close("all")
print("Dendrogram order:", wes_adata.uns["dendrogram_cnv_leiden"]["categories_ordered"])

# %% -- Count how many samples found in each cluster
print(wes_adata.obs["cnv_leiden"].value_counts())

# %% -- Save the Leiden cluster
if RUN_CLUSTERING:
    leiden_cluster_baseline = pd.DataFrame(wes_adata.obs["cnv_leiden"])
    leiden_cluster_baseline.index.name = SAMPLE_COL
    output.save_table(OUT, "leiden_cluster_wes_baseline_samples", leiden_cluster_baseline)

# %% -- Calculate CNV score on the WES profile (Supp Fig 3a)
# Generate UMAP
cnv.tl.umap(wes_adata)
cnv.tl.cnv_score(wes_adata, groupby="cnv_leiden")

# Create composite figure
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(11, 5))

# Plot UMAPs into specific axes
cnv.pl.umap(
    wes_adata,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False
)
cnv.pl.umap(wes_adata, color="cnv_score", ax=ax2, show=False)
output.save(OUT, "umap_wes_cnv", fig, panel="Supp Fig 3a")
plt.show()
plt.close(fig)

# %% -- Calculate CNV scores per WES CNA cluster to report in manuscript
cnv_scores = wes_adata.obs.groupby("cnv_leiden", observed=True).agg({"cnv_score": "mean"})
cnv_scores["published"] = pd.Series(PUBLISHED_CNV_SCORE)
print(cnv_scores.round(4))
output.save_table(OUT, "cnv_score_per_cluster", cnv_scores)

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
TUKEY_RC = {
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}

# %% -- inputs
cluster_labels = wes_adata.obs["cnv_leiden"].astype(str)

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

# %% -- Plotting (Supp Fig 3c)
with plt.rc_context(TUKEY_RC):
    fig_w = max(6, 0.18 * n_x)
    fig_h = 2.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    sc_plot = ax.scatter(
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

    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Δ mean log2ratio")

    plt.tight_layout(pad=1)
    output.save(OUT, "suppfig3c", fig, P, panel="Supp Fig 3c")
    plt.show()
    plt.close(fig)

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

output.save_table(OUT, "cluster2_strict_loci_tukey_long", strict_long, index=False)
output.save_table(OUT, "cluster2_strict_loci_summary", strict_summary, index=False)

published_loci = pd.read_csv(PUBLISHED_STRICT_LOCI)
cmp_cols = ["segment", "chromosome", "start", "end", "n_bins", "min_q", "max_abs_delta", "comps"]
same = strict_summary[cmp_cols].reset_index(drop=True).round(4).equals(published_loci[cmp_cols].round(4))
print(f"Identical to the published {PUBLISHED_STRICT_LOCI.name}: {same}")
print(f"Same segments and bins: "
      f"{strict_summary[['chromosome', 'start', 'end', 'n_bins']].reset_index(drop=True).equals(published_loci[['chromosome', 'start', 'end', 'n_bins']])}")

# %% -- Labels behind the published Supp Fig 3c and CNV scores
# Per-segment delta (difference of cluster means, as Tukey's meandiff) and CNV scores with
# the labels used and with ER_G_673a in cluster 2
RELABELLED = {"ER_G_673a": "2"}
labels_rerun = cluster_labels.copy()
for s_id, lab in RELABELLED.items():
    print(f"{s_id}: cluster {labels_rerun[s_id]} in the labels used here -> {lab}")
    labels_rerun[s_id] = lab

published_long = pd.read_csv(PUBLISHED_STRICT_LONG).set_index(["segment", "comp"])
rows = []
for seg, seg_bins in seg_coords.set_index("segment")["bins"].items():
    seg_bins = seg_bins.split(",")
    for lab_name, labs in (("labels used", cluster_labels), ("ER_G_673a in cluster 2", labels_rerun)):
        mu = cnv_df[seg_bins].groupby(labs.to_numpy()).mean()
        for k in COMPARISONS:
            d_bins = mu.loc[CL2] - mu.loc[k]
            rows.append(dict(segment=seg, comp=f"2 vs {k}", labels=lab_name,
                             delta=d_bins.mean(), abs_delta=d_bins.abs().max()))
seg_check = pd.DataFrame(rows).pivot_table(index=["segment", "comp"], columns="labels",
                                            values=["delta", "abs_delta"])
seg_check.columns = [f"{v} ({lab})" for v, lab in seg_check.columns]
seg_check["delta (published)"] = published_long["delta"]
seg_check["abs_delta (published)"] = published_long["abs_delta"]
print(seg_check.round(4).to_string())

score_check = pd.DataFrame({
    "labels used": cnv_df.abs().groupby(cluster_labels.to_numpy()).mean().mean(axis=1),
    "ER_G_673a in cluster 2": cnv_df.abs().groupby(labels_rerun.to_numpy()).mean().mean(axis=1),
    "published": pd.Series(PUBLISHED_CNV_SCORE)})
print(score_check.round(4))
output.save_table(OUT, "cluster2_labels_check", seg_check)

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

# -- read gtf: genes only (no exon/transcript), useful columns
genes_df = read_gtf_genes(GTF_PATH)

# -- intersection segments x genes
hits_df = genes_overlapping_segments(segments_df, genes_df)

hits_df = hits_df.sort_values(["Segment","Chromosome","Start"])

# -- list of genes per segment
genes_by_segment = (hits_df.groupby("Segment")["gene_name"]
                    .apply(lambda x: sorted(set(x)))
                    .reset_index(name="genes"))

genes_by_segment["n_genes"] = genes_by_segment["genes"].apply(len)

# %% -- Save the gene list inside the previous dataframe
genes_by_segment.rename(columns={"Segment": "segment"}, inplace=True)
genes_by_segment["segment"] = genes_by_segment["segment"].astype(int)
strict_summary["segment"] = strict_summary["segment"].astype(int)

merged_df = pd.merge(strict_summary, genes_by_segment, on="segment")

# %% -- Save the table for the supplementary results
output.save_table(OUT, "segment_cluster_genes", merged_df, index=False)
published_genes = pd.read_csv(PUBLISHED_SEGMENT_GENES)
print(merged_df[["segment", "chromosome", "start", "end", "n_bins", "n_genes"]])
print("n_genes per segment, published table:", published_genes["n_genes"].tolist())
print("Same gene lists as the published table:",
      merged_df["genes"].astype(str).tolist() == published_genes["genes"].tolist())

# %%
print(output.versions(["numpy", "pandas", "anndata", "scanpy", "infercnvpy", "statsmodels"]).to_string(index=False))
