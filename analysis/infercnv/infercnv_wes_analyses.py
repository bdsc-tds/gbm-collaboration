"""Script to run infercnv on WES CNV matrix."""

# To run successfully this script, you need to have pre-computed the
# CNV WES matrix of all samples across centers, and this should have been
# saved in:
# "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
# to do so, run PART 1 of the script "generate_cnv_wes_matrix.py"

# %% -- Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from typing import List
from pathlib import Path
import subprocess
import anndata as ad
import infercnvpy as cnv
import scanpy as sc

# Constants
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/infercnv/")

# %% -- Load the single-cell adata object saved in EFS
# Pre-step: load the adata and keep only the malignant cells
single_cell_obj = ad.read_h5ad("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/scRNAseq/infercnv/infercnv_adata.h5ad")

# append the wes cnv profile
path_wes = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
wes_profile = pd.read_csv(path_wes, index_col=0)

# %% -- Create a pseudobulk object
# Step 1: Get the sample id variable
sample_col = "orig.ident"

# Step 2: Aggregate X by group
# Make sure it's a dense matrix if it's sparse
X = single_cell_obj.X.toarray() if not isinstance(single_cell_obj.X, np.ndarray) else single_cell_obj.X
df = pd.DataFrame(X, index=single_cell_obj.obs[sample_col])

# Sum expression per gene per group
pseudobulk_matrix = df.groupby(df.index).sum().T  # genes x samples

# Step 3: Create a new AnnData object
pseudobulk_adata = ad.AnnData(X=pseudobulk_matrix.values.T)  # samples x genes
pseudobulk_adata.obs[sample_col] = pseudobulk_matrix.columns
pseudobulk_adata.var_names = single_cell_obj.var_names
pseudobulk_adata.obs_names = pseudobulk_matrix.columns
pseudobulk_adata.uns["cnv"] = single_cell_obj.uns["cnv"]
pseudobulk_adata.var = single_cell_obj.var

# %% - Keep only the common ids in wes_profile.index and pseudobulk_adata.obs["orig.ident"]
pb_ids = pseudobulk_adata.obs[sample_col].astype(str)
wes_ids = wes_profile.index.astype(str)

# 1. Check for missing samples
missing = pb_ids[~pb_ids.isin(wes_ids)]
if len(missing) > 0:
    raise ValueError(f"Missing WES profiles for samples: {missing.unique()}")
# Missing sample id ['ER_G_743a'] in WES data

# %% Keep only samples present in WES
keep_mask = pb_ids.isin(wes_ids)

print(f"Dropping {(~keep_mask).sum()} pseudobulk sample(s):",
      pb_ids[~keep_mask].unique())

# Subset pseudobulk
pseudobulk_adata = pseudobulk_adata[keep_mask].copy()

# 2. Check uniqueness
pb_ids = pseudobulk_adata.obs[sample_col].astype(str) # relist
assert pb_ids.is_unique, "Duplicate sample IDs in pseudobulk_adata.obs"
assert wes_profile.index.is_unique, "Duplicate sample IDs in wes_profile.index"

# %% -- Append the WES matrix profile to the object
pseudobulk_adata.obsm["X_cnv"] = wes_profile.loc[pseudobulk_adata.obs[sample_col]].values

# %% -- Start calculation and plotting
cnv.tl.pca(pseudobulk_adata)
cnv.pp.neighbors(pseudobulk_adata)
cnv.tl.leiden(pseudobulk_adata)

# %% -- Calculate dendogram
sc.tl.dendrogram(
    pseudobulk_adata,
    groupby="cnv_leiden",
    use_rep="X_cnv",
    cor_method="pearson",
    linkage_method="single",      # or single or average
    optimal_ordering=True
)

# %% -- Plotting
cnv.pl.chromosome_heatmap_summary(
    pseudobulk_adata,
    groupby="cnv_leiden",
    figsize=(15,5),
    dendrogram=True)

# %% -- Count how many samples found in each cluster
print(pseudobulk_adata.obs["cnv_leiden"].value_counts())
# cnv_leiden
# 0    65
# 1    48
# 2    38
# 3    36
# 4    32
# 5    28
# 6    18

# %% -- Save the Leiden cluster
leiden_cluster_baseline = pd.DataFrame(pseudobulk_adata.obs["cnv_leiden"])
leiden_cluster_baseline.to_csv(RESULTS_DIR / "leiden_cluster_wes_baseline_samples.csv", sep=",", index=True)

# %% -- Calculate CNV score on the WES profile
# Generate UMAP
cnv.tl.umap(pseudobulk_adata)
cnv.tl.cnv_score(pseudobulk_adata, groupby="cnv_leiden")

# Create composite figure
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(11, 5))
#ax4.axis("off")  # turn off the empty subplot

# Plot UMAPs into specific axes
cnv.pl.umap(
    pseudobulk_adata,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False
)
cnv.pl.umap(pseudobulk_adata, color="cnv_score", ax=ax2, show=True)

#####################
###### check previous cluster 1 with chr12 alterations #####
####################

# %%
pseudobulk_adata[pseudobulk_adata.obs["cnv_leiden"]=="1"].obs[sample_col].unique()

# %%
cnv.pl.chromosome_heatmap_summary(
    single_cell_obj[single_cell_obj.obs["orig.ident"].isin(['CH_G_159a', 'CH_G_168a', 'CH_G_169a', 'CH_G_172a', 'CH_G_173a'])],
    groupby="cnv_leiden",
    dendrogram=False,
    figsize=(15,10),
    show=False
)

# %%
single_cell_obj[single_cell_obj.obs["orig.ident"]=="CH_G_159a"].uns["cnv"]

# %%
list_gene_expression = single_cell_obj.var[
    (single_cell_obj.var["chromosome"]=="chr12") &
    (single_cell_obj.var["start"]>57530102) &
    (single_cell_obj.var["end"]<106510198)
]["gene"]

# %%

single_cell_153 = single_cell_obj[single_cell_obj.obs["orig.ident"]=="CH_G_153a"]
mean_expression = single_cell_153[:, list_gene_expression].X.mean(axis=0)
# %%
mean_expression_df = pd.DataFrame(mean_expression.T, index=list_gene_expression, columns=["mean_expr"])
# %%
mean_expression_df.sort_values(by="mean_expr", ascending=False)
# %%

# %% look at which tumor clones is associated with WES CNV profile 4
wes_cnv = "/home/sagemaker-user/gbm-collaboration/notebooks/leiden_cluster_wes_baseline_samples.csv"
wes_cnv_df = pd.read_csv(wes_cnv).set_index("orig.ident")
# list of samples in cluster 4
wes_cnv_df_cl4 = wes_cnv_df[wes_cnv_df["cnv_leiden"]==4].index.tolist()

# %%
single_cell_obj.obs["wes_cl4"] = "0"
single_cell_obj.obs.loc[
    single_cell_obj.obs["orig.ident"].isin(wes_cnv_df_cl4),
    "wes_cl4"
] = "1"

single_cell_obj_cl4 = single_cell_obj[single_cell_obj.obs["wes_cl4"]=="1"]

# %%
cnv.pl.chromosome_heatmap(single_cell_obj_cl4, groupby="cnv_leiden", dendrogram=False)

# %%
# %%
cluster = "25"
chromosome = "chr19"
#single_cell_obj = single_cell_obj.copy()
# Ensure genes have associated chromosome and position information
single_cell_obj.var['chromosome'] = single_cell_obj.var['chromosome'].astype(str)  # Ensure it's a string
# Subset genes on chromosome 2
chr_genes = single_cell_obj.var[single_cell_obj.var['chromosome'] == chromosome].copy()
# Get average expression for the cluster
cluster_cells = single_cell_obj.obs[single_cell_obj.obs['cnv_leiden'] == cluster].index
X_subset = single_cell_obj[cluster_cells, chr_genes.index].layers["LogNormalize"]
# Convert to dense array if sparse
if hasattr(X_subset, "toarray"):
    avg_expr = X_subset.toarray().mean(axis=0)  # Converts sparse to dense before computing mean
else:
    avg_expr = X_subset.mean(axis=0)  # Works for dense arrays

# Ensure it's a 1D array
avg_expr = np.ravel(avg_expr)  # Flattens any matrix into a vector
# Convert to DataFrame
chr_genes['avg_expr'] = pd.Series(avg_expr, index=chr_genes.index)
# Sort by genomic position
chr_genes = chr_genes.sort_values(by='start')

# Plot
plt.figure(figsize=(10, 5))
sns.lineplot(x=chr_genes['start'], y=chr_genes['avg_expr'], marker='o')
plt.ylim(0, 5)
plt.xlabel(f"Genomic Position on {chromosome} (bp)")
plt.ylabel("Average Expression")
plt.title(f"Chromosome {chromosome} Gene Expression in cluster {cluster}")
plt.grid(True)
plt.show()

# %%
# Step 1: Extract CNV matrix and create the DataFrame
cnv_matrix = single_cell_obj.obsm["X_cnv"]
cnv_df = pd.DataFrame.sparse.from_spmatrix(cnv_matrix, index=single_cell_obj.obs_names)

# Add cluster labels to the CNV DataFrame
cnv_df["cnv_leiden"] = single_cell_obj.obs["cnv_leiden"].values

# Step 2: Get bin-chromosome info
chr_pos = single_cell_obj.uns["cnv"]["chr_pos"]

# Initialize the bins_by_chr dictionary
bins_by_chr = {}

# Step 1: Initialize start_idx as 0 for chr1
start_idx = 0
last_bin = cnv_df.shape[1] - 1  # Get the last bin index

# Step 2: Loop through chr_pos and generate bins for each chromosome
chrom_list = list(chr_pos.keys())  # Extract the chromosome names (keys)
for i, chrom in enumerate(chrom_list):
    # Get the start index for the current chromosome
    start_idx = chr_pos[chrom]

    # For the last chromosome, set the end_idx as a reasonable value, e.g., len of all bins
    if i + 1 < len(chrom_list):
        # Get the start index of the next chromosome as the end index
        end_idx = chr_pos[chrom_list[i + 1]]
    else:
        # For the last chromosome, use the last bin index
        end_idx = last_bin

    # Generate bins for the current chromosome (start_idx to end_idx)
    bins_for_chrom = list(range(start_idx, end_idx))

    # Store the bins for the current chromosome
    bins_by_chr[chrom] = bins_for_chrom

# Step 5: Calculate mean CNV per cluster per chromosome
chrom_cluster_means = {}
for chrom, bins in bins_by_chr.items():
    # chrom_cnv = cnv_df.iloc[:, bins] # select bins + the cnv_leiden
    chrom_cnv = cnv_df.iloc[:, bins + [-1]]  # Select CNV values for the bins on this chromosome
    # Check if the selection works
    chrom_cluster_means[chrom] = chrom_cnv.groupby("cnv_leiden").mean().mean(axis=1)  # Mean per cluster

# Step 6: Convert to DataFrame for better visualization
chrom_cluster_df = pd.DataFrame(chrom_cluster_means)

# create an additional annotation for each cluster in index,
# that tells if there is an amplification in each column,
# if the value is above 0.1 or deletion if below -0.1
chrom_cluster_df_annot = chrom_cluster_df.copy()
# flatten the dataframe
chrom_cluster_df_annot = chrom_cluster_df_annot.stack().reset_index()
chrom_cluster_df_annot.columns = ["cnv_leiden", "chromosome", "mean"]
chrom_cluster_df_annot["cnv_type"] = ""

for row, value in chrom_cluster_df_annot.iterrows():
    if value["mean"] > 0.015:
        chrom_cluster_df_annot.loc[row, "cnv_type"] = "amplification"
    elif value["mean"] < -0.015:
        chrom_cluster_df_annot.loc[row, "cnv_type"] = "deletion"
    else:
        chrom_cluster_df_annot.loc[row, "cnv_type"] = "neutral"

# drop row that are neutral
chrom_cluster_df_annot = chrom_cluster_df_annot[chrom_cluster_df_annot["cnv_type"] != "neutral"]

count_cnv_type = chrom_cluster_df_annot.groupby(["cnv_leiden", "cnv_type"]).size().reset_index(name="count")

# %%
# reconstruct genes in bins
# Ensure gene info is in a DataFrame
gene_df = single_cell_obj.var
bin_df = pd.read_csv("/home/sagemaker-user/gbm-collaboration/notebooks/bin_positions.csv")

# Prepare output list
bin_info_list = []

# Loop through each bin and find matching genes
for _, bin_row in bin_df.iterrows():
    chrom = str(bin_row["chromosome"])
    start = bin_row["start"]
    end = bin_row["end"]
    bin_number = bin_row["bin"]

    # Select genes within this bin
    genes_in_bin = gene_df[
        (gene_df["chromosome"].astype(str) == chrom) &
        (gene_df["start"] >= start) &
        (gene_df["end"] <= end)
    ]["gene"].tolist()

    bin_info_list.append({
        "bin_number": bin_number,
        "chromosome": chrom,
        "start": start,
        "end": end,
        "genes": genes_in_bin
    })

# Final DataFrame
bin_gene_df = pd.DataFrame(bin_info_list)
# %%
# save
#bin_gene_df.to_csv("bins_genes.csv")
# %%
# Make sure bin_number is the correct column and type
genes_in_bin_1220 = bin_gene_df.loc[bin_gene_df["bin_number"] == 1220, "genes"].values[0]

# Print the full list (no truncation)
print(genes_in_bin_1220)

# %%
len(genes_in_bin_1220)

# %%
# Plot
plt.figure(figsize=(10, 5))
sns.lineplot(x=chr_genes[chr_genes["gene"].isin(genes_in_bin_1220)]['start'], y=chr_genes[chr_genes["gene"].isin(genes_in_bin_1220)]['avg_expr'], marker='o')
plt.ylim(0, 5)
plt.xlabel(f"Genomic Position on {chromosome} (bp)")
plt.ylabel("Average Expression")
plt.title(f"Chromosome {chromosome} Gene Expression in cluster {cluster}")
plt.grid(True)
plt.show()


# %%
# redo cnv clustering with heatmap grouping by wes
leiden_cluster_wes = pd.read_csv("/home/sagemaker-user/gbm-collaboration/gbmcollaboration/results/infercnv_tests/wes_input/leiden_cluster_wes_baseline_samples.csv", sep=",")
leiden_cluster_wes.rename(columns={"cnv_leiden": "cnv_wes"}, inplace=True)
leiden_cluster_wes["cnv_wes"] = leiden_cluster_wes["cnv_wes"].astype("category")

# %% merge with the single_cell_ob
single_cell_obj.obs = single_cell_obj.obs.merge(leiden_cluster_wes, right_on="orig.ident", left_on="orig.ident", how="outer")

# %%
# create subeet object for those without nans in cnv_leiden_y
single_cell_obj_cnv_wes = single_cell_obj[~single_cell_obj.obs["cnv_wes"].isna()]

# %%
sc.tl.dendrogram(single_cell_obj_cnv_wes, groupby="cnv_wes")

# %%
cnv.pl.chromosome_heatmap_summary(
    single_cell_obj_cnv_wes,
    groupby="cnv_wes",
    figsize=(10,5)
)

# %%
single_cell_obj_cnv_wes.obs["concat_leiden"] = single_cell_obj_cnv_wes.obs["cnv_wes"].astype(str) + "_" + single_cell_obj_cnv_wes.obs["cnv_leiden"].astype(str)
single_cell_obj_cnv_wes.obs["concat_leiden"] = single_cell_obj_cnv_wes.obs["concat_leiden"].astype("category")

# %%
sc.tl.dendrogram(single_cell_obj_cnv_wes, groupby="concat_leiden")
cnv.pl.chromosome_heatmap(
    single_cell_obj_cnv_wes,
    groupby="concat_leiden",
    figsize=(5,20)
)
# %%
pd.set_option('display.max_rows', None)
single_cell_obj_cnv_wes.obs["concat_leiden"].value_counts().sort_values(ascending=False)

# %%
pd.crosstab(single_cell_obj_cnv_wes.obs["cnv_leiden"], single_cell_obj_cnv_wes.obs["cnv_wes"])

# %%
top_n = 5
top_cnv_leiden = single_cell_obj_cnv_wes.obs.groupby("cnv_wes")["cnv_leiden"].value_counts().groupby(level=0).head(top_n).reset_index(name="count")
# Pivot the DataFrame to get cnv_wes as rows and cnv_leiden as columns
top_cnv_leiden_pivot = top_cnv_leiden.pivot(index="cnv_wes", columns="cnv_leiden", values="count").fillna(0)
# Display the pivoted DataFrame
top_cnv_leiden_pivot = top_cnv_leiden_pivot.astype(int)  # Convert to integer type for better readability
print(top_cnv_leiden_pivot)

# %% same but to avoid zeros in the other clusters where actually there are cells
# manually select the top clusters
top_cnv_leiden = single_cell_obj_cnv_wes.obs.groupby("cnv_wes")["cnv_leiden"].value_counts().reset_index(name="count")
top_cnv_leiden = top_cnv_leiden[top_cnv_leiden["cnv_leiden"].isin(["0","1","2","4","5","3","7","15","21","9"])]
# Pivot the DataFrame to get cnv_wes as rows and cnv_leiden as columns
top_cnv_leiden_pivot = top_cnv_leiden.pivot(index="cnv_wes", columns="cnv_leiden", values="count").fillna(0)
# Display the pivoted DataFrame
top_cnv_leiden_pivot = top_cnv_leiden_pivot.astype(int)  # Convert to integer type for better readability
print(top_cnv_leiden_pivot)


# %%
# do a visualization of the top cnv_leiden clusters in each cnv_wes cluster
import seaborn as sns
import matplotlib.pyplot as plt
# Create a heatmap of the top cnv_leiden clusters in each cnv_wes cluster
plt.figure(figsize=(12, 8))
sns.heatmap(
    top_cnv_leiden_pivot,
    annot=True,
    fmt="d",
    cmap="viridis",
    cbar_kws={'label': 'Count'},
    linewidths=0.5,
    linecolor='black'
)
plt.title("Top CNV Single-cell Clusters in Each CNV WES Cluster")
plt.xlabel("CNV Single-cell Clusters")
plt.ylabel("CNV WES Clusters")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# %%
# same figure but normalize for each cnv_wes cluster
top_cnv_leiden_normalized = top_cnv_leiden.copy()
# Normalize the counts by dividing by the sum of counts in each cnv_wes cluster
top_cnv_leiden_normalized["count"] = top_cnv_leiden_normalized.groupby("cnv_wes")["count"].transform(lambda x: x / x.sum())
# Pivot the DataFrame to get cnv_wes as rows and cnv_leiden as columns
top_cnv_leiden_normalized_pivot = top_cnv_leiden_normalized.pivot(index="cnv_wes", columns="cnv_leiden", values="count").fillna(0)
# Display the normalized pivoted DataFrame
top_cnv_leiden_normalized_pivot = top_cnv_leiden_normalized_pivot.astype(float)  # Convert to float type for better readability
# Create a heatmap of the normalized top cnv_leiden clusters in each cnv_wes cluster
plt.figure(figsize=(12, 8))
sns.heatmap(
    top_cnv_leiden_normalized_pivot,
    annot=True,
    fmt=".2f",
    cmap="viridis",
    cbar_kws={'label': 'Normalized Count'},
    linewidths=0.5,
    linecolor='black'
)
plt.title("Normalized Top CNV Leiden Clusters in Each CNV WES Cluster")
plt.xlabel("CNV Leiden Clusters")
plt.ylabel("CNV WES Clusters")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# %%
# now I am interested in particular for the leiden clusters 18, 28, 17
# 30, 14, 1, 19.
# Create a subset of the DataFrame for the specified CNV Leiden clusters
# and recompute the number of those clusters in each CNV WES cluster
selected_leiden_clusters = ["18", "28", "17", "30", "14", "1", "19"]
#single_cell_obj_cnv_wes.obs["cnv_leiden"] = single_cell_obj_cnv_wes.obs["cnv_leiden"].astype(int)

# restart from the beginning, not the top_cnv_leiden
subset_top_cnv_leiden = single_cell_obj_cnv_wes.obs[
    single_cell_obj_cnv_wes.obs["cnv_leiden"].isin(selected_leiden_clusters)
].groupby(["cnv_wes", "cnv_leiden"]).size().reset_index(name="count")
print(subset_top_cnv_leiden)

subset_top_cnv_leiden = subset_top_cnv_leiden[subset_top_cnv_leiden["cnv_leiden"].isin(selected_leiden_clusters)]

# %% Pivot the DataFrame to get CNV WES clusters as rows and CNV Leiden clusters as columns
subset_top_cnv_leiden_pivot = subset_top_cnv_leiden.pivot(index="cnv_wes", columns="cnv_leiden", values="count").fillna(0)
# Convert to integer type for better readability
subset_top_cnv_leiden_pivot = subset_top_cnv_leiden_pivot.astype(int)
# Display the pivoted DataFrame
print(subset_top_cnv_leiden_pivot)

# %%
# Create a heatmap of the selected CNV Leiden clusters in each CNV WES cluster
plt.figure(figsize=(12, 8))
sns.heatmap(
    subset_top_cnv_leiden_pivot,
    annot=True,
    fmt="d",
    cmap="viridis",
    cbar_kws={'label': 'Count'},
    linewidths=0.5,
    linecolor='black'
)
plt.title("Selected CNV Leiden Clusters in Each CNV WES Cluster")
plt.xlabel("CNV Leiden Clusters")
plt.ylabel("CNV WES Clusters")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
# %%
