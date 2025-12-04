"""Script to explore the bulkRNAseq data."""

# %% -- Load package data
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from sklearn.decomposition import PCA
from umap import UMAP


from gbmcollaboration import CHUVDataCenter, UKERDataCenter, AllDataCenters

# %% -- Define constants
# Choose your normalization method
NORMALIZATION_METHOD = "TPM"  # "TPM" or "normalized counts" or "raw counts"
clinical_variables_of_interest = ["Sample collection chronology", "MGMT promoter methylation", "os_years", "os_censor", "Sample source", "GB - Primary tumour side", "Number of brain tumour sites", "GB-Largest diameter of the primary tumour (mm)", "years_sample_harvest_dearchival"]

# Set the display options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

# %% -- Load the data across centers
all_centers = AllDataCenters

# or just one center at a time
chuv_center = CHUVDataCenter
uker_center = UKERDataCenter

chuv_dict = chuv_center.load_tabular(exclude_idh1_mut=True)
chuv_counts=chuv_dict["bulk_rna"][NORMALIZATION_METHOD]

uker_dict = uker_center.load_tabular(exclude_idh1_mut=True)
uker_counts=uker_dict["bulk_rna"][NORMALIZATION_METHOD]

print(f"CHUV df size: {chuv_counts.shape} and UKER df size: {uker_counts.shape}")

# Merge the counts data
counts_df = chuv_counts.merge(uker_counts, left_index=True, right_index=True)
counts_df = counts_df.T
counts_df.index = counts_df.index.str.replace("_mRNA", "")

# Add a column to indicate the center
counts_df["center"] = counts_df.index.str.split("_").str[0]

# %% -- Load the clinical data and merge with the gene expression data
chuv_clin=chuv_dict["clinical"]["processed gbm clinical"].reset_index().set_index("Tumour block ID for MOSAIC")
chuv_clin=chuv_clin[clinical_variables_of_interest]
uker_clin=uker_dict["clinical"]["processed gbm clinical"].reset_index().set_index("Tumour block ID for MOSAIC")
uker_clin=uker_clin[clinical_variables_of_interest]

clin_df = pd.concat([chuv_clin, uker_clin])
print(f"Clinical data size: {clin_df.shape}")

# Merge the clin info
counts_clinical = counts_df.merge(clin_df, left_index=True, right_index=True)
print(f"Merged data size: {counts_clinical.shape}")

# %% -- Explore the data with some PCA and UMAP plots
pca = PCA(n_components=6)

# drop the center and clinical variables of interest
drop_cols = ["center"] + clinical_variables_of_interest
data_pca = pca.fit_transform(counts_clinical.drop(columns=drop_cols))

# Calculate % of variance explained by each component
print(f"Explained variance: {pca.explained_variance_ratio_*100}")

# Plot the first two components
plt.scatter(
    data_pca[:, 0],
    data_pca[:, 1],
    s=10,
    c=counts_clinical["center"].map({"CH": "red", "ER": "blue"}))
plt.legend(
    handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='CH'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='ER'),
    ],
    title="Center")
plt.xlabel('First PCA component')
plt.ylabel('Second PCA component')
plt.show()

# %% -- Get the genes that have the highest loading on the first PCA axis
loading = pca.components_[0]
loading_df = pd.DataFrame(loading, index=counts_clinical.drop(columns=drop_cols).columns, columns=["loading"])
loading_df.sort_values("loading", ascending=False).head(10)

# %% -- UMAP
umap_model = UMAP(
    n_neighbors=50,        # Try larger values for tighter clusters, e.g., 30-100
    min_dist=0.01,         # Smaller min_dist for tighter clusters
    spread=1.2,            # Spread can help separate clusters
    n_epochs=500,          # More epochs for thorough optimization
    init='spectral',       # Spectral init for better stability
    random_state=42        # Ensures reproducibility
)

# Fit UMAP to the data
umap_embedding = umap_model.fit_transform(counts_clinical.drop(columns=drop_cols))

# Plot the UMAP embedding
plt.scatter(
    umap_embedding[:, 0],
    umap_embedding[:, 1],
    s=10,
    c=counts_clinical["center"].map({"CH": "red", "ER": "blue"}))
plt.legend(
    handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='CH'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='ER'),
    ],
    title="Center"
)

plt.show()

# %% -- Take top 100 variable genes and calculate correlation to see if there are clusters
top100 = counts_clinical.drop(columns=drop_cols).std().sort_values(ascending=False).head(100)
counts_df_100 = counts_clinical.loc[:, top100.index.tolist() + drop_cols]

# Calculate correlation of gene expression between samples
correlation = counts_df_100.drop(columns=drop_cols).T.corr()

# cluster map of the correlation
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    # color by center
    row_colors=counts_df_100["center"].map({"CH": "red", "ER": "blue"})
)

# %% -- Clustermap with MGMT promoter methylation
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=counts_df_100["MGMT promoter methylation"].map({"Methylated": "red", "Unmethylated": "blue"})
)
# %% -- Clustermap with baseline vs recurrent
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=counts_df_100["Sample collection chronology"].map({"At recurrence": "red", "Baseline": "blue"})
)
# %% -- Clustermap with OS
# Normalize OS years for color mapping
norm_OS = plt.Normalize(counts_df_100["os_years"].min(), counts_df_100["os_years"].max())
os_colors = sns.color_palette("viridis", as_cmap=True)(norm_OS(counts_df_100["os_years"]))

# Convert os_colors to a 1D list of hex colors
os_colors = [to_hex(color) for color in os_colors]

# Map os_censor to discrete colors
censorship_palette = {1: "blue", 0: "red"}  # 1 for non-censored, 0 for censored
censorship_colors = counts_df_100["os_censor"].map(censorship_palette)

# Combine both row colors into a DataFrame
row_colors = pd.DataFrame({
    "OS (years)": os_colors,
    "Censorship": censorship_colors
})

# Plot the clustermap
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=row_colors
)
plt.show()

# %% -- Clustermap with biopsy vs surgical specimen
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=counts_df_100["Sample source"].map({"Biopsy": "red", "Surgical specimen": "blue"})
)

# %% -- Clustermap with primary tumour side
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=counts_df_100["GB - Primary tumour side"].map({"Right": "red", "Left": "blue", "Bilateral": "green"})
)

# %% -- Clustermap with number of brain tumour sites
sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    # color by cancer type
    row_colors=counts_df_100["Number of brain tumour sites"].map({3: "red", 1: "blue", 2: "green"})
)
# %% -- Clustermap with tumour size
norm_size = plt.Normalize(counts_df_100["GB-Largest diameter of the primary tumour (mm)"].min(), counts_df_100["GB-Largest diameter of the primary tumour (mm)"].max())
row_colors = sns.color_palette("viridis", as_cmap=True)(norm_size(counts_df_100["GB-Largest diameter of the primary tumour (mm)"]))

sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=row_colors)

# %% -- Clustermap with years between sample harvest and dearchival
norm_dearchival = plt.Normalize(counts_df_100["years_sample_harvest_dearchival"].min(), counts_df_100["years_sample_harvest_dearchival"].max())
row_colors = sns.color_palette("viridis", as_cmap=True)(norm_dearchival(counts_df_100["years_sample_harvest_dearchival"]))

sns.clustermap(
    correlation,
    figsize=(10, 10),
    method="ward",
    metric="euclidean",
    row_cluster=True,
    col_cluster=True,
    row_colors=row_colors)
