"""Script to generate genomic bins and compute average log2 copy ratio per bin from WES CNV data."""

# --------------------PART 1--------------------------------
# %% -- Imports packages, load constants, functions and class
import pandas as pd
import numpy as np
import os
from typing import List
from pathlib import Path
import subprocess
from gbmcollaboration.definitions import (
    PATIENT_EXCLUSION,
    RELAPSE_EXCLUSION,
    PATIENTS_IDH1_MUT,
)

# Constants
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/")
PATH_BINS_SINGLECELL = RESULTS_DIR / "scRNAseq/infercnv/genomic_windows.csv"
BUCKET_NAME="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias"
PREFIX_UKER="3205eea0-5076-4294-8876-08cf98474603"
PREFIX_CHUV="a3857b67-3f28-40e6-b1f9-4842a2ad49af"
PREFIX_CR="b0d0f0c7-a904-4331-ad1c-f86ca0527120"

# create a class that will perform the binning and averaging
class CNVBinAverager:
    """Class to bin CNV data and compute average log2 copy ratio per bin.

    Arguments:
    ----------
    bin_path : str
        Path to the CSV file containing bin definitions with columns:
        'chromosome', 'start', 'end', 'bin'.

    Returns:
    --------
    pd.DataFrame
        DataFrame with samples as rows and bins as columns containing average log2 copy ratios.
    """
    def __init__(self, bin_path: str):
        self.bins = pd.read_csv(bin_path)
        self.bins["chromosome"] = self.bins["chromosome"].astype(str)

    def _process_sample(self, sample_path: str) -> pd.Series:
        sample_id = os.path.basename(sample_path).replace(".cnr", "")
        cnr = pd.read_csv(sample_path, sep="\t")
        cnr = cnr[cnr["gene"] != "Antitarget"].copy()

        cnr["mid"] = (cnr["start"] + cnr["end"]) // 2
        cnr["chromosome"] = cnr["chromosome"].astype(str)
        cnr["bin_id"] = np.nan

        # Assign each cnr row to a bin
        for _, row in self.bins.iterrows():
            mask = (
                (cnr["chromosome"] == row["chromosome"]) &
                (cnr["start"] >= row["start"]) &
                (cnr["end"] <= row["end"])
            )
            cnr.loc[mask, "bin_id"] = row["bin"]

        # Drop and average
        cnr_binned = cnr.dropna(subset=["bin_id"]).copy()
        cnr_binned["bin_id"] = cnr_binned["bin_id"].astype(int)
        avg_log2 = cnr_binned.groupby("bin_id")["log2"].mean()

        # Return a dataframe with bin_id and log2 and sample_id
        return avg_log2.rename(sample_id)

    def process_samples(self, sample_paths: List[str]) -> pd.DataFrame:
        matrix = []
        for path in sample_paths:
            series = self._process_sample(path)
            matrix.append(series)

        # Combine and reindex columns to cover all bins
        result_df = pd.DataFrame(matrix)
        result_df = result_df.sort_index(axis=1)
        result_df.columns.name = "bin_id"
        return result_df

def list_wes_cnr_files(bucket_name, search_path=".cnr", prefix=""):
    try:
        command = f"aws s3 ls s3://{bucket_name}/{prefix}/WES/CNVs_CNVkit/ --recursive"
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        matching_files = []
        for line in result.stdout.splitlines():
            parts = line.strip().split()
            if len(parts) == 4:  # aws s3 ls --recursive gives lines like: date time size path
                filepath = parts[-1]
                if filepath.endswith(search_path):
                    full_path = f"s3://{bucket_name}/{filepath}"
                    matching_files.append(full_path)
        return matching_files
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        print(e.output)
        return []

# %% -- Generate sample bins path
sample_paths_uker = list_wes_cnr_files(BUCKET_NAME, prefix=PREFIX_UKER)
sample_paths_chuv = list_wes_cnr_files(BUCKET_NAME, prefix=PREFIX_CHUV)
sample_paths_cr = list_wes_cnr_files(BUCKET_NAME, prefix=PREFIX_CR)
# combine all sample paths
sample_paths = sample_paths_uker + sample_paths_chuv + sample_paths_cr

# %% -- Run CNVBinAverager
# record time
import time
time_start = time.time()
averager = CNVBinAverager(PATH_BINS_SINGLECELL)
log2_matrix = averager.process_samples(sample_paths)
end_time = time.time()
print(f"Time taken: {end_time - time_start:.2f} seconds")
# It takes approximately 73 minutes to run

# %% -- Post processing of the constructed log2 matrix
# drop "_S000021" on each sample id in index
log2_matrix.index = log2_matrix.index.str.replace("_S000021", "", regex=False)

# drop patients or samples using pre-defined list
patients_set = (
    set(PATIENT_EXCLUSION)
    | set(RELAPSE_EXCLUSION)
    | set(PATIENTS_IDH1_MUT)
)
pattern = "|".join(patients_set)

idx_matches = log2_matrix.index[
    log2_matrix.index.astype(str).str.contains(pattern)
]

if not idx_matches.empty:
    print(f"Excluding {len(idx_matches)} rows.")
    log2_matrix = log2_matrix.drop(index=idx_matches)

# load bins for exploring
bins = pd.read_csv(PATH_BINS_SINGLECELL)
# Exclude bins that are unique in some chromosome
bins.chromosome.value_counts()
# here chromosome 21 and Y only have one bin, and add noise
# in the heatmap.

# %%
bins_to_exclude = bins[bins["chromosome"].isin(["chr21","chrY"])][["bin"]].index.tolist()
log2_matrix_curated = log2_matrix.drop(columns=bins_to_exclude)

# %% -- Save or explore
log2_matrix_curated.to_csv(
    "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
    )

# --------------------PART 2--------------------------------
# %% -- Exploration of the data to find any issue
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

# load data
log2_matrix = pd.read_csv("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv")
log2_matrix = log2_matrix.set_index("Unnamed: 0")

bins = pd.read_csv(PATH_BINS_SINGLECELL)
bins.drop(columns=["Unnamed: 0"], inplace=True)

# %% -- Calculate and Plot the distribution of CNV scores
cnv_score = log2_matrix.abs().sum(axis=1)

plt.figure(figsize=(10, 5))
plt.hist(cnv_score, bins=50, color="blue", alpha=0.7)
plt.title("Distribution of CNV scores")
plt.xlabel("CNV score")
plt.ylabel("Frequency")
plt.axvline(cnv_score.mean(), color="red", linestyle="--", label="Mean")
plt.axvline(cnv_score.median(), color="green", linestyle="--", label="Median")
plt.legend()
plt.show()

# %% -- Plot a PCA of the CNV scores
# scale the data first
scaler = StandardScaler()
log2_matrix_scaled = scaler.fit_transform(log2_matrix)
# Perform PCA
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(log2_matrix_scaled)

plt.scatter(pca_coords[:, 0], pca_coords[:, 1],
            c=cnv_score, cmap='viridis', s=50, alpha=0.7)
plt.title("PCA of WES-derived CNV profiles")
plt.xlabel("PC1")
plt.ylabel("PC2")

# %% find who has PC1 loading below -20
df_pca = pd.DataFrame(pca_coords, columns=["PC1", "PC2"], index=log2_matrix.index)
df_pca["cnv_score"] = cnv_score

sample_below_threshold = df_pca[df_pca["PC1"] < -20].index.tolist()
# now look for those with high PCA2
df_pca.loc[sample_below_threshold]

# %% -- Explore the CNV profiles of some outliers on PCA
# plot the heatmap for sample CR_G_315a
sample_id = "CH_G_192a" #"CR_G_366a" "ER_G_196a" "CR_G_315a" "CH_G_192a"
sample_log2 = log2_matrix.loc[sample_id].values
sample_bins = log2_matrix.columns
sample_bins_df = bins[bins["bin"].isin(sample_bins)].copy()
sample_bins_df["log2"] = sample_log2
sample_bins_df = sample_bins_df.sort_values(["bin"])
# Create a matrix for heatmap: single row with log2 per bin
heatmap_data = pd.DataFrame([sample_bins_df["log2"].values], index=[sample_id])
plt.figure(figsize=(16, 1.5))
sns.heatmap(
    heatmap_data,
    cmap="coolwarm",  # diverging colormap
    center=0,         # 0 log2 = neutral
    cbar_kws={"label": "log2 copy ratio"},
    xticklabels=False,
    yticklabels=True
)
# Add vertical lines for chromosome boundaries
chrom_bounds = sample_bins_df.groupby("chromosome")["bin"].agg(["min", "max"])
for idx in chrom_bounds["min"].values:
    plt.axvline(x=idx, color="black", linestyle="--", linewidth=0.5)
# Set chromosome labels
plt.xticks(
    ticks=chrom_bounds["min"] + (chrom_bounds["max"] - chrom_bounds["min"]) // 2,
    labels=chrom_bounds.index,
    rotation=90
)
plt.title(f"CNV heatmap per bin ({sample_id})")
plt.tight_layout()
plt.show()

# %% -- Same plot but with scatter
plt.figure(figsize=(16, 4))
sns.scatterplot(data=sample_bins_df, x="bin", y="log2", color="black", s=10, linewidth=0)
# Add vertical lines for chromosome boundaries
for idx in chrom_bounds["min"].values:
    plt.axvline(x=idx, color="grey", linestyle="--", alpha=0.4)
# Chromosome ticks as labels
plt.xticks(chrom_bounds["min"].values, chrom_bounds.index, rotation=90)
# Formatting
plt.title(f"Binned CNV signal across genome from WES ({sample_id})")
plt.xlabel("Chromosome")
plt.ylabel("Average log2 copy ratio")
plt.legend().remove()  # 🔥 Remove legend
plt.tight_layout()
# draw a horizontal line at 0
plt.axhline(y=0, color='blue', linestyle='--', linewidth=1)
plt.axhline(y=-1, color='red', linestyle='--', linewidth=1)
#plt.axhline(y=1, color='red', linestyle='--', linewidth=1)
plt.show()

# %% -- Plot CNV heatmap grouping samples by profiles
def plot_cnv_clustermap(
    log2_matrix,
    figsize=(12, 12),
    cmap="coolwarm",
    center=0,
    q=(0.02, 0.98),          # quantiles pour l'échelle
    metric="correlation",    # souvent mieux que euclidean pour profils
    method="average",
    col_cluster=False,
    row_cluster=True,
):
    # bornes robustes
    lo, hi = np.nanquantile(log2_matrix.to_numpy(), q)
    # symétriser autour de 0 pour une palette divergente
    M = max(abs(lo-center), abs(hi-center))
    vmin, vmax = center - M, center + M

    g = sns.clustermap(
        log2_matrix,
        metric=metric,
        method=method,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        col_cluster=col_cluster,
        row_cluster=row_cluster,
        xticklabels=False,
        yticklabels=True,
        figsize=figsize,
        cbar_kws={"label": "log2 ratio"},
    )
    g.ax_heatmap.set_title("Hierarchical Clustering of CNV Profiles", pad=20)
    plt.show()

# exemple
plot_cnv_clustermap(log2_matrix, figsize=(10, 14), q=(0.02, 0.98))

# %% -- Quality control
qc = pd.DataFrame(index=log2_matrix.index)

# Variabilité globale
qc["sd"] = log2_matrix.std(axis=1)
# Amplitude extrême
qc["mad"] = (log2_matrix - log2_matrix.median(axis=1).values[:, None]).abs().median(axis=1)
# Proportion de segments extrêmes
qc["prop_extreme"] = (log2_matrix.abs() > 1).mean(axis=1)

qc.describe()

# %% -- Find outlier samples
qc_z = qc.apply(zscore, nan_policy="omit")
qc["outlier_score"] = qc_z.abs().sum(axis=1)

# seuil simple
threshold = qc["outlier_score"].quantile(0.95)
outliers = qc[qc["outlier_score"] > threshold]

outliers

# %% -- Plot the outlier
qc = qc.copy()
qc["is_outlier"] = qc["outlier_score"] > threshold   # bool

sns.pairplot(
    qc,
    vars=["sd", "mad", "prop_extreme"],
    hue="is_outlier",
    plot_kws={"alpha": 0.6, "s": 20},
    diag_kind="hist"
)
plt.show()

# %%
