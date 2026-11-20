"""Pipeline: WES CNV matrix (frozen as results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv).

Output: for each WES sample and each genomic bin of the single-cell inferCNV
(genomic_windows.csv, 200-gene windows), the mean CNVkit log2 copy ratio of the target bins
(.cnr, antitargets removed) lying entirely inside the bin. Frozen file: 272 samples x 1,344
bins, read by fig3b_suppfig3ac_infercnv_wes_analyses.py and the data release.

Origin: analysis/infercnv/generate_cnv_wes_matrix.py, Part 1 (Part 2, exploration of the
matrix, is not part of the pipeline).

Inputs: the CNVkit .cnr files of the three centres on S3 (<prefix>/WES/CNVs_CNVkit/, needs
AWS credentials), the bins written by fig3b_0_infercnv_sc.py
(results/scRNAseq/infercnv/genomic_windows.csv).

Runtime: about 73 minutes (loop over the 1,346 bins for each of ~300 files), < 10 GB RAM.
Not re-run for the paper.

Comparison with the frozen matrix:
- The code drops chrX and chrY (`bins_to_exclude`). The frozen matrix lacks chr21 (bin 1274)
  and chrY (bin 1345) and keeps chrX (bins 1296-1344), so it was written by an earlier
  version of this line. chrX is dropped again before clustering and plotting (Fig 3b), so
  chr21 is absent from the clustering and from Fig 3b.
- Samples are excluded by substring match of PATIENT_EXCLUSION, RELAPSE_EXCLUSION and
  PATIENTS_IDH1_MUT; `_S000021` is stripped from the IDs.
"""

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
from gbmcollaboration.paper import output

# Constants
OUT = output.setup("fig3b_0_cnv_wes_matrix")
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
# in the heatmap. But we keep 21 because it's difficult to justify the removal
# we also want to exclude chrX

# %%
bins_to_exclude = bins[bins["chromosome"].isin(["chrY", "chrX"])][["bin"]].index.tolist()
log2_matrix_curated = log2_matrix.drop(columns=bins_to_exclude)

# %% -- Save or explore
log2_matrix_curated.to_csv(OUT / "cnv_wes_matrix_full_cohorts.csv")
