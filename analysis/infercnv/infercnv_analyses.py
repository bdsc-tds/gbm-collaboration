"""Script to infer CNV using inferCNVpy."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import infercnvpy as cnv
from pathlib import Path

from gbmcollaboration import AllDataCenters, CHUVDataCenter
from gbmcollaboration.utils.infercnv_utils import InferCNVAnalyzer
from gbmcollaboration import (
    PATIENTS_IDH1_MUT,
)


# %% Params
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/scRNAseq/infercnv")
REFERENCE_COLUMN = "celltype_level1_scanvi"

# %% Create results dir
RESULTS_DIR.mkdir(parents=True,exist_ok=True)

# %%-- Load data
data = AllDataCenters.load_singlecell(exclude_idh1_mut=True)
# data['mosaic_cr'].obs['celltype_level4_curated'] = data['mosaic_cr'].obs['celltype_level4_scanvi']
adata = sc.concat(data,label='center')
del data

# ensure raw counts are in .X.
# clear log1p status if present
if 'log1p' in adata.uns:
    del adata.uns['log1p']
adata.X = adata.layers['unnormalized_counts']

# %% -- Initialize the InferCNVAnalyzer class
infercnv_obj = InferCNVAnalyzer(adata=adata)

# add a "reference" column to add "low_genes" and "high_counts"
# for downstream filtering
infercnv_obj.add_reference_column(cell_type_annotation_key=REFERENCE_COLUMN)

# %%
# Define the reference categories as the ones that are not malignant
list_reference_cat = adata[
    ~adata.obs[REFERENCE_COLUMN].isin(["Malignant", "Doublet", "Other", "Unknown"])
].obs[REFERENCE_COLUMN].unique()

# %% -- Run infercnv inference
import time
start_time = time.time()
infercnv_obj.run_infercnv(list_reference_cat=list_reference_cat)

infercnv_obj.run_leiden_clustering(resolution=0.4)

end_time = time.time()
print(f"Total time running: {round((end_time-start_time)/60)} minutes")

# %% -- Assign each cnv cluster a CNV score
infercnv_obj.compute_cnv_score()

# %% -- Perform visualizations
infercnv_obj.composite_plots(save_path=RESULTS_DIR / "cnv_plots.png")

# %% -- Create the annotation file
infercnv_obj.tumor_clone_annotation(save_path=RESULTS_DIR / "tumor_clone_annotation.csv")

# %% -- Save the summary heatmap
cnv.pl.chromosome_heatmap_summary(
    infercnv_obj.adata,
    groupby="cnv_leiden",
    dendrogram=True
)
fig = plt.gcf()
fig.savefig(RESULTS_DIR / "cnv_heatmap.png", dpi=300, bbox_inches='tight')
plt.close(fig)

# %% -- Create gene expression UMAPs colored by cell type or by cluster
# infercnv_obj.plot_gene_exp_umaps(
#     cell_type_lvl="celltype_level2_scanvi",
#     save_path=RESULTS_DIR / "umaps_lvl2_CNVcluster.png",
#     reset_umap=False # reset to False the second time
# )

# %% -- Recompute genomic window and derive start-end position
infercnv_obj.get_bin_positions(
    save_path=RESULTS_DIR / "genomic_windows.csv"
)

# %%
infercnv_obj.adata.write(RESULTS_DIR / "infercnv_adata.h5ad")
