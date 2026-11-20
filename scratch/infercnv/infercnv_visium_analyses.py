"""Script to infer CNV using inferCNVpy."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy
import infercnvpy as cnv
import gc
import time

from pathlib import Path
from gbmcollaboration import AllDataCenters, CHUVDataCenter, UKERDataCenter, CRDataCenter
from gbmcollaboration.utils.infercnv_utils import InferCNVAnalyzer


# %% Params
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/visium/infercnv")
METACONTROL_PATH = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/meta_control_v3.h5ad")
CELL_TYPE_KEY = "deconv_cohort_level_4"
REFERENCE_COLUMN = "malignant_status"
BATCH_KEY = 'sample'

THRESHOLD_N_GENES = 100
THRESHOLD_N_COUNTS = 1_000_000

# %% Create results dir
RESULTS_DIR.mkdir(parents=True,exist_ok=True)

# %% -- Load data
adata = AllDataCenters.load_spatial_adata(
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
    list_sample_id=None,
    normalization="raw_count",
    radius=1,
    deconvolution_slot_name=CELL_TYPE_KEY,
)

# ensure raw counts are in .X.
# clear log1p status if present
if 'log1p' in adata.uns:
    del adata.uns['log1p']
adata.X = adata.layers['raw_count']

# %% Load meta-control spots info
adata_ctrl = sc.read_h5ad(METACONTROL_PATH)
control_spots = adata_ctrl.obs_names
control_spots_filtered = adata_ctrl.obs_names[~adata_ctrl.obs_names.isin(adata.obs_names)]
adata = sc.concat([adata,adata_ctrl[control_spots_filtered]])

adata.obs[REFERENCE_COLUMN] = 'malignant'
adata.obs.loc[control_spots,REFERENCE_COLUMN] = 'control'
print(f"Using {len(control_spots)} control spots")
# %% -- Initialize the InferCNVAnalyzer class
infercnv_obj = InferCNVAnalyzer(adata=adata)

# add a "reference" column to add "low_genes" and "high_counts"
# for downstream filtering
infercnv_obj.add_reference_column(
    cell_type_annotation_key=REFERENCE_COLUMN)

# %%
# Define the reference categories as the ones that are not malignant
list_reference_cat = ['control']

# %% -- Run infercnv inference
start_time = time.time()
infercnv_obj.run_infercnv(list_reference_cat=list_reference_cat)

infercnv_obj.run_leiden_clustering(resolution=0.4)

end_time = time.time()
print(f"Total time running: {round((end_time-start_time)/60)} minutes")

# %% -- Assign each cnv cluster a CNV score
infercnv_obj.run_cnv_score()

# %% -- Run ITH scoring
start_time = time.time()
infercnv_obj.run_ith_scores(groupby=BATCH_KEY, save_path=RESULTS_DIR / "ith_scores.csv")
end_time = time.time()
print(f"Total time running: {round((end_time-start_time)/60)} minutes")
# %% -- Perform visualizations
# infercnv_obj.composite_plots(save_path=RESULTS_DIR / "cnv_plots.png")
# infercnv_obj.adata.obsm['X_umap']=infercnv_obj.adata.obsm['X_cnv_umap']

# %% -- Create the annotation file
infercnv_obj.tumor_clone_annotation(save_path=RESULTS_DIR / "tumor_clone_annotation.csv")

# %% -- Save the summary heatmap
cnv.pl.chromosome_heatmap_summary(
    infercnv_obj.adata,
    groupby="cnv_leiden",
    dendrogram=True,
    show=False
)
fig = plt.gcf()
fig.savefig(RESULTS_DIR / "cnv_heatmap.png", dpi=300, bbox_inches='tight')
plt.close(fig)

# %% -- Create gene expression UMAPs colored by cell type or by cluster
# infercnv_obj.plot_gene_exp_umaps(
#     cell_type_lvl=CELL_TYPE_KEY,
#     save_path=RESULTS_DIR / f"umaps_{CELL_TYPE_KEY}.png",
#     reset_umap=False # reset to False the second time
# )

# %% -- Recompute genomic window and derive start-end position
infercnv_obj.get_bin_positions(
    save_path=RESULTS_DIR / "genomic_windows.csv"
)

# %%
infercnv_obj.adata.layers.clear()
infercnv_obj.adata.X = scipy.sparse.csr_matrix(infercnv_obj.adata.X.shape)
infercnv_obj.adata.write(RESULTS_DIR / "infercnv_adata.h5ad")
