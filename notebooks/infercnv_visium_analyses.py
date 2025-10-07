"""Script to infer CNV using inferCNVpy."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import infercnvpy as cnv
import gc

from pathlib import Path
from gbmcollaboration import CHUVDataCenter, UKERDataCenter
from gbmcollaboration.utils.infercnv_utils_fixed import InferCNVAnalyzer

def load_visium_adata(
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
    list_sample_id=None,
    normalization="log_norm_count",
    radius=1,
    deconvolution_slot_name="deconv_cohort_level_4_2steps",
    ):
    """Load entire Visium cohort as anndata."""

    # Define data centers and common loading parameters
    data_centers = [CHUVDataCenter, UKERDataCenter]

    # Load all data
    loaded_objects = []
    for data_center in data_centers:
        print('Loading',data_center)

        loaded_objects.append(data_center.load_spatial(
            exclude_low_quality_visium=exclude_low_quality_visium,
            exclude_idh1_mut=exclude_idh1_mut,
            list_sample_id=list_sample_id,
            normalization=normalization,
            radius=radius,
            deconvolution_slot_name=deconvolution_slot_name,
        ))

    # Concatenate the dataframes
    counts_df = pd.concat([obj.genecounts_df for obj in loaded_objects])
    obs_df = pd.concat([obj.celltype_df for obj in loaded_objects])

    del loaded_objects
    gc.collect()

    # Align obs to counts and create the 'sample' column
    obs_df = obs_df.loc[counts_df.index]
    obs_df['sample'] = obs_df.index.get_level_values(0)

    # Flatten the MultiIndex for AnnData compatibility
    flat_index = counts_df.index.map(lambda idx: f"{idx[0]}_{idx[1]}")
    obs_df.index = flat_index
    counts_df.index = flat_index

    counts_df.sum(axis=0).sort_values()

    # Create the final AnnData object
    adata = sc.AnnData(
        obs=obs_df,
        var=pd.DataFrame(index=counts_df.columns),
        layers={normalization:counts_df.values}
    )

    del counts_df, obs_df
    gc.collect()

    # --- Output ---
    print(" Final adata created:")
    print(" • shape:", adata.shape)
    print(" • layers:", list(adata.layers.keys()))
    print(" • obs columns:", adata.obs.columns.tolist())
    print(" • samples:", adata.obs['sample'].unique().tolist())

    return adata

# %% Params
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/visium/infercnv")
METACONTROL_PATH = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/meta_control_v2.h5ad")
CELL_TYPE_KEY = "deconv_cohort_level_4_2steps"
REFERENCE_COLUMN = "malignant_status"

# /!\ to adjust later
THRESHOLD_N_GENES = 100
THRESHOLD_N_COUNTS = 1_000_000

# %% Create results dir
RESULTS_DIR.mkdir(parents=True,exist_ok=True)

# %% -- Load data
adata = load_visium_adata(
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

genes_to_keep = adata.var_names[~np.isnan(adata.X).any(0)]
genes_to_drop = adata.var_names[np.isnan(adata.X).any(0)]
adata._inplace_subset_var(genes_to_keep)
print(f"Dropped nan genes:",*genes_to_drop)

# %% Load meta-control spots info
adata_ctrl = sc.read_h5ad(METACONTROL_PATH)
healthy_spots = adata_ctrl.obs_names

adata.obs[REFERENCE_COLUMN] = 'malignant'
adata.obs.loc[healthy_spots,REFERENCE_COLUMN] = 'healthy'
print(f"Using {len(healthy_spots)} healthy spots")


# %% -- Initialize the InferCNVAnalyzer class
infercnv_obj = InferCNVAnalyzer(adata=adata)

# add a "reference" column to add "low_genes" and "high_counts"
# for downstream filtering
infercnv_obj.add_reference_column(
    cell_type_annotation_key=REFERENCE_COLUMN)

# %%
# Define the reference categories as the ones that are not malignant
list_reference_cat = ['healthy']

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
# infercnv_obj.adata.obsm['X_umap']=infercnv_obj.adata.obsm['X_cnv_umap']

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
infercnv_obj.plot_gene_exp_umaps(
    cell_type_lvl=CELL_TYPE_KEY,
    save_path=RESULTS_DIR / f"umaps_{CELL_TYPE_KEY}.png",
    reset_umap=False # reset to False the second time
)

# %% -- Recompute genomic window and derive start-end position
infercnv_obj.get_bin_positions(
    save_path=RESULTS_DIR / "genomic_windows.csv"
)
