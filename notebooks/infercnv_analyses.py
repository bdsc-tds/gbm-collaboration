"""Script to infer CNV using inferCNVpy."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import infercnvpy as cnv

from gbmcollaboration import CHUVDataCenter, UKERDataCenter
from gbmcollaboration.utils.infercnv_utils import InferCNVAnalyzer
# from gbmcollaboration.viz import plot_sankey_plot

# %% -- Load data
centers = ["CHUV", "UKER"]
for center in centers:
    if center == "UKER":
        data_center = UKERDataCenter
        single_cell_obj_uker = data_center.load_singlecell()
        single_cell_obj_uker.obs["center"] = "UKER"
        # drop samples
        sample_to_exclude = ["ER_G_272a", "ER_G_273a", "ER_G_292a", "ER_G_292b"]
        single_cell_obj_uker = single_cell_obj_uker[~single_cell_obj_uker.obs["orig.ident"].isin(sample_to_exclude)]
    else:
        data_center = CHUVDataCenter
        single_cell_obj_chuv = data_center.load_singlecell()
        single_cell_obj_chuv.obs["center"] = "CHUV"

single_cell_obj = single_cell_obj_uker.concatenate(single_cell_obj_chuv)

# %% -- Initialize the InferCNVAnalyzer class
infercnv_obj = InferCNVAnalyzer(adata=single_cell_obj)
# add a "reference" column to add "low_genes" and "high_counts"
# for downstream filtering
infercnv_obj.add_reference_column(cell_type_annotation_key="celltype_level2_scanvi")

# %% -- Run infercnv inference
import time
start_time = time.time()
infercnv_obj.run_infercnv(list_reference_cat=[
    "MoMac", "Microglia", 'Astrocyte', "Muscle", "T_NK",
    "Endothelial", "Oligodendrocyte", "B_cell", "Granulocyte",
    'Excitatory_Glutamatergic', "DC", "Neuron", 'Inhibitory_Cortical_GABAergic',
    'Plasma', 'Fibroblast', 'Mast'
    ])

infercnv_obj.run_leiden_clustering(resolution=0.3)

end_time = time.time()
print(f"Total time running: {round((end_time-start_time)/60)} minutes")

# %% -- Assign each cnv cluster a CNV score
infercnv_obj.compute_cnv_score()

# %% -- Perform visualizations
infercnv_obj.composite_plots(run_tl_umap=True, save_path="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/infercnv_runs/cnv_plots.png")

# %% -- Create the annotation file
infercnv_obj.tumor_clone_annotation(save_path="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/infercnv_runs/tumor_clone_annotation.csv")

# %% -- Save the summary heatmap
cnv.pl.chromosome_heatmap_summary(
    infercnv_obj.adata,
    groupby="cnv_leiden",
    save="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/infercnv_runs/cnv_heatmap.png",
    dendrogram=True
)

# %% -- Recompute genomic window and derive start-end position
infercnv_obj.get_bin_positions(
    save_path="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/infercnv_runs/genomic_windows.csv"
)
