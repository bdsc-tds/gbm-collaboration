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

# %% -- Run inferCNV with/without layer option
def prepare_and_run_infercnv(adata, list_reference_cat):
    infercnv = InferCNVAnalyzer(adata=adata)
    infercnv.add_reference_column(cell_type_annotation_key="celltype_level2_scanvi")
    infercnv.run_infercnv(list_reference_cat=list_reference_cat)
    infercnv.run_leiden_clustering(resolution=0.4)
    infercnv.compute_cnv_score()
    return infercnv

# make subsample
adata = sc.pp.subsample(single_cell_obj, n_obs=10_000, copy=True)
list_reference_cat = adata[
    ~adata.obs["celltype_level1_scanvi"].isin(["Malignant", "Doublet", "Other", "Unknown"])
].obs["celltype_level2_scanvi"].unique()

# copy to test if results change across runs
adata2 = adata.copy()
# copy to test if results change with layer option
adata3 = adata.copy()

# Run inferCNV
infercnv_obj1 = prepare_and_run_infercnv(adata, list_reference_cat)
infercnv_obj2 = prepare_and_run_infercnv(adata2, list_reference_cat)

adata3.X = adata.layers['LogNormalize']
infercnv_obj3 = prepare_and_run_infercnv(adata3, list_reference_cat)


# %%
assert np.all(infercnv_obj1.adata.obs['cnv_score']==infercnv_obj2.adata.obs['cnv_score']), 'Different results found across infercnv runs'
assert np.all(infercnv_obj1.adata.obs['cnv_score']==infercnv_obj3.adata.obs['cnv_score']), 'Different results found with infercnv layer parameter'
