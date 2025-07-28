"""Script to infer CNV using inferCNVpy."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import infercnvpy as cnv

from gbmcollaboration import CHUVDataCenter, UKERDataCenter
from gbmcollaboration import PATH_GTF

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
def prepare_and_run_infercnv(adata, list_reference_cat, layer=None):
    #infercnv
    cnv.io.genomic_position_from_gtf(PATH_GTF, adata)
    cnv.tl.infercnv(
        adata,
        reference_key="celltype_level2_scanvi",
        reference_cat=list_reference_cat,
        layer=layer
        )
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata, resolution=0.4)
    cnv.tl.cnv_score(adata)
    return adata

# make subsample
adata = sc.pp.subsample(single_cell_obj, n_obs=10_000, copy=True)

# ensure raw counts are in adata.X
adata.X = adata.layers['unnormalized_counts']

# lognorm
adata.layers["LogNormalize"] = adata.layers["unnormalized_counts"].copy()
sc.pp.normalize_total(adata, layer="LogNormalize", target_sum=1e4)
sc.pp.log1p(adata, layer="LogNormalize")

# copy to test if results change across runs
adata2 = adata.copy()
# copy to test if results change with layer option
adata3 = adata.copy()
adata3.X = adata.layers['LogNormalize']


# define reference cats
list_reference_cat = adata[
    ~adata.obs["celltype_level1_scanvi"].isin(["Malignant", "Doublet", "Other", "Unknown"])
].obs["celltype_level2_scanvi"].unique()

# Run inferCNV
infercnv_obj1 = prepare_and_run_infercnv(adata, list_reference_cat,layer="LogNormalize")
infercnv_obj2 = prepare_and_run_infercnv(adata2, list_reference_cat,layer="LogNormalize")
infercnv_obj3 = prepare_and_run_infercnv(adata3, list_reference_cat,layer=None)


# %%
assert np.all(adata.obs['cnv_score']==adata2.obs['cnv_score']), 'Different results found across infercnv runs'
assert np.all(adata.obs['cnv_score']==adata3.obs['cnv_score']), 'Different results found with infercnv layer parameter'

# %%
