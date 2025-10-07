"""Script to explore the single-cell data."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import infercnvpy as cnv

from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter
from matplotlib.pyplot import rc_context

sc.set_figure_params(dpi=100, color_map="viridis_r")
sc.settings.verbosity = 0

# %% -- Load data
center = "CR" # or "CHUV"

if center == "UKER":
    data_center = UKERDataCenter
elif center == "CR":
    data_center = CRDataCenter
else:
    data_center = CHUVDataCenter

# %% -- Load the single-cell data
single_cell_obj = data_center.load_singlecell()

# %%
# run infercnv
path_gtf = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/infercnv/gencode.v32.annotation.gtf"
cnv.io.genomic_position_from_gtf(path_gtf, single_cell_obj)

# create a "reference" column in the obs to extract the healthy reference cells
mapping_reference = {
    'Malignant_gbm': "malignant",
    'Endothelial': "endothelial_mural",
    'Smooth_muscle': "endothelial_mural",
    'MoMac': "healthy",
    'T_NK': "healthy",
    'Oligodendrocyte' : "healthy",
    'Granulocyte': "healthy",
    'Neuron': "healthy"
}
single_cell_obj = add_reference_column(
    single_cell_obj,
    cell_type_annotation_key="celltype_level2_scanvi",
    celltypes_to_reference=mapping_reference,
    threshold_n_genes=400,
    threshold_n_counts=10000
)
# create a LogNormalize layer
single_cell_obj.layers["LogNormalize"] = single_cell_obj.layers["unnormalized_counts"].copy()
sc.pp.normalize_total(
    single_cell_obj, layer="LogNormalize", target_sum=1e4
)
sc.pp.log1p(single_cell_obj, layer="LogNormalize")

# %% Plot the UMAP with the reference color
sc.pl.umap(single_cell_obj, color="reference")
# sc.pl.umap(single_cell_obj, color="orig.ident")
