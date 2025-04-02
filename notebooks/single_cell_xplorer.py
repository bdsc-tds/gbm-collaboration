"""Script to explore the single-cell data."""

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import infercnvpy as cnv

from gbmcollaboration import CHUVDataCenter, UKERDataCenter
from gbmcollaboration.utils.infercnv_utils import add_reference_column

from matplotlib.pyplot import rc_context

sc.set_figure_params(dpi=100, color_map="viridis_r")
sc.settings.verbosity = 0
sc.logging.print_header()

# %% -- Load data
center = "CHUV" # or "CHUV"

if center == "UKER":
    data_center = UKERDataCenter
else:
    data_center = CHUVDataCenter

# %% -- Load the single-cell data
single_cell_obj = data_center.load_singlecell()

# %%
# run infercnv
# most recent version of gencode gtf file was saved in the omics abstra dataset,
# other versions can be found at https://genome.ucsc.edu or https://www.gencodegenes.org/
path_gtf = "/home/owkin/data/dataset-5abe643c-1c3a-4117-b5f5-35783a73b3c7/infercnv/gencode.v47.annotation.gtf"
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

# %% run inferCNV
cnv.tl.infercnv(
    single_cell_obj,
    reference_key="reference",
    reference_cat=["healthy"],
    window_size=200,
    layer="LogNormalize"
)

# %%
cnv.pl.chromosome_heatmap(single_cell_obj, groupby="reference")
cnv.tl.pca(single_cell_obj)
cnv.pp.neighbors(single_cell_obj)
cnv.tl.leiden(single_cell_obj)

# %%
# rerun the following if you changed the leiden parameters, before
# running the next cell for chromosome heatmap
sc.tl.dendrogram(single_cell_obj, groupby="cnv_leiden")

# %%
cnv.pl.chromosome_heatmap(single_cell_obj, groupby=["cnv_leiden"], dendrogram=True)
cnv.pl.chromosome_heatmap_summary(single_cell_obj, groupby="cnv_leiden", dendrogram=True)

# %%
cnv.tl.umap(single_cell_obj)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    single_cell_obj,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(single_cell_obj, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(single_cell_obj, color="cnv_leiden", ax=ax3)
