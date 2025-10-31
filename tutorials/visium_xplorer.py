# # Script to explore the visium data.

# %% -- Import packages and prepare functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter, AllDataCenters
from gbmcollaboration.viz.visium_heatmap import (
    plot_celltype_abun_visium,
    plot_gene_expression_visium
)

# %% -- Constants
# Constants
center = "CHUV" # or "CHUV"
normalization = "raw_count"  # Can be "raw_count", "log_norm_count", "norm_count", "SCTransform"
deconvolution_slot_name = "deconv_cohort_level_4"
# other deconvolution available: 'deconv_GBMatlas_2steps',
# 'deconv_sample_level_custom1', 'deconv_cohort_level_2',
# 'deconv_GBMatlas', 'deconv_hipe'

# %% -- Load spatial data as anndata for all centers

adata = AllDataCenters.load_spatial_adata(
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
    list_sample_id=None,
    normalization="raw_count",
    radius=1,
    deconvolution_slot_name=CELL_TYPE_KEY,
)

# %% Load spatial data for a given center
if center == "UKER":
    center_data = UKERDataCenter
else:
    center_data = CHUVDataCenter

# %% -- Load spatial data for the full cohort
# Load spatial data for the full cohort
# You can also use .load_spatial_adata
spatial_object = center_data.load_spatial(
    normalization=normalization,
    deconvolution_slot_name=deconvolution_slot_name,
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
)

# you can now access the genecounts_df attribute or celltype_df attribute as such:
spatial_counts_df = spatial_object.genecounts_df
cellabun_df = spatial_object.celltype_df
cellabun_df.head()

# the coordinates of the spots can be found as the last two columns of the dataframes

# %% -- Load spatial data for a single sample
# Load spatial data for a single sample
sample_id = "ER_G_195a_vis"
normalization = "raw_count"  # Available methods are "raw_count", "log_norm_count", "norm_count", "SCTransform"
spatial_object = center_data.load_spatial(
    normalization=normalization,
    deconvolution_slot_name=deconvolution_slot_name,
    list_sample_id=[sample_id],
)

spatial_counts_df = spatial_object.genecounts_df
cellabun_df = spatial_object.celltype_df
cellabun_df.head()

# %% -- Plotting the cell type abundance on Visium image
# Plotting the cell type abundance on Visium image
plot_celltype_abun_visium(cellabun_df, "ER_G_195a_vis", "Mural cell")

# %% -- Plotting the biomarker expression on Visium image
# Plotting the biomarker expression on Visium image
plot_gene_expression_visium(spatial_counts_df, "ER_G_195a_vis", "CD34")
# -

# ## Colocalisation of cell types

# %% -- Colocalisation of cell types
# First calculate the mean expression to potentially visualize necrosis area
spatial_counts_df["mean_counts"] = spatial_counts_df.drop(columns=["imagerow", "imagecol"]).mean(axis=1)

# Merge the mean counts with cell abundance
cellabun_meancounts_df = cellabun_df.merge(spatial_counts_df["mean_counts"], left_index=True, right_index=True)
spatial_counts_df.drop(columns="mean_counts", inplace=True)

# Reset index for ploting
cellabun_meancounts_df.reset_index(inplace=True)
cellabun_meancounts_df.set_index('sample', inplace=True)

# Define if you want to save the plots in a pdf
SAVE_PDF = False  # by default
# Define the cell types to plot along with the mean counts
CELL_TYPES = ["AC-like", "MES-like", "OPC-like", "NPC-like"]


# %% -- Plotting the colocalisation of cell types
# Function to plot colocalisation
def _plot_colocalisation(df, sample_id, cell_types):
    fig, axs = plt.subplots(1, len(cell_types), figsize=(20, 5))
    for i, celltype in enumerate(cell_types):
        plot_celltype_abun_visium(df, sample_id, celltype, ax=axs[i])
    plt.tight_layout()
    return fig

if SAVE_PDF:
    with PdfPages(f"spatial_localisation_subtypes_{center}.pdf") as pdf:
        for sample in cellabun_meancounts_df.index.unique():
            fig = _plot_colocalisation(cellabun_meancounts_df, sample, CELL_TYPES + ["mean_counts"])
            pdf.savefig(fig)
            plt.close(fig)
else:
    for sample in cellabun_meancounts_df.index.unique():
        _plot_colocalisation(cellabun_meancounts_df, sample, CELL_TYPES + ["mean_counts"])
        plt.show()
# -
