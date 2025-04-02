"""Script to compare deconvolution methods."""

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go

from gbmcollaboration import UKERDataCenter, CHUVDataCenter
from gbmcollaboration.utils.celltype_discordance_fcts import (
    process_spatial_data,
    discordant_deconvolution,
    dominant_celltype_frac,
    violin_plot_celltype_abundance,
    pdf_violin_plot_across_cohort,
)

# print full table
pd.set_option('display.max_columns', None)

# %% --- Loading data ---
# Start with the spatial object from CHUV and UKER
uker_spatial_object_step1 = UKERDataCenter.load_spatial(
    normalization="raw_count",
    deconvolution_slot_name="deconv_GMBatlas",
    exclude_low_quality_visium=True)
uker_spatial_object_step2 = UKERDataCenter.load_spatial(
    normalization="raw_count",
    deconvolution_slot_name="deconv_GMBatlas_2steps",
    exclude_low_quality_visium=True)

chuv_spatial_object_step1 = CHUVDataCenter.load_spatial(
    normalization="raw_count",
    deconvolution_slot_name="deconv_GMBatlas",
    exclude_low_quality_visium=True)
chuv_spatial_object_step2 = CHUVDataCenter.load_spatial(
    normalization="raw_count",
    deconvolution_slot_name="deconv_GMBatlas_2steps",
    exclude_low_quality_visium=True)

# %%
# Example usage for UKER
discordant_samples_uker, dominant_celltype_freq_uker, concat_abun_uker_step1, concat_abun_uker_step2 = process_spatial_data(
    spatial_object_step1=uker_spatial_object_step1,
    spatial_object_step2=uker_spatial_object_step2
)

# %%
# Example usage for CHUV
discordant_samples_chuv, dominant_celltype_freq_chuv, concat_abun_chuv_step1, concat_abun_chuv_step2 = process_spatial_data(
    spatial_object_step1=chuv_spatial_object_step1,
    spatial_object_step2=chuv_spatial_object_step2
)

# %% -- Average discordance between the two deconvolution methods --
pd.DataFrame(discordant_samples_chuv.mean())

# %% --- Plotting cell abundance for one sample ---
# eg 'ER_G_195a_vis'
sample = 'ER_G_195a_vis'
violin_plot_celltype_abundance(concat_abun_uker_step2, sample)

# %% --- Plotting cell abundance for all samples and save a pdf ---
# and save it in the current directory
pdf_violin_plot_across_cohort(concat_abun_uker_step2, 'violin_plt_celltype_abundance_uker.pdf')

# %% --- Plotting cell abundance for chuv ---
pdf_violin_plot_across_cohort(concat_abun_chuv_step2, 'violin_plt_celltype_abundance_chuv.pdf')

# %% --- Summary table for CHUV and several samples
dominant_celltype_frac(concat_abun_chuv_step2, chuv_spatial_object_step2.common_celltypes, freq=True)

# %% --- Summary table for CHUV with cell type abundance, not fraction
cell_abun_chuv = dominant_celltype_frac(concat_abun_chuv_step2, chuv_spatial_object_step2.common_celltypes, freq=False)
cell_abun_chuv = cell_abun_chuv.drop(columns=["dominant_celltype"]).reset_index().drop(columns="spot").groupby("sample").agg("mean")
# order we want:
order = ["AC-like","MES-like","NPC-like","OPC-like","Astrocyte","Neuron","OPC","Oligodendrocyte","RG","Endothelial","Mural cell","Mast","Mono","NK","DC","B cell","CD4/CD8","Neutrophil","Plasma B","TAM-BDM","TAM-MG"]
cell_abun_chuv[order]

# %% --- Summary table for UKER and several samples
dominant_celltype_frac(concat_abun_uker_step2, uker_spatial_object_step2.common_celltypes, freq=True)

# %% --- Summary table for UKER with cell type abundance, not fraction
cell_abun_uker = dominant_celltype_frac(concat_abun_uker_step2, uker_spatial_object_step2.common_celltypes, freq=False)
cell_abun_uker = cell_abun_uker.drop(columns=["dominant_celltype"]).reset_index().drop(columns="spot").groupby("sample").agg("mean")
# order we want:
order = ["AC-like","MES-like","NPC-like","OPC-like","Astrocyte","Neuron","OPC","Oligodendrocyte","RG","Endothelial","Mural cell","Mast","Mono","NK","DC","B cell","CD4/CD8","Neutrophil","Plasma B","TAM-BDM","TAM-MG"]
cell_abun_uker[order]
