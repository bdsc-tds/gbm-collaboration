# ## Notebook to explore and visualize the MOSAIC data modalities

# ## Data loaders

# The current implementation allows to rapidly visualize the data sources available for each cohort (ie. CHUV and UKER).

# %% Load packages and classes
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import openslide
from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter, AllDataCenters

# %% Available constants that can be used to filter the data
from gbmcollaboration import (
    PATIENTS_IDH1_MUT,
    RELAPSE_EXCLUSION,
    PATIENT_EXCLUSION,
    BAD_QUALITY_UKER_VISIUM,
    BAD_QUALITY_CHUV_VISIUM,
    BAD_QUALITY_CR_VISIUM
)

# # Access data for all the cohort
all_centers = AllDataCenters

# Choose your center of interest, here it is CHUV
# This returns a dictionary with the data sources in the keys and the list of files (and path) used to return the data.
chuv_center = CHUVDataCenter
chuv_center.sources.keys()

# example with UKER
uker_center = UKERDataCenter
uker_center.sources.keys()

# example with CT
cr_center = CRDataCenter
cr_center.sources.keys()

# # How to load the tabular and clinical data in memory
# Available tabular data are:
# - data_dict: dictionary created by Owkin to keep tracks of available features
#     and their data types
# - clinical:
#     - original clinical: raw clinical data that have been processed by mosaic pipeline
#     - key events: summary of the date clinical data processed by mosaic pipeline
#     - processed gbm clinical: further processing in this repository
#     - endpoints: outcome data processed in this repository
# - bulk_rna:
#     - TPM: Transcripts Per Million
#     - normalized counts: DESeq2 normalization
# - wes:
#     - WES mutations: snv and indels
#     - WES CNV deletion
#     - WES CNV amplification
# - he:
#     - HE files: path of the H&E slides
#     - H-Bioptimus-0 features: path to the tiles features

# Load the tabular data and clinical data
# Note that you can exclude patients with IDH1 mutation by setting exclude_idh1_mut=True
source_dict = cr_center.load_tabular(
    exclude_idh1_mut=True,
    exclude_patient_ids=True,
    exclude_relapse=True)

# %%
# Clinical data
clin_df = source_dict["clinical"]["processed gbm clinical"]
clin_df.head(2)

# %%
# Same thing for UKER data
uker_center.load_tabular(exclude_idh1_mut=True)

# %% Example with wes loading
CNVamp = source_dict["wes"]["WES CNV amplification"]
CNVamp.head(2)

# %% Or the loading of the TMB values
TMB = source_dict["wes"]["WES TMB"]
TMB.head(2)

# %% Access and load Visium spatial data for the given cohort
cr_spatial_object = cr_center.load_spatial(
    normalization="raw_count",
    deconvolution_slot_name="deconv_cohort_level_4",
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
    exclude_patient_ids=True,
    exclude_relapse=True)
# cr_spatial_object is a dictionary with the samples ID as keys and the counts in each spot as values

# Obtain the counts
cr_spatial_counts = cr_spatial_object.load_counts_cohort()

# Information about the barcodes and spatial coordinates
cr_spatial_object.coords

# %% load the cell type abundance
# if you wanted cell fraction instead, use the argument frac=True
cr_celltype_abund = cr_spatial_object.load_deconvolution_cohort()
common_celltypes = cr_spatial_object.common_celltypes

# %% Load the single-cell data
single_cell_obj = cr_center.load_singlecell(
    exclude_idh1_mut = True,
    exclude_patient_ids = True,
    exclude_relapse = True
)

# look at the single_cell_xplorer.py script to explore more!

# %%
