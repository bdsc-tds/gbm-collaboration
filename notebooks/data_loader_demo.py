# ## Notebook to explore and visualize the MOSAIC data modalities

# ## Data loaders

# The current implementation allows to rapidly visualize the data sources available for each cohort (ie. CHUV and UKER).

# + Load packages and classes
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from openslide import open_slide

from gbmcollaboration import CHUVDataCenter, UKERDataCenter

# + Available constants that can be used to filter the data
from gbmcollaboration import (
    PATIENTS_IDH1_MUT,
    BAD_QUALITY_UKER_VISIUM,
    BAD_QUALITY_CHUV_VISIUM
)
# -

# # Access data for the given cohort
# This returns a dictionary with the data sources in the keys and the list of files (and path) used to return the data.

# Choose your center of interest, here it is CHUV
chuv_center = CHUVDataCenter
chuv_center.sources.keys()

# example with UKER
uker_center = UKERDataCenter
uker_center.sources.keys()

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
source_dict = chuv_center.load_tabular(exclude_idh1_mut=True)

# H&E slides
he = source_dict["he"]["HE files"]
he.head(2)

# +
# Display a H&E slide

slide_idx = 12  # pick a H&E slide to display

slide = open_slide(he.path.values[slide_idx])

# Display the slide
slide_img = slide.get_thumbnail(slide.level_dimensions[-2])
fig, ax = plt.subplots()
ax.imshow(slide_img)
ax.set_axis_off()
plt.show()
# -

# H&E AI features
he_fea = source_dict["he"]["H0 features"]
he_fea.head(2)

# Clinical data
clin_df = source_dict["clinical"]["processed gbm clinical"]
clin_df.head(2)

uker_center.load_tabular(exclude_idh1_mut=True)

# + Example with wes loading
CNVamp = source_dict["wes"]["WES CNV amplification"]
CNVamp.head(2)
# -

# ### Access and load Visium spatial data for the given cohort

# +
chuv_spatial_object = chuv_center.load_spatial(
    normalization="raw_count",
    deconvolution_slot_name="deconv_GMBatlas_2steps",
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True)

chuv_spatial_counts = chuv_spatial_object.load_counts_cohort()

# This returns a dictionary with the samples ID as keys and the counts in each spot as values
# -

# Information about the barcodes and spatial coordinates
chuv_spatial_object.coords

# + load the cell type abundance
# if you wanted cell fraction instead, use the argument frac=True
chuv_celltype_abund = chuv_spatial_object.load_deconvolution_cohort()
common_celltypes = chuv_spatial_object.common_celltypes
