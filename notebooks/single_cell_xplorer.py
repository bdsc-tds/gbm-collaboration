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
elif center =='CHUV':
    data_center = CHUVDataCenter
else:
    raise ValueError('undefined center name')
# %% -- Load the single-cell data
single_cell_obj = data_center.load_singlecell()

# create a LogNormalize layer
single_cell_obj.layers["LogNormalize"] = single_cell_obj.layers["unnormalized_counts"].copy()
sc.pp.normalize_total(
    single_cell_obj, layer="LogNormalize", target_sum=1e4
)
sc.pp.log1p(single_cell_obj, layer="LogNormalize")

# %% Plot the UMAP with the reference color
sc.pl.umap(single_cell_obj, color="reference")
# sc.pl.umap(single_cell_obj, color="orig.ident")
