"""Script to check the IDH and H3 mutations status across samples."""

# %% Load packages and classes
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import openslide
from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter, AllDataCenters

# # Access data for all the cohort
all_centers = AllDataCenters

source_dict = all_centers.load_tabular(
    exclude_idh1_mut=True,
    exclude_patient_ids=True,
    exclude_relapse=True)

# H3F3A # ENSG00000163041 # but couldn't find this gene
# HISTH3B/C > H3C2 and H3C3 > ENSG00000286522 and ENSG00000287080
# could only find H3C2

# %% Load WES data

centers = ["mosaic_chuv", "mosaic_uker", "mosaic_cr"]

for center in centers:
    df_wes = source_dict[center]["wes"]["WES mutations"]
    for gene in ["H3C2", "IDH1", "IDH2"]:
        print(f"For the {center} center, there is")
        print(df_wes[gene].sum())
        print(f"for gene {gene}")

# H3F3A # ENSG00000163041 # but couldn't find this gene
# HISTH3B/C > H3C2 and H3C3 > ENSG00000286522 and ENSG00000287080
# could only find H3C2

# %%
# for IDH2, there is one mutant:
df_wes_uker = source_dict["mosaic_uker"]["wes"]["WES mutations"]

df_wes_uker_idh2 = df_wes_uker[df_wes_uker["IDH2"]].index
df_wes_uker_idh2
# ER_G_619a has an IDH2 mutation (R140Q) but it's not part of the canonical glioma IDH spectrum

# %%
