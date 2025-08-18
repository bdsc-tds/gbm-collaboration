# %% load packages, function and data

import pandas as pd
import numpy as np

# load clinical data
from gbmcollaboration import CHUVDataCenter, UKERDataCenter

## load one center
Center = "UKER"
Project = "Glimmer" #"BoC" # or GLIMMER
date = "" #last used "2025-08-18"

if Center == "UKER":
    center = UKERDataCenter
    # preselected samples for BoC:
else:
    center = CHUVDataCenter

source_dict = center.load_tabular(exclude_idh1_mut=True)
clin_df = source_dict["clinical"]["processed gbm clinical"]

BoC_samples = [
    # training/testing samples
    "ER_G_195a", "ER_G_196a", "ER_G_197a", "ER_G_198a", "ER_G_201a",
    "ER_G_202a", "ER_G_206a", "ER_G_207a", "ER_G_307a", "ER_G_205a",
    "CH_G_152a", "CH_G_153a", "CH_G_154a",
    "CH_G_155a", "CH_G_156a", "CH_G_158a", "CH_G_159a", "CH_G_160a",
    "CH_G_161a"," CH_G_162a", "CH_G_163a", "CH_G_164a", "CH_G_165a",
    "CH_G_166a", "CH_G_167a", "CH_G_168a", "CH_G_169a", "CH_G_170a",
    "CH_G_171a", "CH_G_172a", "CH_G_173a", "CH_G_174a", "CH_G_175a",
    "CH_G_176a", "CH_G_177a", "CH_G_181a", "CH_G_182a", "CH_G_183a",
    "CH_G_184a", "CH_G_185a", "CH_G_186a", "CH_G_187a", "CH_G_191a",
    # validation samples
    'CH_G_180a', 'CH_G_189a', 'CH_G_190a', 'CH_G_191a', 'CH_G_192a',
    'CH_G_193a', 'CH_G_195a', 'CH_G_196a', 'CH_G_197a', 'CH_G_198a',
    'CH_G_199a', 'CH_G_200a', 'CH_G_201a', 'CH_G_202a', 'CH_G_203a',
    'CH_G_204a', 'CH_G_205a', 'CH_G_206a', 'CH_G_207a', 'ER_G_269a',
    'ER_G_270a', 'ER_G_276a', 'ER_G_280a', 'ER_G_284b', 'ER_G_285a',
    'ER_G_294a', 'ER_G_295a', 'ER_G_298a', 'ER_G_299a', 'ER_G_300a',
    'ER_G_306a', 'ER_G_308a', 'ER_G_310a']

if Project == "BoC":
    clin_df = clin_df[clin_df["Tumour block ID for MOSAIC"].isin(BoC_samples)]
if Project == "Glimmer":
    clin_df = clin_df[clin_df["Sample collection chronology"] == "Baseline"]


print("Shape of Clinical df is:", clin_df.shape)

# %%
from gbmcollaboration.data.utils.descriptive_stats import create_demographic_table

BASE_DATA_DIR="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/descriptive_stats"

demographic_table = create_demographic_table(
    clin_df, f"{BASE_DATA_DIR}", verbose=False
)
data_path = f"{BASE_DATA_DIR}/processed/demographic_table_{date}_{Project}_{Center}.csv"

demographic_table.to_csv(data_path)
