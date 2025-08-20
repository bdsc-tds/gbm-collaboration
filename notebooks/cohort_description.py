# %% load packages, function and data
import pandas as pd
import numpy as np

# load clinical data
from gbmcollaboration import CHUVDataCenter, UKERDataCenter
from gbmcollaboration.data.utils.descriptive_stats import create_demographic_table

# constants
Project = "BoC"  # or "BoC"
date = "2025-08-19"  # e.g., "2025-08-18" if you want it in the filename
BASE_DATA_DIR = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/descriptive_stats"
ADD_OVERALL_ROW = True  # set False if you don't want an overall summary

# BoC sample list
BoC_samples = [
    # training/testing samples
    "ER_G_195a", "ER_G_196a", "ER_G_197a", "ER_G_198a", "ER_G_201a",
    "ER_G_202a", "ER_G_206a", "ER_G_207a", "ER_G_307a", "ER_G_205a",
    "CH_G_151a", "CH_G_152a", "CH_G_153a", "CH_G_154a",
    "CH_G_155a", "CH_G_156a", "CH_G_158a", "CH_G_159a", "CH_G_160a",
    "CH_G_161a", "CH_G_162a", "CH_G_163a", "CH_G_164a", "CH_G_165a",
    "CH_G_166a", "CH_G_167a", "CH_G_168a", "CH_G_169a", "CH_G_170a",
    "CH_G_171a", "CH_G_172a", "CH_G_173a", "CH_G_174a", "CH_G_175a",
    "CH_G_176a", "CH_G_177a", "CH_G_181a", "CH_G_182a", "CH_G_183a",
    "CH_G_184a", "CH_G_185a", "CH_G_186a", "CH_G_187a", "CH_G_191a",
    # validation samples
    'CH_G_180a', 'CH_G_189a', 'CH_G_190a', 'CH_G_192a',
    'CH_G_193a', 'CH_G_195a', 'CH_G_196a', 'CH_G_197a', 'CH_G_198a',
    'CH_G_199a', 'CH_G_200a', 'CH_G_201a', 'CH_G_202a', 'CH_G_203a',
    'CH_G_204a', 'CH_G_205a', 'CH_G_206a', 'CH_G_207a', 'ER_G_269a',
    'ER_G_270a', 'ER_G_276a', 'ER_G_280a', 'ER_G_284b', 'ER_G_285a',
    'ER_G_294a', 'ER_G_295a', 'ER_G_298a', 'ER_G_299a', 'ER_G_300a',
    'ER_G_306a', 'ER_G_308a', 'ER_G_310a'
]

# features to exclude from clinical df:
excluded_features = [
    "Ethnicity", "Specify minority - White", "EGFR amplification",
    "Glioblastoma - mutational status", "Chromosome status", "TERT promoter mutation",
    "GB - BRAF mutation", "Glioblastoma - other mutation?", "Glioblastoma - other mutation",
    "OLIG2", "p16 / CDKN2A (IHC)", "GB - EGFR mutation", "ATRX mutation", "BRAF V600E expression (IHC)",
    "GB - PTEN mutation", "FGFR3-TACC3 fusion", "NTRK fusion", "Histone H3 K27M mutation (sequencing)",
    "Histone H3 G34 mutation (sequencing)", "Baseline whole blood count performed 98 - Code.1",
    "GB - Baseline neutrophils count (/mm3)",
    "Baseline lymphocyte count (/mm3)", "Alanine transaminase (UI/L)", "Aspartate aminotransferase (UI/L)",
    "Alkaline phosphatase (UI/L)", "Gamma-glutamyltransferase (UI/L)", "Total bilirubin (mg/dL)",
    "Unconjugated bilirubin (mg/dL)", "Specify history of other cancer",
    "Has there been any progression(s)/recurrence(s) before MOSAIC?", "Subject history of auto-immune diseases",
    "Specify auto-immune disease(s)", "Subject history of immunosuppression 0", "Hematologic disorders", "Diabetes",
    "Smoking status", "Alcohol intake", "Other medical / surgical history", "Other medical conditions: specify",
    "clinical_trial_participation", "Participation in clinical trial(s): specify", "previous_cancer",
    "Other hematologic disorders: specify", "Metabolic syndrome", "Family history of cancer (first degree)",
    "Family history of cancer: type 99", "Other type of cancer", "Cause of death",
    "Performance status (ECOG) at cancer diagnosis", "Weight at follow up (kg)", "Height at follow up (cm)",
    "BMI at follow up (kg/m2)", "Vital status at inclusion", "violates_cohort_inclusion_criteria_cohort1",
    "cohort_violations", "Surgery adjunct 1", "Surgery adjunct 6", "concurrent_treatments",
    "TLS observed from same date / same site", "metastatic_status_ok", "OS", "PFS",
    "cancer_death", "cancer_death_with_progression_dates", "death_unrelated_to_cancer",
    "Sample collection chronology", "sample_classification",
    "Systemic therapy implant", "os_years", "pfs_years",
]

# load both centers
centers = [
    ("UKER", UKERDataCenter),
    ("CHUV", CHUVDataCenter),
]

def _load_center_df(center_cls, project):
    """Load and filter a center's clinical dataframe according to your rules."""
    source_dict = center_cls.load_tabular(exclude_idh1_mut=True)
    df = source_dict["clinical"]["processed gbm clinical"]

    if project == "BoC":
        df = df[df["Tumour block ID for MOSAIC"].isin(BoC_samples)]
    elif project == "Glimmer":
        df = df[df["Sample collection chronology"] == "Baseline"]
        # drop the duplicated (low tumor content) samples if found
        exclude_sample = "ER_G_290a"
        df = df[df["Tumour block ID for MOSAIC"]!= exclude_sample]

    return df

# %% build a single demographic table with center under each feature
tables = []
for center_name, center_cls in centers:
    clin_df = _load_center_df(center_cls, Project)
    print(f"[{center_name}] Clinical df shape:", clin_df.shape)

    demo = create_demographic_table(clin_df, f"{BASE_DATA_DIR}", verbose=False)
    demo.insert(0, "Center", center_name)  # make Center the first column
    tables.append(demo)

combined = pd.concat(tables, ignore_index=True)

# add a total row per feature by recomputing on the stacked data
if ADD_OVERALL_ROW:
    # Create a pooled baseline dataframe across centers (same filters already applied)
    pooled = []
    for _, center_cls in centers:
        pooled.append(_load_center_df(center_cls, Project))
    pooled_df = pd.concat(pooled, ignore_index=True)

    overall = create_demographic_table(pooled_df, f"{BASE_DATA_DIR}", verbose=False)
    overall.insert(0, "Center", "Both Centers")

    combined = pd.concat([combined, overall], ignore_index=True)

# group by variable then center for easy reading in publications
combined = combined.sort_values(["Variable", "Center"], kind="stable").reset_index(drop=True)

# drop features
combined = combined[~combined["Variable"].isin(excluded_features)]

# saving
data_path = f"{BASE_DATA_DIR}/processed/demographic_table_{date}_{Project}.csv"
combined.to_csv(data_path, index=False)
print("Saved:", data_path)

# %%
