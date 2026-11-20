"""Table 1b-1c: demographic and clinical description of the cohort, per centre and overall.

Origin: analysis/cohort_desc/cohort_description.py.
Output: demographic_table.csv (273 baseline patients: UKER 135, CHUV 65, Charite 73),
identical to the published descriptive_stats/processed/demographic_table_2026-06-10_Glimmer.csv.
"""

# %% load packages, function and data
import pandas as pd
import numpy as np

# load clinical data
from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter, AllDataCenters
from gbmcollaboration.data.utils.descriptive_stats import create_demographic_table
from gbmcollaboration.paper import output

# constants
OUT = output.setup("table1bc_cohort_description")
BASE_DATA_DIR = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/descriptive_stats"
ADD_OVERALL_ROW = True  # set False if you don't want an overall summary

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

# load all centers
centers = [
    ("UKER", UKERDataCenter),
    ("CHUV", CHUVDataCenter),
    ("CR", CRDataCenter)
]

def _load_center_df(center_cls):
    """Load and filter a center's clinical dataframe according to your rules."""
    source_dict = center_cls.load_tabular(
        exclude_idh1_mut=True,
        exclude_patient_ids=True,
        exclude_relapse=True)
    df = source_dict["clinical"]["processed gbm clinical"]

    df = df[df["Sample collection chronology"] == "Baseline"]
    # drop the duplicated (low tumor content) samples if found
    exclude_sample = "ER_G_290a"
    if "Tumour block ID for MOSAIC" == df.index.name:
        col_ = df.index
    else:
        col_ = df["Tumour block ID for MOSAIC"]
    df = df[col_ != exclude_sample]

    # drop outlier values for one sample:
    # need to find back this sample before excluding
    # df = df[df["GB - Baseline hemoglobin level (g/dL)"] < 1000]

    return df

# %% build a single demographic table with center under each feature
tables = []
for center_name, center_cls in centers:
    clin_df = _load_center_df(center_cls)
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
        pooled.append(_load_center_df(center_cls))
    pooled_df = pd.concat(pooled, ignore_index=True)

    overall = create_demographic_table(pooled_df, f"{BASE_DATA_DIR}", verbose=False)
    overall.insert(0, "Center", "All Centers")

    combined = pd.concat([combined, overall], ignore_index=True)

# group by variable then center for easy reading in publications
combined = combined.sort_values(["Variable", "Center"], kind="stable").reset_index(drop=True)

# drop features
combined = combined[~combined["Variable"].isin(excluded_features)]

# saving
output.save_table(OUT, "demographic_table", combined, panel="Table 1b-c", index=False)
print("Saved:", OUT / "demographic_table.csv")

# %%
