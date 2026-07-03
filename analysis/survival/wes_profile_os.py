"""
Survival analysis between CNA clusters,
using covariates. It outputs the figure 3d of the manuscript.
"""

# %% -- Import packages and load data
import pandas as pd
import numpy as np
from pathlib import Path
from sksurv.ensemble import RandomSurvivalForest
from sklearn.ensemble import RandomForestRegressor
from sksurv.util import Surv
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from lifelines.utils import concordance_index
# import shap
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from statsmodels.stats.multitest import fdrcorrection

from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter
from gbmcollaboration import (
    PATIENT_EXCLUSION,
    RELAPSE_EXCLUSION,
    PATIENTS_IDH1_MUT,
    SURVIVAL_EXCLUSION
)
from gbmcollaboration.definitions import AVG_NUM_DAYS_PER_YEAR
from gbmcollaboration.utils.clinical_functions import (
    univariate,
    plot_forest,
    plot_correlations
)
from gbmcollaboration.utils.viz_survival import kmplot
from gbmcollaboration.viz import set_mpl_rcparams, savefig

set_mpl_rcparams()

FIGURES_DIR = Path('/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/')

# -- Loop over the three data centers and merge the clinical data
merged_df = []

for center in [CHUVDataCenter, UKERDataCenter, CRDataCenter]:
    center_data_dict = center.load_tabular(
        exclude_idh1_mut=True,
        exclude_patient_ids=True,
        exclude_relapse=True
    )
    clin_dict = center_data_dict['clinical']
    df = clin_dict['processed gbm clinical']
    # Exclude CR samples for which survival data is pending
    if center == CRDataCenter:
        df = df.drop(index=SURVIVAL_EXCLUSION)
    # To be sure check if we only have Baseline samples
    print(f"---Loading the clinical data from {center.name}---")
    print(df['Sample collection chronology'].value_counts())
    print(f"Found {len(df)} primary cases")
    # append to the merged dataframe
    merged_df.append(df)

# merge the three dataframes
df_concat = pd.concat(merged_df, ignore_index=False)
print("Dataframes from the three centers have been merged.")

################################################################
################ Survival analysis of CNA cluster ################
################################################################

# %% --- Load Leiden clusters (autosomes-only) previously computed with
# generate_cnv_wes_matrix.py and infercnv_wes_analyses.py
df_wes_leiden_cluster = pd.read_csv(
    "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
)
df_wes_leiden_cluster.set_index("orig.ident", inplace=True)

# -- Prepare dataframe
cluster_col = "cnv_leiden"
covariate_cols = ["Age (years)", "MGMT promoter methylation"]#, "Administrative gender", "Site code"]

# Pull OS + covariates from merged clinical df
df_surv_cov = df_concat[["OS"] + covariate_cols].copy()

# Map categorical covariates
MAPPING = {
    "MGMT promoter methylation": {"Unmethylated": 0, "Methylated": 1},
    "Administrative gender": {"Male": 0, "Female": 1},
    "Site code": {"CH": 0, "ER": 1, "CR": 2}
}
df_surv_cov = df_surv_cov.replace(MAPPING)

# Merge clusters with survival/covariates
df_cluster = df_wes_leiden_cluster[[cluster_col]].merge(
    df_surv_cov,
    how="inner",
    left_index=True,
    right_index=True
).copy()

df_cluster[cluster_col] = df_cluster[cluster_col].astype(int)

print("Merged df shape:", df_cluster.shape)
print("Cluster counts:\n", df_cluster[cluster_col].value_counts().sort_index())

# Build time/event in lifelines format
df_cluster["OS.time"] = df_cluster["OS"].abs()
df_cluster["OS.event"] = df_cluster["OS"] > 0

# %% -- Adjusted Cox: each cluster vs all others (one model per cluster)
clusters = sorted(df_cluster[cluster_col].unique())
df_cluster_vs_rest_rows = []

for cl in clusters:
    tmp = df_cluster.copy()
    tmp["is_cluster"] = (tmp[cluster_col] == cl).astype(int)

    use_cols = ["OS.time", "OS.event", "is_cluster"] + covariate_cols
    tmp = tmp[use_cols].dropna()

    n_in = int(tmp["is_cluster"].sum())
    n_out = int((1 - tmp["is_cluster"]).sum())

    cph = CoxPHFitter()
    cph.fit(tmp, duration_col="OS.time", event_col="OS.event", show_progress=False)

    # extract the adjusted effect for "is_cluster"
    s = cph.summary.loc["is_cluster"]
    df_cluster_vs_rest_rows.append({
        "cluster": cl,
        "n_in": n_in,
        "n_out": n_out,
        "HR": float(s["exp(coef)"]),
        "CI95_low": float(s["exp(coef) lower 95%"]),
        "CI95_high": float(s["exp(coef) upper 95%"]),
        "p": float(s["p"])
    })

    # forest plot for this adjusted model (cluster term + covariates)
    df_sum = cph.summary.reset_index().rename(columns={"index": "covariate"})
    _, df_sum["fdr"] = fdrcorrection(df_sum["p"].values)

    plot_forest(
        df_all=df_sum,
        distances=[0.45, 0.15, 0.35],
        figsize=(10, 3.8),
        outcome="OS",
        significance="FDR"  # or "FDR"
    )
    plt.title(f"Adjusted Cox: cluster {cl} vs all others")
    plt.show()

df_cox_cluster_vs_rest_adj = pd.DataFrame(df_cluster_vs_rest_rows).sort_values("p")
print("\nAdjusted Cox (cluster vs rest) summary table:")
print(df_cox_cluster_vs_rest_adj)

#######################################################
#################### FIGURE 3d ########################
#######################################################
# %% -- Adjusted categorical Cox (fixed reference cluster)
# Choose reference cluster 0 since it's the largest
ref_cluster = 0 if 0 in clusters else int(df_cluster[cluster_col].value_counts().idxmax())
print(f"\nReference cluster for categorical Cox: {ref_cluster}")

df_cat = df_cluster[[cluster_col, "OS.time", "OS.event"] + covariate_cols].dropna().copy()
df_cat[cluster_col] = pd.Categorical(df_cat[cluster_col])

# reorder categories so reference is first
new_cats = [ref_cluster] + [c for c in df_cat[cluster_col].cat.categories if c != ref_cluster]
df_cat[cluster_col] = df_cat[cluster_col].cat.reorder_categories(new_cats)

df_cat_dum = pd.get_dummies(df_cat, columns=[cluster_col], drop_first=True)

cph_cat = CoxPHFitter()
cph_cat.fit(df_cat_dum, duration_col="OS.time", event_col="OS.event", show_progress=False)

df_sum_cat = cph_cat.summary.reset_index().rename(columns={"index": "covariate"})
_, df_sum_cat["fdr"] = fdrcorrection(df_sum_cat["p"].values)

print("\nAdjusted categorical Cox summary (top hits):")
print(df_sum_cat.sort_values("p").head(20))

# Forest plot for categorical model (clusters + covariates)
# Increase global font sizes
mpl.rcParams.update({
    "font.size": 20,
    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})
df_sum_cat["covariate"] = df_sum_cat["covariate"].map(
    {
        "Age (years)" : "Age (years)",
        "Administrative gender" : "Administrative gender",
        "MGMT promoter methylation" : "MGMT promoter methylation",
        "Site code" : "Site code",
        "cnv_leiden_1" : "CNA cluster 1",
        "cnv_leiden_2" : "CNA cluster 2",
        "cnv_leiden_3" : "CNA cluster 3",
        }
    )
# Forest plot
plt = plot_forest(
    df_all=df_sum_cat,
    distances=[0.60, 0.10, 0.50],
    figsize=(10, 8),
    outcome="OS",
    significance="FDR"
)
plt.title(f"Adjusted categorical Cox (reference cluster : {ref_cluster})")
savefig(
    df = df_sum_cat,
    out_file=FIGURES_DIR/"figure3/fig3d.pdf",
)
plt.show()

# Optional: extract only cluster terms into a clean table
cluster_prefix = f"{cluster_col}_"
df_cat_clusters = df_sum_cat[df_sum_cat["covariate"].str.startswith(cluster_prefix)].copy()
df_cat_clusters["cluster"] = df_cat_clusters["covariate"].str.replace(cluster_prefix, "", regex=False)
df_cat_clusters["reference"] = ref_cluster
df_cat_clusters = df_cat_clusters[["cluster", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p", "reference"]]
df_cat_clusters = df_cat_clusters.rename(columns={
    "exp(coef)": "HR",
    "exp(coef) lower 95%": "CI95_low",
    "exp(coef) upper 95%": "CI95_high"
}).sort_values("p")

print("\nAdjusted categorical Cox (cluster terms only):")
print(df_cat_clusters)

# %% -- KM curves
df_cluster['OS.time.days'] = df_cluster['OS.time'] * AVG_NUM_DAYS_PER_YEAR
df_km = df_cluster[["OS.time.days", cluster_col]].copy()
df_km[cluster_col] = df_km[cluster_col].astype(str)

fig, ax = plt.subplots(1, 1, figsize=(9, 7))
kmplot(
    data=df_km,
    x="OS.time.days",
    hue=cluster_col,
    hue_order=sorted(df_km[cluster_col].unique(), key=lambda x: int(x)),
    convert_days_to_months=True,
    add_hazard_ratio=False,      # only works for 2 groups
    add_logrank_pvalue=False,    # global logrank isn't supported in that helper
    at_risk_counts=False,
    ax=ax,
    title="Overall survival by WES CNV Leiden cluster (autosomes only)"
)
plt.tight_layout()
savefig(
    df = df_km,
    out_file=FIGURES_DIR/"figure3/fig3c_km.pdf",
)
plt.show()

# %% KM curves with multivariate_logrank_test
fig, ax = plt.subplots(1, 1, figsize=(9, 7))
kmplot(
    data=df_km,
    x="OS.time.days",
    hue=cluster_col,
    hue_order=sorted(df_km[cluster_col].unique(), key=lambda x: int(x)),
    convert_days_to_months=True,
    add_hazard_ratio=False,
    add_logrank_pvalue=False,   # keep False; we’ll add global p after
    at_risk_counts=False,
    ax=ax,
    title="Overall survival by WES CNV Leiden cluster (autosomes only)"
)

df_lr = df_cluster[[ "OS", cluster_col ]].dropna().copy()
df_lr["time"] = df_lr["OS"].abs()
df_lr["event"] = df_lr["OS"] > 0

res = multivariate_logrank_test(
    event_durations=df_lr["time"],
    groups=df_lr[cluster_col].astype(str),
    event_observed=df_lr["event"]
)

p = float(res.p_value)
p_txt = "Global log-rank: p<0.001" if p < 0.001 else f"Global log-rank: p={p:.3g}"

ax.text(0.02, 0.02, p_txt, transform=ax.transAxes, fontsize=10, va="bottom")
plt.tight_layout()
savefig(
    df = df_km,
    out_file=FIGURES_DIR/"figure3/fig3c_km_mulvariate_test.pdf",
)
plt.show()

# %%
