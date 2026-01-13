"""Script to associate WES CNV profiles to OS."""

# %% -- Import packages and load data
import pandas as pd
import numpy as np
from sksurv.ensemble import RandomSurvivalForest
from sklearn.ensemble import RandomForestRegressor
from sksurv.util import Surv
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines.utils import concordance_index
import shap
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter
from gbmcollaboration import (
    PATIENT_EXCLUSION,
    RELAPSE_EXCLUSION,
    PATIENTS_IDH1_MUT,
    SURVIVAL_EXCLUSION
)
from gbmcollaboration.utils.clinical_functions import (
    univariate,
    plot_forest,
    plot_correlations
)
from gbmcollaboration.utils.viz_survival import kmplot

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
################ PART 1: survival per bin ################
################################################################

# %% -- Load the WES CNV profiles, previously generated
df_wes_profile = pd.read_csv("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv")
df_wes_profile.rename(columns={"Unnamed: 0":"Tumour block ID for MOSAIC"}, inplace=True)
df_wes_profile.set_index("Tumour block ID for MOSAIC", inplace=True)

event_df = df_concat[["OS"]]
df_wes_profile = df_wes_profile.merge(event_df, how="inner", left_index=True, right_index=True)

# %% -- Run survival analysis using Random Forest, to see which bins are prognostic
# 1. Prepare features and clinical outcome
X = df_wes_profile.drop(columns=["OS"])
df_clinical = df_wes_profile[["OS"]]
df_clinical["OS_event"] = df_clinical["OS"] > 0
df_clinical["OS_time"] = df_clinical["OS"].abs()
df_clinical = df_clinical.drop(columns=["OS"])

X_scaled = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns, index=X.index)

# 2. Prepare survival outcome
y = Surv.from_dataframe("OS_event", "OS_time", df_clinical)

# 3. Fit Random Survival Forest
rsf = RandomSurvivalForest(
    n_estimators=1000,
    min_samples_split=10,
    min_samples_leaf=5,
    max_features="sqrt",
    n_jobs=-1,
    random_state=42
)
rsf.fit(X_scaled, y)

# %% Now use SHAP for feature importance
# 1. Scale CNV features
X_scaled = pd.DataFrame(StandardScaler().fit_transform(X), columns=X.columns, index=X.index)

# 2. Prepare survival object
y = Surv.from_dataframe("OS_event", "OS_time", df_clinical)

# 3. Fit Random Survival Forest
rsf = RandomSurvivalForest(
    n_estimators=1000,
    min_samples_split=10,
    min_samples_leaf=5,
    max_features="sqrt",
    n_jobs=-1,
    random_state=42
)
rsf.fit(X_scaled, y)

# 4. Predict risk score (negative risk = lower survival time = worse prognosis)
rsf_risk_scores = rsf.predict(X_scaled)

# 5. Train a regression model to mimic RSF risk scores
surrogate = RandomForestRegressor(n_estimators=100, random_state=42)
surrogate.fit(X_scaled, rsf_risk_scores)

# 6. Apply SHAP
explainer = shap.TreeExplainer(surrogate)
shap_values = explainer.shap_values(X_scaled)

# 7. Summary plot
shap.summary_plot(shap_values, X_scaled, max_display=20)

# %% -- Extract the top 20 most important features
# and do a univariate model to get the HR and KM
# Calculate mean absolute SHAP values for each feature
shap_importance = pd.DataFrame({
    "feature": X_scaled.columns,
    "mean_abs_shap": np.abs(shap_values).mean(axis=0)
}).sort_values(by="mean_abs_shap", ascending=False)

top_20_features = shap_importance.head(20)["feature"].tolist()
print("Top 20 SHAP-ranked CNV bins:")
print(top_20_features)

X_top20 = X_scaled[top_20_features]

# Make sure we are using SHAP values from the surrogate model
# You’ve already done this:
# explainer = shap.TreeExplainer(surrogate)
# shap_values = explainer.shap_values(X_scaled)

# For each top feature, make a dependence plot
for feat in top_20_features:
    shap.dependence_plot(feat, shap_values, X_scaled, show=True,
    interaction_index=feat)

# %% -- Load the clinical feature to also compare with known biomarker
known_bmk = df_concat[['Age (years)', 'Administrative gender', "Site code", 'MGMT promoter methylation', 'OS']]
X_top20_bmk = X_top20.merge(known_bmk, how="inner", left_index=True, right_index=True)

# -- Do a bit of processing for the modeling
var_names_all = {
    'Age (years)': 'continuous',
    'Administrative gender': 'categorical',
    'MGMT promoter methylation': 'categorical',
    'Site code': 'categorical',
    '1326': 'continuous',
    '1296': 'continuous',
    '475': 'continuous',
    '1304': 'continuous',
    '1301': 'continuous',
    '659': 'continuous',
    '1319': 'continuous',
    '922': 'continuous',
    '1217': 'continuous',
    '925': 'continuous',
    '935': 'continuous',
    '916': 'continuous',
    '981': 'continuous',
    '1333': 'continuous',
    '1334': 'continuous',
    '1216': 'continuous',
    '1310': 'continuous',
    '465': 'continuous',
    '1315': 'continuous',
    '448': 'continuous'
}
var_list = list(var_names_all.keys())

MAPPING = {
    "MGMT promoter methylation": {"Unmethylated": 0, "Methylated": 1},
    "Administrative gender": {"Male": 0, "Female": 1},
    "Site code": {"CH": 0, "ER": 1, "CR": 2}
}

df_map = X_top20_bmk.replace(MAPPING).astype(float)
df_map = df_map.reset_index().drop(columns=["Tumour block ID for MOSAIC"])
df = df_map[var_list + ["OS"]]
df.head()

# %% -- Plot correlations
plot_correlations(df_map, var_list=var_list, values='correlations')
# We see that profiles in some bins are highly correlated. The number follows
# so those bins are close genomically, which explains the correlation (Linkage disequilibrium)

# %% -- Run the univariate analyses
df_uni = univariate(df_map, endpoint="OS", var_names_all=var_names_all)
df_uni = df_uni.reset_index()

plot_forest(df_all=df_uni,
            distances=[0.4, 0.1, 0.3],
            figsize=(10, 7),
            outcome='OS',
            significance='FDR')

# Compute the C-index fo each clinical variable
for x_col in var_list:
    tmp = df_map[['OS', x_col]]
    tmp = tmp[~tmp[x_col].isna()]
    c_index_0 = concordance_index(
        event_times=tmp['OS'].abs(),
        predicted_scores=tmp[x_col],
        event_observed=tmp['OS'] >= 0,
    )
    # Try the other direction
    c_index_1 = concordance_index(
        event_times=tmp['OS'].abs(),
        predicted_scores=(-1) * tmp[x_col],
        event_observed=tmp['OS'] >= 0,
    )
    c_index = max(c_index_0, c_index_1)
    print('%s: C-index=%.3f (n=%d)' % (x_col, c_index, len(tmp)))

# %% -- Plot the KM survival plots
for col in ["1326", "981"]: # or df.columns:
    if col == "OS":
        continue
    print(col)
    df_kmplot = df[df[col].notna()]
    if col == "Age (years)":
        median = round(df_kmplot["Age (years)"].median())
        df_kmplot["Age category"] = df_kmplot["Age (years)"].map(lambda x: f"age < {median} yo" if x < median else f"age > {median} yo")
        col = "Age category"
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        fig = kmplot(df_kmplot, x="OS", hue=col, add_hazard_ratio=True, add_logrank_pvalue=True, ax=ax)
    if col == "1326":
        lim_log2ratio = -1
        df_kmplot["1326 category"] = df_kmplot["1326"].map(lambda x: f"log2ratio < {lim_log2ratio}" if x < lim_log2ratio else f"log2ratio > {lim_log2ratio}")
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        fig = kmplot(df_kmplot, x="OS", hue="1326 category", add_hazard_ratio=True, add_logrank_pvalue=True, ax=ax)
    if col == "981":
        lim_log2ratio = 0.2
        df_kmplot["981 category"] = df_kmplot["981"].map(lambda x: f"log2ratio < {lim_log2ratio}" if x < lim_log2ratio else f"log2ratio > {lim_log2ratio}")
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        fig = kmplot(df_kmplot, x="OS", hue="981 category", add_hazard_ratio=True, add_logrank_pvalue=True, ax=ax)

################################################################
################ PART 2: survival per leiden cluster ################
################################################################

# %% --- Load Leiden clusters (autosomes-only) previously computed
df_wes_leiden_cluster = pd.read_csv(
    "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
)
df_wes_leiden_cluster.set_index("orig.ident", inplace=True)

# %% Prepare dataframe
cluster_col = "cnv_leiden"
covariate_cols = ["Age (years)", "Administrative gender", "MGMT promoter methylation", "Site code"]

# Pull OS + covariates from merged clinical df
df_surv_cov = df_concat[["OS"] + covariate_cols].copy()

# Map categorical covariates
MAPPING = {
    "MGMT promoter methylation": {"Unmethylated": 0, "Methylated": 1},
    "Administrative gender": {"Male": 0, "Female": 1},
    "Site code": {"CH": 0, "ER": 1, "CR": 2},
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
    from statsmodels.stats.multitest import fdrcorrection
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
from statsmodels.stats.multitest import fdrcorrection
_, df_sum_cat["fdr"] = fdrcorrection(df_sum_cat["p"].values)

print("\nAdjusted categorical Cox summary (top hits):")
print(df_sum_cat.sort_values("p").head(20))

# Forest plot for categorical model (clusters + covariates)
plot_forest(
    df_all=df_sum_cat,
    distances=[0.55, 0.15, 0.40],
    figsize=(10, 6),
    outcome="OS",
    significance="FDR"  # or "FDR"
)
plt.title(f"Adjusted categorical Cox (ref cluster = {ref_cluster})")
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
df_km = df_cluster[["OS", cluster_col]].copy()
df_km[cluster_col] = df_km[cluster_col].astype(str)

fig, ax = plt.subplots(1, 1, figsize=(9, 7))
kmplot(
    data=df_km,
    x="OS",
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
plt.show()

# %% KM curves with multivariate_logrank_test
fig, ax = plt.subplots(1, 1, figsize=(9, 7))
kmplot(
    data=df_km,
    x="OS",
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
plt.show()

# %%
