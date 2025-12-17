"""Script to associate WES CNV profiles to OS."""

# %% -- Import packages and load data
import pandas as pd
import numpy as np
from sksurv.ensemble import RandomSurvivalForest#, RandomForestRegressor
from sksurv.util import Surv
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines.utils import concordance_index
import shap

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
#X_scaled = X.copy()

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

# %% -- Load the WES CNV profiles, previously generated
df_wes_leiden_cluster = pd.read_csv(
    "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/infercnv/leiden_cluster_wes_baseline_samples_exclX.csv"
)
df_wes_leiden_cluster.set_index("orig.ident", inplace=True)

# %% -- Merge clusters with OS and build survival columns
# OS is encoded as signed time: event if OS>0, time = abs(OS)
df_surv = df_concat[["OS"]].copy()
df_surv["OS_event"] = df_surv["OS"] > 0
df_surv["OS_time"] = df_surv["OS"].abs()

df_cluster_surv = df_wes_leiden_cluster.merge(df_surv[["OS_time", "OS_event"]], how="inner",
                                              left_index=True, right_index=True)

# --- infer cluster column name robustly
cluster_candidates = [c for c in df_cluster_surv.columns if "leiden" in c.lower() or "cluster" in c.lower()]
if len(cluster_candidates) == 0:
    # fallback: take the first non-survival column
    cluster_col = [c for c in df_cluster_surv.columns if c not in ["OS_time", "OS_event"]][0]
else:
    cluster_col = cluster_candidates[0]

df_cluster_surv[cluster_col] = df_cluster_surv[cluster_col].astype(str)

print(f"Using cluster column: {cluster_col}")
print("Clusters:", sorted(df_cluster_surv[cluster_col].unique()))

# %% -- Cox: each cluster vs all others (univariate)
from lifelines import CoxPHFitter
import numpy as np

cox_rows = []
clusters = sorted(df_cluster_surv[cluster_col].unique())

for cl in clusters:
    tmp = df_cluster_surv[["OS_time", "OS_event", cluster_col]].copy()
    tmp["is_cluster"] = (tmp[cluster_col] == cl).astype(int)

    # Skip if degenerate
    n_in = int(tmp["is_cluster"].sum())
    n_out = int((1 - tmp["is_cluster"]).sum())
    if n_in < 5 or n_out < 5:
        cox_rows.append({
            "cluster": cl, "n_in": n_in, "n_out": n_out,
            "HR": np.nan, "CI95_low": np.nan, "CI95_high": np.nan, "p": np.nan
        })
        continue

    cph = CoxPHFitter()
    cph.fit(tmp[["OS_time", "OS_event", "is_cluster"]],
            duration_col="OS_time", event_col="OS_event")

    s = cph.summary.loc["is_cluster"]
    cox_rows.append({
        "cluster": cl,
        "n_in": n_in,
        "n_out": n_out,
        "HR": float(s["exp(coef)"]),
        "CI95_low": float(s["exp(coef) lower 95%"]),
        "CI95_high": float(s["exp(coef) upper 95%"]),
        "p": float(s["p"])
    })

df_cox_cluster_vs_rest = pd.DataFrame(cox_rows).sort_values("p")
print("\nCox (cluster vs all others):")
print(df_cox_cluster_vs_rest)

# %% -- KM curves: all clusters on the same plot
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(9, 7))

kmf = KaplanMeierFitter()
for cl in clusters:
    mask = (df_cluster_surv[cluster_col] == cl)
    if mask.sum() == 0:
        continue
    kmf.fit(
        durations=df_cluster_surv.loc[mask, "OS_time"],
        event_observed=df_cluster_surv.loc[mask, "OS_event"],
        label=f"Cluster {cl} (n={mask.sum()})"
    )
    kmf.plot_survival_function(ax=ax, ci_show=False)

ax.set_title("Overall survival by Leiden cluster (KM)")
ax.set_xlabel("Time")
ax.set_ylabel("Survival probability")
ax.legend(title="Leiden clusters", loc="best")
plt.tight_layout()
plt.show()

# %% -- Forest plot
################################################################
################ FOREST PLOT: CLUSTER VS REST ###################
################################################################

import matplotlib.pyplot as plt
import numpy as np

# Clean dataframe
df_forest = df_cox_cluster_vs_rest.copy()
df_forest = df_forest.dropna(subset=["HR"])
df_forest = df_forest.sort_values("HR")

# Y positions
y_pos = np.arange(len(df_forest))

fig, ax = plt.subplots(figsize=(8, 0.6 * len(df_forest) + 2))

# Plot HRs with CI
ax.errorbar(
    x=df_forest["HR"],
    y=y_pos,
    xerr=[
        df_forest["HR"] - df_forest["CI95_low"],
        df_forest["CI95_high"] - df_forest["HR"]
    ],
    fmt="o",
    color="black",
    ecolor="black",
    capsize=3
)

# Reference line
ax.axvline(1, color="grey", linestyle="--")

# Y-axis labels with sample size
ax.set_yticks(y_pos)
ax.set_yticklabels(
    [
        f"Cluster {cl} (n={n})"
        for cl, n in zip(df_forest["cluster"], df_forest["n_in"])
    ]
)

ax.set_xlabel("Hazard ratio (cluster vs all others)")
ax.set_title("Overall survival – Leiden clusters")

# Annotate p-values
for i, p in enumerate(df_forest["p"]):
    ax.text(
        x=df_forest["CI95_high"].max() * 1.05,
        y=i,
        s=f"p = {p:.3g}",
        va="center",
        fontsize=9
    )

plt.tight_layout()
plt.show()

################################################################
################ CATEGORICAL COX: LEIDEN CLUSTERS ################
################################################################

# %% -- Run a categorical CoxPH
from lifelines import CoxPHFitter

df_cox_cat = df_cluster_surv[[cluster_col, "OS_time", "OS_event"]].copy()

# Make clusters categorical
df_cox_cat[cluster_col] = df_cox_cat[cluster_col].astype("category")

# Choose reference cluster (largest by default = more stable)
ref_cluster = df_cox_cat[cluster_col].value_counts().idxmax()
print(f"Reference cluster: {ref_cluster}")

# Reorder categories so reference is first
new_cats = [ref_cluster] + [c for c in df_cox_cat[cluster_col].cat.categories if c != ref_cluster]
df_cox_cat[cluster_col] = df_cox_cat[cluster_col].cat.reorder_categories(new_cats)

# One-hot encode (drop reference)
df_cox_dummies = pd.get_dummies(df_cox_cat, columns=[cluster_col], drop_first=True)

# %%
cph_cat = CoxPHFitter()
cph_cat.fit(
    df_cox_dummies,
    duration_col="OS_time",
    event_col="OS_event"
)

cph_cat.print_summary()

# %%
summary = cph_cat.summary.copy()
summary["cluster"] = summary.index.str.replace(f"{cluster_col}_", "")
summary = summary.reset_index(drop=True)

summary = summary[[
    "cluster",
    "exp(coef)",
    "exp(coef) lower 95%",
    "exp(coef) upper 95%",
    "p"
]].rename(columns={
    "exp(coef)": "HR",
    "exp(coef) lower 95%": "CI95_low",
    "exp(coef) upper 95%": "CI95_high"
})

summary["reference"] = ref_cluster
summary = summary.sort_values("p")

print(summary)

# %%
plt.figure(figsize=(7, 0.5 * len(summary) + 2))

plt.errorbar(
    x=summary["HR"],
    y=summary["cluster"],
    xerr=[
        summary["HR"] - summary["CI95_low"],
        summary["CI95_high"] - summary["HR"]
    ],
    fmt="o"
)

plt.axvline(1, color="black", linestyle="--")
plt.xlabel("Hazard ratio (vs reference cluster)")
plt.ylabel("Leiden cluster")
plt.title(f"Cox model – Leiden clusters (ref = {ref_cluster})")
plt.tight_layout()
plt.show()

# %%
# CHECK ASSOCIATION WITH GENDER
merge_leiden_sex = df_cluster_surv.merge(df_concat[["Administrative gender"]], left_index=True, right_index=True)
# %%
