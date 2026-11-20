"""
ECM neighborhood vs pathologist annotation association analysis.

Produces Fig 4f (Fig_residual_heatmap: standardized residuals of ECM neighborhood x pathology
annotation, dots at BH q < 0.05) and the chi-square / batch-adjusted test statistics
(figures/paper/fig4f_pathologist_annotation_analyses/, Fig 4f also in the panels folder).

Origin: analysis/hne/pathologist_annotation_analyses.py (Owkin). Comparison: the Fig 4f table is
compared at the end with the published one (results/hne/ecm_pathology_figures/
Fig_residual_heatmap.csv): same rows, columns and order, identical values.

Pipeline:
  1) Load ECM table + annotations from S3
  2) Clean / standardize annotations, build spot_id
  3) Merge + add center from sample_name (E->UKER, C->CHUV)
  4) Statistics:
      - Chi-square test of independence + Cramér's V
      - Standardized residuals + per-cell FDR (BH)
      - Batch-adjusted interaction test controlling for center (Negative Binomial GLM)
  5) Figures:
      - Publication-style residual heatmap with significance dots
  6) Optional metrics:
      - Shannon diversity of annotations within each ECM neighborhood (and vice versa)
      - Jaccard similarity of "enriched associations" between centers

"""
# %% -- Load packages
import os
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import chi2_contingency, chi2, norm
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import statsmodels.formula.api as smf

from gbmcollaboration.paper import config, output
from gbmcollaboration.utils import pathology_annotation_loader as patho
from gbmcollaboration.utils.pathology_annotation_loader import AnalysisConfig, DataConfig

OUT = output.setup("fig4f_pathologist_annotation_analyses")

plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
})

# -- Functions

def build_contingency(df: pd.DataFrame, ecm_col: str, annot_col: str, min_spots: int = 0) -> pd.DataFrame:
    ct = pd.crosstab(df[ecm_col], df[annot_col])
    if min_spots and min_spots > 0:
        keep = ct.sum(axis=0)
        keep = keep[keep > min_spots].index
        ct = ct[keep]
    return ct


def chi2_independence(ct: pd.DataFrame) -> Dict[str, object]:
    chi2_stat, p, dof, expected = chi2_contingency(ct)
    expected_df = pd.DataFrame(expected, index=ct.index, columns=ct.columns)
    return {"chi2": float(chi2_stat), "p": float(p), "dof": int(dof), "expected": expected_df}


def cramers_v(ct: pd.DataFrame, chi2_stat: Optional[float] = None) -> float:
    n = ct.to_numpy().sum()
    r, k = ct.shape
    if chi2_stat is None:
        chi2_stat = chi2_contingency(ct)[0]
    phi2 = chi2_stat / n
    return float(np.sqrt(phi2 / min(k - 1, r - 1)))


def standardized_residuals(ct: pd.DataFrame, expected_df: pd.DataFrame) -> pd.DataFrame:
    expected_df = expected_df.clip(lower=1e-12)
    return (ct - expected_df) / np.sqrt(expected_df)


def residuals_fdr(resid: pd.DataFrame) -> pd.DataFrame:
    # approximate N(0,1) under H0
    pvals = 2 * norm.sf(np.abs(resid.to_numpy()))
    pvals = pd.DataFrame(pvals, index=resid.index, columns=resid.columns)
    qflat = multipletests(pvals.to_numpy().ravel(), method="fdr_bh")[1]
    qvals = pd.DataFrame(qflat.reshape(pvals.shape), index=resid.index, columns=resid.columns)
    return qvals


def batch_adjusted_interaction_test(
    df: pd.DataFrame,
    ecm_col: str,
    annot_col: str,
    center_col: str = "center",
    min_spots: int = 0,
    family: str = "negbin",
) -> Dict[str, object]:
    """
    Log-linear model on aggregated counts, testing ECM:annotation interaction while controlling for center.
    Compare:
      H0: ECM + annot + center + ECM:center + annot:center
      H1: H0 + ECM:annot
    Uses Negative Binomial by default (more robust than Poisson for overdispersion).
    """

    tmp = df[[ecm_col, annot_col, center_col]].dropna().copy()

    if min_spots and min_spots > 0:
        keep = tmp[annot_col].value_counts()
        keep = keep[keep > min_spots].index
        tmp = tmp[tmp[annot_col].isin(keep)]

    counts = tmp.groupby([ecm_col, annot_col, center_col]).size().reset_index(name="count")

    # expand grid to include zeros
    ecm_levels = counts[ecm_col].unique()
    annot_levels = counts[annot_col].unique()
    center_levels = counts[center_col].unique()
    full = pd.MultiIndex.from_product([ecm_levels, annot_levels, center_levels],
                                      names=[ecm_col, annot_col, center_col]).to_frame(index=False)
    counts_full = full.merge(counts, on=[ecm_col, annot_col, center_col], how="left")
    counts_full["count"] = counts_full["count"].fillna(0).astype(int)

    f0 = (
        f"count ~ C({ecm_col}) + C({annot_col}) + C({center_col})"
        f" + C({ecm_col}):C({center_col}) + C({annot_col}):C({center_col})"
    )
    f1 = f0 + f" + C({ecm_col}):C({annot_col})"

    fam = sm.families.NegativeBinomial() if family.lower().startswith("neg") else sm.families.Poisson()

    m0 = smf.glm(formula=f0, data=counts_full, family=fam).fit()
    m1 = smf.glm(formula=f1, data=counts_full, family=fam).fit()

    lr_stat = 2 * (m1.llf - m0.llf)
    df_diff = int(m1.df_model - m0.df_model)
    p = chi2.sf(lr_stat, df_diff)

    return {"lr_stat": float(lr_stat), "df": df_diff, "p": float(p), "model_null": m0, "model_full": m1}


# -- Optional metrics: Shannon + Jaccard

def shannon_from_counts(counts: np.ndarray) -> float:
    """Shannon entropy (natural log)."""
    counts = counts.astype(float)
    total = counts.sum()
    if total <= 0:
        return np.nan
    p = counts / total
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def shannon_diversity_by_group(df: pd.DataFrame, group_col: str, label_col: str) -> pd.Series:
    """
    Shannon diversity of label distribution inside each group.
    Example: for each ECM neighborhood, diversity of pathology labels.
    """
    out = {}
    for g, sub in df.groupby(group_col):
        c = sub[label_col].value_counts().to_numpy()
        out[g] = shannon_from_counts(c)
    return pd.Series(out).sort_index()


def enriched_set(resid: pd.DataFrame, qvals: pd.DataFrame, q_thr: float = 0.05, resid_thr: float = 2.0) -> set:
    """Return set of (row, col) pairs that are significant & strong."""
    mask = (qvals < q_thr) & (np.abs(resid) >= resid_thr)
    pairs = list(zip(*np.where(mask.to_numpy())))
    return {(resid.index[i], resid.columns[j]) for i, j in pairs}


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return np.nan
    return len(a & b) / len(a | b)


def fig_residual_heatmap(
    resid: pd.DataFrame,
    qvals: pd.DataFrame,
    title: str,
    clip: float,
    q_thr: float,
    out_pdf: str,
    order_by: Optional[str] = "tumor",
    figsize: Tuple[float, float] = (6.5, 5.5),
    panel: Optional[str] = None,
):

    Z = resid.copy()
    if order_by is not None and order_by in Z.columns:
        Z = Z.loc[Z[order_by].sort_values(ascending=False).index]
        qvals = qvals.loc[Z.index, Z.columns]

    Zc = Z.clip(-clip, clip)
    sig = (qvals < q_thr).to_numpy()

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.imshow(Zc.to_numpy(), aspect="auto", interpolation="nearest", vmin=-clip, vmax=clip)

    ax.set_title(title, pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("ECM neighborhood")

    ax.set_yticks(np.arange(Zc.shape[0]))
    ax.set_yticklabels(Zc.index)
    ax.set_xticks(np.arange(Zc.shape[1]))
    ax.set_xticklabels(Zc.columns, rotation=45, ha="right")

    ys, xs = np.where(sig)
    ax.scatter(xs, ys, s=8)  # default color (black)

    # cell grid
    ax.set_xticks(np.arange(-.5, Zc.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, Zc.shape[0], 1), minor=True)
    ax.grid(which="minor", linewidth=0.2)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Standardized residual (enrichment ↔ depletion)")

    output.save(os.path.dirname(out_pdf), os.path.splitext(os.path.basename(out_pdf))[0], fig, Zc, panel=panel)
    plt.show()
    return Zc


def run_analysis(df: pd.DataFrame, cfg: DataConfig, acfg: AnalysisConfig) -> Dict[str, object]:
    """
    Run main tests and return key objects.
    """
    ct = build_contingency(df, cfg.ecm_col, cfg.annot_col, min_spots=acfg.min_spots_per_annotation)
    chi = chi2_independence(ct)
    V = cramers_v(ct, chi2_stat=chi["chi2"])
    resid = standardized_residuals(ct, chi["expected"])
    qvals = residuals_fdr(resid)

    batch = batch_adjusted_interaction_test(
        df,
        ecm_col=cfg.ecm_col,
        annot_col=cfg.annot_col,
        center_col="center",
        min_spots=acfg.min_spots_per_annotation,
        family="negbin",
    )

    # Shannon diversity summaries
    sh_ecm = shannon_diversity_by_group(df, group_col=cfg.ecm_col, label_col=cfg.annot_col)
    sh_annot = shannon_diversity_by_group(df, group_col=cfg.annot_col, label_col=cfg.ecm_col)

    # Jaccard similarity of enriched pairs between centers (optional sanity check)
    jac = np.nan
    if set(df["center"].unique()) >= {"UKER", "CHUV"}:
        # compute enriched sets per center using same table columns/rows (for comparability)
        sub_u = df[df["center"] == "UKER"]
        sub_c = df[df["center"] == "CHUV"]

        # align both to overall table structure
        ct_u = build_contingency(sub_u, cfg.ecm_col, cfg.annot_col, min_spots=0).reindex(
            index=ct.index, columns=ct.columns, fill_value=0
        )
        ct_c = build_contingency(sub_c, cfg.ecm_col, cfg.annot_col, min_spots=0).reindex(
            index=ct.index, columns=ct.columns, fill_value=0
        )
        chi_u = chi2_independence(ct_u)
        chi_c = chi2_independence(ct_c)
        resid_u = standardized_residuals(ct_u, chi_u["expected"])
        resid_c = standardized_residuals(ct_c, chi_c["expected"])
        q_u = residuals_fdr(resid_u)
        q_c = residuals_fdr(resid_c)

        A = enriched_set(resid_u, q_u, q_thr=acfg.q_threshold, resid_thr=2.0)
        B = enriched_set(resid_c, q_c, q_thr=acfg.q_threshold, resid_thr=2.0)
        jac = jaccard(A, B)

    print("\n=== Overall (no batch adjustment) ===")
    print(f"Chi2 = {chi['chi2']:.2f}, dof = {chi['dof']}, p = {chi['p']:.3e}")
    print(f"Cramér's V = {V:.3f}")

    print("\n=== Batch-adjusted (controls for center; NB GLM) ===")
    print(f"LR stat = {batch['lr_stat']:.2f}, df = {batch['df']}, p = {batch['p']:.3e}")

    print("\n=== Optional metrics ===")
    print("Shannon diversity of pathology within each ECM neighborhood (higher = more mixed):")
    print(sh_ecm.sort_values(ascending=False).round(3))
    print("\nShannon diversity of ECM neighborhoods within each pathology label:")
    print(sh_annot.sort_values(ascending=False).round(3))
    if not np.isnan(jac):
        print(f"\nJaccard similarity of enriched associations (UKER vs CHUV) = {jac:.3f}")

    return {
        "ct": ct,
        "chi2": chi,
        "cramers_v": V,
        "residuals": resid,
        "qvals": qvals,
        "batch": batch,
        "shannon_by_ecm": sh_ecm,
        "shannon_by_annot": sh_annot,
        "jaccard_enriched_centers": jac,
    }


def make_all_figures(df: pd.DataFrame, cfg: DataConfig, acfg: AnalysisConfig, outdir: str = "."):
    os.makedirs(outdir, exist_ok=True)

    # main heatmap
    ct = build_contingency(df, cfg.ecm_col, cfg.annot_col, min_spots=acfg.min_spots_per_annotation)
    chi = chi2_independence(ct)
    resid = standardized_residuals(ct, chi["expected"])
    qvals = residuals_fdr(resid)

    fig_residual_heatmap(
        resid=resid,
        qvals=qvals,
        title="",
        clip=acfg.clip_residuals,
        q_thr=acfg.q_threshold,
        out_pdf=os.path.join(outdir, "Fig_residual_heatmap.pdf"),
        order_by="tumor",
        panel="Fig 4f",
    )


def load_and_prepare(cfg: DataConfig, acfg: AnalysisConfig) -> pd.DataFrame:
    ecm = pd.read_csv(cfg.ecm_cluster_path).set_index(cfg.spot_id_col)

    # annotations from S3: standardised, bad spots / nan / arachnoid removed, optional shared subset
    df_annot = patho.load_and_prepare(cfg, acfg)

    merged = pd.merge(df_annot, ecm, how="inner", left_index=True, right_index=True)
    merged = patho.add_center_from_sample(merged, cfg.sample_col)
    return merged

# %% -- Run code to derive the plots
cfg = DataConfig(
    ecm_cluster_path=str(config.ECM_ANNOTATIONS_CSV),
    s3_bucket="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
    s3_prefixes=(
        "fd319842-9157-4f04-bd4f-7a36ce6a0ba7/Non_mosaic_upload/Mosaic Annotations/",  # CHUV
        "2ca49430-8b5c-4504-ad96-344019c90e84/",  # UKER
    ),
)

acfg = AnalysisConfig(
    min_spots_per_annotation=0,
    keep_only_shared_annotations=False,
    clip_residuals=3.0,
    q_threshold=0.05,
)

df = load_and_prepare(cfg, acfg)
results = run_analysis(df, cfg, acfg)
make_all_figures(df, cfg, acfg, outdir=str(OUT))

# Compare Fig 4f with the published data (results/hne/ecm_pathology_figures/)
pub = pd.read_csv(config.EFS / "results/hne/ecm_pathology_figures/Fig_residual_heatmap.csv", index_col=0)
new = pd.read_csv(OUT / "Fig_residual_heatmap.csv", index_col=0)
print(f"Fig 4f: same rows/columns/order as published: {pub.index.equals(new.index) and pub.columns.equals(new.columns)}; "
      f"max |difference| = {np.abs(pub.to_numpy() - new.to_numpy()).max():.2e}")

# %%
