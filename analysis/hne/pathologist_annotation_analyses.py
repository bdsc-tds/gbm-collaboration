"""
ECM neighborhood vs pathologist annotation association analysis.

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
      - 3-panel replication (All / UKER / CHUV) with shared scale
      - Dotplot (color=residual, size=significance)
      - Sparse/top-effects heatmap (only |resid|>=2 & q<0.05)
  6) Optional metrics:
      - Shannon diversity of annotations within each ECM neighborhood (and vice versa)
      - Jaccard similarity of "enriched associations" between centers

"""
# %%
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import chi2_contingency, chi2, norm
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import statsmodels.formula.api as smf

# External util in your environment (from your original script)
from gbmcollaboration.utils import concat_pathology_annotations_from_s3


# ----------------------------
# Config
# ----------------------------

@dataclass(frozen=True)
class DataConfig:
    ecm_cluster_path: str
    s3_bucket: str
    s3_prefixes: Tuple[str, ...]
    sample_col: str = "sample_name"
    barcode_col: str = "Barcode"
    annot_col: str = "pathology_annotations"
    ecm_col: str = "ECM_neighborhood"
    spot_id_col: str = "spot_id"


@dataclass(frozen=True)
class AnalysisConfig:
    # Filtering
    min_spots_per_annotation: int = 0  # 0 keeps all
    keep_only_shared_annotations: bool = True
    shared_annotations: Tuple[str, ...] = (
        "blood vessels",
        "cortex with tumor infiltration",
        "cortex without visible tumor infiltration",
        "infiltration zone",
        "necrosis",
        "tumor",
    )

    # plotting
    clip_residuals: float = 3.0
    q_threshold: float = 0.05


# ----------------------------
# Load + preprocessing
# ----------------------------

def load_ecm_cluster(path: str, index_col: str = "spot_id") -> pd.DataFrame:
    ecm = pd.read_csv(path)
    if index_col in ecm.columns:
        ecm = ecm.set_index(index_col)
    else:
        raise ValueError(f"ECM file must contain column '{index_col}'.")
    return ecm


def load_annotations_from_s3(bucket: str, prefixes: Iterable[str]) -> pd.DataFrame:
    dfs = []
    for p in prefixes:
        dfs.append(concat_pathology_annotations_from_s3(bucket=bucket, prefix=p))
    out = pd.concat(dfs, ignore_index=False)
    return out


def standardize_annotations(df: pd.DataFrame, annot_col: str) -> pd.DataFrame:
    """Lowercase + fix known typos/variants from your script."""
    out = df.copy()

    out[annot_col] = out[annot_col].astype(str).str.lower()

    mapping = {
        "necrose": "necrosis",
        "white matter with tumor": "white matter with tumor infiltration",
        "white matter with tumor ": "white matter with tumor infiltration",
        "blood vessel": "blood vessels",
        "cortex without visible tumor cell infiltration": "cortex without visible tumor infiltration",
        "cortex": "cortex without visible tumor infiltration",
        "gray matter (cortex) with tumor infiltration": "cortex with tumor infiltration",
        "non-tumoral (gray matter)": "cortex without visible tumor infiltration",
        "tumour (microvascular proliferation)": "tumor (microvascular proliferation)",
        "tumor (microvascular prolifeation)": "tumor (microvascular proliferation)",
        "tumor (mircovascular proliferation)": "tumor (microvascular proliferation)",
        "tumor (microvascular proliferations)": "tumor (microvascular proliferation)",
        "tumor (microvacular proliferation)": "tumor (microvascular proliferation)",
        "tumor (with microvascular proliferation)": "tumor (microvascular proliferation)",
    }
    out[annot_col] = out[annot_col].replace(mapping)
    return out


def exclude_bad_spots(df: pd.DataFrame, annot_col: str) -> pd.DataFrame:
    """Exclude folds/shifted/no tissue items from your script."""
    bad = {
        "tumor (section fold)",
        "tumor (attention: spots seem shifted compared to the h&e)",
        "no tissue",
        "empty (detachment, maybe tumor in vis)",
        "potentially tumor in vis (fold in he) ",
        "tumor (attention: spots seem shifted compared to h&e ) on h&e the smaller part is superimposed on the bigger one",
        "no tissue (attention: spots seem shifted compared to h&e)",
        "potentially tumor in vis (fold in he)",
    }
    out = df[~df[annot_col].isin(bad)].copy()
    out = out.dropna(subset=[annot_col])
    return out


def build_spot_id_and_index(df: pd.DataFrame, sample_col: str, barcode_col: str, spot_id_col: str) -> pd.DataFrame:
    out = df.copy()

    # Your one-off correction
    out[sample_col] = out[sample_col].replace({"CH_G_1789c_vis": "CH_G_189c_vis"})

    out[spot_id_col] = out[sample_col].astype(str) + "_" + out[barcode_col].astype(str)
    out = out.set_index(spot_id_col)
    return out


def add_center_from_sample(df: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    mapping = {"E": "UKER", "C": "CHUV"}
    out = df.copy()
    first = out[sample_col].astype(str).str[0].fillna("")
    out["center"] = first.map(mapping).fillna("Other")
    return out


def merge_annotation_ecm(df_annot: pd.DataFrame, ecm: pd.DataFrame) -> pd.DataFrame:
    return pd.merge(df_annot, ecm, how="inner", left_index=True, right_index=True)


def filter_to_shared_annotations(df: pd.DataFrame, annot_col: str, shared: Iterable[str]) -> pd.DataFrame:
    shared = set(shared)
    return df[df[annot_col].isin(shared)].copy()


# ----------------------------
# Core stats
# ----------------------------

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


# ----------------------------
# Optional metrics: Shannon + Jaccard
# ----------------------------

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


# ----------------------------
# Figures
# ----------------------------

def _mpl_pub_defaults():
    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })


def fig_residual_heatmap(
    resid: pd.DataFrame,
    qvals: pd.DataFrame,
    title: str,
    clip: float,
    q_thr: float,
    out_pdf: str,
    out_png: str,
    order_by: Optional[str] = "tumor",
    figsize: Tuple[float, float] = (6.5, 5.5), #(10.8, 3.8),
):
    _mpl_pub_defaults()

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

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.show()


def fig_three_panel_centers(
    df: pd.DataFrame,
    ecm_col: str,
    annot_col: str,
    clip: float,
    q_thr: float,
    out_pdf: str,
    out_png: str,
    min_spots: int = 0,
    figsize: Tuple[float, float] = (10, 4.5), #(11.0, 3.6),
):
    _mpl_pub_defaults()

    # ALL reference ordering
    ct_all = build_contingency(df, ecm_col, annot_col, min_spots=min_spots)
    chi = chi2_independence(ct_all)
    resid_all = standardized_residuals(ct_all, chi["expected"])
    q_all = residuals_fdr(resid_all)

    if "tumor" in resid_all.columns:
        resid_all = resid_all.loc[resid_all["tumor"].sort_values(ascending=False).index]
        q_all = q_all.loc[resid_all.index, resid_all.columns]

    cols = resid_all.columns
    rows = resid_all.index

    panels = [("Both centers", df), ("UKER", df[df["center"] == "UKER"]), ("CHUV", df[df["center"] == "CHUV"])]

    fig, axes = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
    ims = []

    for ax, (name, sub) in zip(axes, panels):
        ct = build_contingency(sub, ecm_col, annot_col, min_spots=0).reindex(index=rows, columns=cols, fill_value=0)
        chi = chi2_independence(ct)
        resid = standardized_residuals(ct, chi["expected"]).clip(-clip, clip)
        qvals = residuals_fdr(standardized_residuals(ct, chi["expected"])).loc[rows, cols]
        sig = (qvals < q_thr).to_numpy()

        im = ax.imshow(resid.to_numpy(), aspect="auto", interpolation="nearest", vmin=-clip, vmax=clip)
        ims.append(im)

        ax.set_title(name, pad=6)

        ax.set_xticks(np.arange(len(cols)))
        ax.set_xticklabels(cols, rotation=45, ha="right")

        ax.set_yticks(np.arange(len(rows)))
        if name == "Both centers":
            ax.set_yticklabels(rows)
        else:
            ax.set_yticklabels([])

        ax.set_xlabel("")
        ys, xs = np.where(sig)
        ax.scatter(xs, ys, s=8)

        ax.set_xticks(np.arange(-.5, len(cols), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
        ax.grid(which="minor", linewidth=0.2)
        ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(ims[0], ax=axes, fraction=0.02, pad=0.02)
    cbar.set_label("Standardized residual")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.show()


def fig_dotplot(
    df: pd.DataFrame,
    ecm_col: str,
    annot_col: str,
    clip: float,
    out_pdf: str,
    out_png: str,
    min_spots: int = 0,
    figsize: Tuple[float, float] = (6.5, 5.5), #(11.0, 3.6),
):
    _mpl_pub_defaults()

    ct = build_contingency(df, ecm_col, annot_col, min_spots=min_spots)
    chi = chi2_independence(ct)
    resid = standardized_residuals(ct, chi["expected"])
    qvals = residuals_fdr(resid)

    # order by tumor if available
    if "tumor" in resid.columns:
        resid = resid.loc[resid["tumor"].sort_values(ascending=False).index]
        qvals = qvals.loc[resid.index, resid.columns]

    long = resid.clip(-clip, clip).stack().rename("resid").reset_index()
    long.columns = [ecm_col, annot_col, "resid"]
    long["q"] = qvals.stack().values
    long["mlog10q"] = (-np.log10(long["q"].clip(lower=1e-300))).clip(upper=10)

    xcats = list(resid.columns)
    ycats = list(resid.index)
    xmap = {k: i for i, k in enumerate(xcats)}
    ymap = {k: i for i, k in enumerate(ycats)}
    long["x"] = long[annot_col].map(xmap)
    long["y"] = long[ecm_col].map(ymap)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    sc = ax.scatter(
        long["x"], long["y"],
        c=long["resid"],
        s=8 + 10 * long["mlog10q"],
        vmin=-clip, vmax=clip
    )

    ax.set_xticks(range(len(xcats)))
    ax.set_xticklabels(xcats, rotation=45, ha="right")
    ax.set_yticks(range(len(ycats)))
    ax.set_yticklabels(ycats)
    ax.set_xlabel("")
    ax.set_ylabel("ECM neighborhood")
    ax.set_title("")

    cbar = fig.colorbar(sc, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Standardized residual")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.show()


def fig_sparse_top_effects(
    df: pd.DataFrame,
    ecm_col: str,
    annot_col: str,
    clip: float,
    q_thr: float,
    resid_thr: float,
    out_pdf: str,
    out_png: str,
    min_spots: int = 0,
    figsize: Tuple[float, float] = (6.5, 5.5), #(11.0, 3.6),
):
    _mpl_pub_defaults()

    ct = build_contingency(df, ecm_col, annot_col, min_spots=min_spots)
    chi = chi2_independence(ct)
    resid = standardized_residuals(ct, chi["expected"])
    qvals = residuals_fdr(resid)

    if "tumor" in resid.columns:
        resid = resid.loc[resid["tumor"].sort_values(ascending=False).index]
        qvals = qvals.loc[resid.index, resid.columns]

    keep = (qvals < q_thr) & (np.abs(resid) >= resid_thr)
    Z = resid.clip(-clip, clip).where(keep, other=np.nan)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.imshow(Z.to_numpy(), aspect="auto", interpolation="nearest", vmin=-clip, vmax=clip)
    ax.set_facecolor("0.92")

    ax.set_xticks(np.arange(Z.shape[1]))
    ax.set_xticklabels(Z.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(Z.shape[0]))
    ax.set_yticklabels(Z.index)
    ax.set_title(f"Top associations (|resid|≥{resid_thr}, q<{q_thr})")

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Standardized residual")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.show()


# ----------------------------
# Orchestration
# ----------------------------

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

    # main heatmap (publication style)
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
        out_png=os.path.join(outdir, "Fig_residual_heatmap.png"),
        order_by="tumor",
    )

    # replication panel
    fig_three_panel_centers(
        df=df,
        ecm_col=cfg.ecm_col,
        annot_col=cfg.annot_col,
        clip=acfg.clip_residuals,
        q_thr=acfg.q_threshold,
        min_spots=acfg.min_spots_per_annotation,
        out_pdf=os.path.join(outdir, "Fig_three_panel_centers.pdf"),
        out_png=os.path.join(outdir, "Fig_three_panel_centers.png"),
    )

    # dot plot
    fig_dotplot(
        df=df,
        ecm_col=cfg.ecm_col,
        annot_col=cfg.annot_col,
        clip=acfg.clip_residuals,
        min_spots=acfg.min_spots_per_annotation,
        out_pdf=os.path.join(outdir, "Fig_dotplot.pdf"),
        out_png=os.path.join(outdir, "Fig_dotplot.png"),
    )

    # sparse top effects
    fig_sparse_top_effects(
        df=df,
        ecm_col=cfg.ecm_col,
        annot_col=cfg.annot_col,
        clip=acfg.clip_residuals,
        q_thr=acfg.q_threshold,
        resid_thr=2.0,
        min_spots=acfg.min_spots_per_annotation,
        out_pdf=os.path.join(outdir, "Fig_sparse_top_effects.pdf"),
        out_png=os.path.join(outdir, "Fig_sparse_top_effects.png"),
    )


def load_and_prepare(cfg: DataConfig, acfg: AnalysisConfig) -> pd.DataFrame:
    ecm = load_ecm_cluster(cfg.ecm_cluster_path, index_col=cfg.spot_id_col)

    df_annot = load_annotations_from_s3(cfg.s3_bucket, cfg.s3_prefixes)
    df_annot = build_spot_id_and_index(df_annot, cfg.sample_col, cfg.barcode_col, cfg.spot_id_col)
    df_annot = standardize_annotations(df_annot, cfg.annot_col)
    df_annot = exclude_bad_spots(df_annot, cfg.annot_col)

    merged = merge_annotation_ecm(df_annot, ecm)
    merged = add_center_from_sample(merged, cfg.sample_col)

    if acfg.keep_only_shared_annotations:
        merged = filter_to_shared_annotations(merged, cfg.annot_col, acfg.shared_annotations)

    return merged


def main():
    # ---- edit these paths/prefixes for your environment ----
    cfg = DataConfig(
        ecm_cluster_path="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/Analysis_Repo/Notebooks/ecm_neighborhood_annotations_v4.csv",
        s3_bucket="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        s3_prefixes=(
            "fd319842-9157-4f04-bd4f-7a36ce6a0ba7/Non_mosaic_upload/Mosaic Annotations/",  # CHUV
            "2ca49430-8b5c-4504-ad96-344019c90e84/",  # UKER
        ),
    )

    acfg = AnalysisConfig(
        min_spots_per_annotation=0,
        keep_only_shared_annotations=True,
        clip_residuals=3.0,
        q_threshold=0.05,
    )

    df = load_and_prepare(cfg, acfg)
    results = run_analysis(df, cfg, acfg)
    make_all_figures(df, cfg, acfg, outdir="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/hne/ecm_pathology_figures")


if __name__ == "__main__":
    main()
# %%
