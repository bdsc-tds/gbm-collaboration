"""The paper cohort: 273 patients (one tumour block each), 260 in the survival analyses.

The clinical tables come from the regular loader, AllDataCenters.load_clinical. This module
only adds the cohort definition used throughout the paper: loader exclusions, block
ER_G_290a dropped (ER_G_290 has two baseline blocks), and `survival_eligible` = False for the
13 Charité patients of SURVIVAL_EXCLUSION (immature follow-up). `code_covariates` codes MGMT,
p53 and sex numerically, either keeping unknown values missing (complete-case fits) or
counting them as the reference level, as several published models did.
"""

import functools

import pandas as pd

from gbmcollaboration import SURVIVAL_EXCLUSION, AllDataCenters

from . import config as C
from .ids import block_of, patient_of

COVARIATE_LEVELS = {"mgmt": ("MGMT promoter methylation", "Methylated", "Unmethylated"),
                    "p53_pos": ("p53", "Positive", "Negative"),
                    "female": ("Administrative gender", "Female", "Male"),
                    "male": ("Administrative gender", "Male", "Female")}


@functools.lru_cache
def clinical() -> pd.DataFrame:
    """The paper cohort: one row per patient (index = block), 273 rows."""
    d = AllDataCenters.load_clinical(**C.LOADER_FLAGS).drop(index=C.DUPLICATE_BLOCK_DROP, errors="ignore").copy()
    d.index = block_of(d.index).rename("block")
    d.insert(0, "patient", patient_of(d.index).to_numpy())
    assert d["patient"].is_unique
    d["survival_eligible"] = ~d.index.isin(SURVIVAL_EXCLUSION)
    return d


def survival_cohort() -> pd.DataFrame:
    """The 260 patients of the survival analyses."""
    d = clinical()
    return d[d["survival_eligible"]]


def code_covariates(df: pd.DataFrame, missing: str = "keep") -> pd.DataFrame:
    """Add numeric covariates: mgmt (1 = methylated), p53_pos, female, male, age.

    missing="keep": unknown values stay NaN, so a model on them fits complete cases.
    missing="as_reference": (x == level), so unknown MGMT / p53 / sex count as
    unmethylated / negative / male, as in the source cells.
    """
    out = df.copy()
    for name, (col, pos, neg) in COVARIATE_LEVELS.items():
        if missing == "as_reference":
            out[name] = (df[col] == pos).astype(float)
        elif missing == "keep":
            out[name] = df[col].map({pos: 1.0, neg: 0.0})
        else:
            raise ValueError(missing)
    out["age"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce")
    return out


def n_coded_as_reference(df: pd.DataFrame, covariates=("mgmt", "p53_pos", "female")) -> pd.Series:
    """Rows of df whose unknown value the as-ran coding turned into the reference level."""
    return pd.Series({c: int(df[COVARIATE_LEVELS[c][0]].isna().sum()) for c in covariates})
