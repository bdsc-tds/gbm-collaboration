"""Helpers shared by the export scripts."""

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import config as C


def out(*parts) -> Path:
    """Path in the full-cohort build."""
    p = C.BUILD_DIR.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def release_out(*parts) -> Path:
    """Path in the released folder."""
    p = C.OUT_DIR.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


# Identifiers, the cohort and TIL groups are shared with the paper analyses
from gbmcollaboration.paper import cohort as _cohort  # noqa: E402
from gbmcollaboration.paper import til as _til  # noqa: E402
from gbmcollaboration.paper.ids import block_of, patient_of  # noqa: E402,F401


def centre_code(ids) -> pd.Index:
    return pd.Index(ids).astype(str).str[:2].map(C.CENTRE_CODES)


def definitions_list(name: str) -> list:
    """A list constant of gbmcollaboration/definitions.py (e.g. SURVIVAL_EXCLUSION)."""
    from gbmcollaboration import definitions

    return list(getattr(definitions, name))


def annotation_key() -> pd.DataFrame:
    key = pd.read_csv(C.ANNOTATION_KEY_CSV).drop_duplicates(C.CELLTYPE_SOURCE_COL)
    return key.set_index(C.CELLTYPE_SOURCE_COL).astype(str)


def load_clinical() -> pd.DataFrame:
    """Processed clinical tables of the three centres with the paper's loader exclusions."""
    from gbmcollaboration import AllDataCenters

    return AllDataCenters.load_clinical(**C.LOADER_FLAGS)


def til_labels() -> pd.DataFrame:
    """Per-Visium-section TIL stratification (gbmcollaboration.paper.til), indexed by block.

    Columns: visium_section, n_spots, n_til_pos, til_fraction, til_group (TIL-high = top
    quartile of TIL+ spot fraction over the 228 sections; TIL-low = the others, including
    the 23 without TIL+ spots; D3), til_group_zones (NA for sections without TIL+ spots,
    left out of TIL-zone comparisons). The cut-off is in .attrs["cutoff"].
    """
    g = _til.til_groups()
    label = {"TIL_high": "TIL-high", "TIL_low": "TIL-low"}
    out = pd.DataFrame({
        "visium_section": g.index,
        "n_spots": g["Total_Spots"].to_numpy(),
        "n_til_pos": g["TIL_spots"].to_numpy(),
        "til_fraction": g["TIL_fraction"].to_numpy(),
        "til_group": g["TIL_group"].map(label).to_numpy(),
        "til_group_zones": g["TIL_group_zones"].map(label).to_numpy(),
    }, index=pd.Index(g["block"], name="sample_id"))
    assert out.index.is_unique
    out.attrs["cutoff"] = g.attrs["cutoff"]
    return out


def sum_by_group(X, groups: pd.Series):
    """Sum rows of X (cells/spots x genes) per group label. Returns (matrix, group_index)."""
    codes, uniques = pd.factorize(groups, sort=True)
    ind = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (codes, np.arange(len(codes)))),
        shape=(len(uniques), len(codes)),
    )
    X = X if sp.issparse(X) else sp.csr_matrix(X)
    return (ind @ X).tocsr(), pd.Index(uniques)


def nnz_by_group(X, groups: pd.Series):
    """Number of genes detected (count>0 in at least one cell) per group."""
    X = X if sp.issparse(X) else sp.csr_matrix(X)
    Xb = X.copy()
    Xb.data = (Xb.data > 0).astype(np.float32)
    S, idx = sum_by_group(Xb, groups)
    return pd.Series(np.asarray((S > 0).sum(axis=1)).ravel(), index=idx)
