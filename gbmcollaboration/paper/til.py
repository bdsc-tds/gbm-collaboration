"""TIL+ spots, TIL zones and TIL-high / TIL-low sections.

The spot calls and zones are frozen in TIL_ZONES_CSV (pipeline: paper/pipelines/
p07_til_zones.py). A section is TIL-high when its TIL+ spot fraction reaches the upper
quartile over all 228 sections (cut-off 0.0304, 57 sections) and TIL-low otherwise (171),
including the 23 sections without any TIL+ spot. Comparisons between TIL zones leave
those 23 out, because border zones are defined around TIL+ spots (57 vs 148) (D3).
"""

import functools

import numpy as np
import pandas as pd

from . import config as C
from gbmcollaboration import apply_exclusions

from .ids import parse_ids

_BOOL = {"true": True, "false": False, "1": True, "0": False, "1.0": True, "0.0": False}


@functools.lru_cache
def til_zones() -> pd.DataFrame:
    """Spot table (index spot_id): sample (Visium section), patient, block, TIL_positive, TIL_zone.

    The zone file was written from the Visium cohort loaded with the standard exclusions,
    so applying them again removes nothing today; it keeps the spot set in line with the
    loader if the lists in definitions.py change.
    """
    z = pd.read_csv(C.TIL_ZONES_CSV, usecols=["spot_id", "sample", "TIL_positive", "TIL_zone"],
                    dtype={"spot_id": str, "sample": str, "TIL_zone": str})
    z["TIL_positive"] = z["TIL_positive"].astype(str).str.strip().str.lower().map(_BOOL)
    if z["TIL_positive"].isna().any():
        raise ValueError("Spots without a TIL_positive call")
    sections = pd.DataFrame(index=pd.Index(z["sample"].unique()))
    kept = apply_exclusions(sections, exclude_low_quality_visium=True, **C.LOADER_FLAGS).index
    z = z[z["sample"].isin(kept)]
    ids = parse_ids(z["sample"])
    z["patient"], z["block"] = ids["patient"].to_numpy(), ids["block"].to_numpy()
    return z.set_index("spot_id")


@functools.lru_cache
def til_groups(quantile: float = C.TIL_HIGH_QUANTILE) -> pd.DataFrame:
    """One row per Visium section (index = section name).

    Columns: patient, block, center, TIL_spots, Total_Spots, TIL_fraction, TIL_group
    (TIL_high / TIL_low over all sections), TIL_group_zones (NaN for sections without
    TIL+ spots). The cut-off is in .attrs["cutoff"].
    """
    z = til_zones()
    sec = z.groupby("sample").agg(TIL_spots=("TIL_positive", "sum"), Total_Spots=("TIL_positive", "size"))
    sec["TIL_spots"] = sec["TIL_spots"].astype(int)
    sec["TIL_fraction"] = sec["TIL_spots"] / sec["Total_Spots"]
    cutoff = sec["TIL_fraction"].quantile(quantile)
    sec["TIL_group"] = np.where(sec["TIL_fraction"] >= cutoff, "TIL_high", "TIL_low")
    sec["TIL_group_zones"] = sec["TIL_group"].where(sec["TIL_spots"] > 0)
    ids = parse_ids(sec.index)
    sec.insert(0, "patient", ids["patient"].to_numpy())
    sec.insert(1, "block", ids["block"].to_numpy())
    sec.insert(2, "center", sec["patient"].str[:2])
    sec.attrs["cutoff"] = float(cutoff)
    return sec


def zones_with_groups() -> pd.DataFrame:
    """Spot table with the section's TIL_group_zones attached (for zone comparisons)."""
    z = til_zones().copy()
    z["TIL_group"] = z["sample"].map(til_groups()["TIL_group_zones"])
    return z
