"""MOSAIC identifiers.

A patient is <centre>_G_<number> (CH_G_151), a tumour block adds one letter (CH_G_151a)
and a Visium section adds a suffix: _vis for the main section, _vbu, _vb2, _bu2 or _bu3
for the 21 backup sections. Parsing the whole ID keeps backup sections attached to their
block, which stripping a fixed number of characters does not.
"""

import numpy as np
import pandas as pd

ID_PATTERN = r"^([A-Z]{2}_G_\d+)([a-z]?)(_[A-Za-z0-9]+)?$"


def parse_ids(values) -> pd.DataFrame:
    """patient, block and suffix of each ID (NaN where an ID does not parse)."""
    s = pd.Series(np.asarray(values)).astype(str)
    m = s.str.extract(ID_PATTERN)
    return pd.DataFrame({"patient": m[0], "block": m[0] + m[1].fillna(""), "suffix": m[2].fillna("")})


def block_of(values) -> pd.Index:
    b = parse_ids(values)["block"]
    if b.isna().any():
        raise ValueError(f"Unparseable IDs: {list(pd.Series(values)[b.isna().to_numpy()])[:5]}")
    return pd.Index(b)


def patient_of(values) -> pd.Index:
    return pd.Index(parse_ids(values)["patient"])
