"""Visium helpers shared by the analysis notebooks.

Loading itself goes through the regular loaders: AllDataCenters.load_spatial_adata (with
`genes=` to keep only the genes an analysis needs) and AllDataCenters.load_cell2location.
`sections()` lists the paper's 228 sections; `add_c2l_to_obs` is the "Merge C2L data" cell
that several notebooks repeated.
"""


import pandas as pd



C2L_PREFIX = "meanscell_abundance_w_sf_"


def sections() -> list[str]:
    """The 228 Visium sections of the paper (those of the frozen TIL zone file)."""
    from .til import til_zones

    return sorted(til_zones()["sample"].unique())


def add_c2l_to_obs(adata, c2l: pd.DataFrame) -> None:
    """Join the abundances onto adata.obs by (sample, barcode), as the 'Merge C2L data' cells did."""
    keys = pd.MultiIndex.from_arrays([adata.obs["sample"].to_numpy(), adata.obs_names.str.rsplit("_", n=1).str[-1]],
                                     names=["sample", "barcode"])
    add = c2l.reindex(keys)
    add.index = adata.obs.index
    adata.obs = adata.obs.join(add[[c for c in add.columns if c not in adata.obs.columns]])
