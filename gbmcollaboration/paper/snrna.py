"""snRNA-seq helpers shared by fig2a_suppfig2_figure2.py and fig2bc_suppfig2_marker_dotplots.ipynb.

The loader object (`AllDataCenters.load_singlecell`) holds log1p(counts) in .X and the counts
in layers["unnormalized_counts"]; obsm X_pca / X_umap are the frozen embeddings of the
annotation pipeline.
"""

import pandas as pd

REFERENCE_COLUMN = "celltype_level4_scanvi"
LEVELS = ["level1", "level2", "level3", "level4"]


def add_curated_levels(adata, key: pd.DataFrame, levels=LEVELS) -> None:
    """celltype_<level>_curated columns from the level-4 scANVI labels (figure2.py and
    marker_dotplots.ipynb). key: indexed by celltype_level4_scanvi, as annotation.celltype_key()
    or the marker workbook's key sheet; labels missing from the key are kept as they are."""
    for lvl in levels:
        print(f"Mapping {lvl}")
        ct_mapping = key[f"celltype_{lvl}_curated"]
        adata.obs[f"celltype_{lvl}_curated"] = adata.obs[REFERENCE_COLUMN].astype(str).replace(ct_mapping)
