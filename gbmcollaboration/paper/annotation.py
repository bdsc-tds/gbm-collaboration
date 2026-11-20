"""Cell-type annotation key: level-4 scANVI labels to curated and figure labels.

Replaces the "Load custom cell annotations from file" cell of the notebooks. They read an
EFS copy of the key (Analysis_Repo/Branch_infercnv_visium_backup/analysis/
SingleCell_Annotation_Key.csv), identical to the repo's SingleCell_Annotation_Key_v3b.csv
except for 10 celltype_level3_figure labels; `celltype_key(efs_copy=True)` reads that copy
for panels drawn with celltype_level3_figure.
"""

import pandas as pd

from . import config as C

LEVEL4_KEY = "celltype_level4_scanvi"
NEW_COLS = ["celltype_level4_curated", "celltype_level3_curated", "celltype_level3_figure",
            "celltype_level2_curated", "celltype_level2_figure", "celltype_level1_curated"]
EFS_KEY_CSV = C.EFS / "Analysis_Repo/Branch_infercnv_visium_backup/analysis/SingleCell_Annotation_Key.csv"


def celltype_key(efs_copy: bool = False) -> pd.DataFrame:
    """Key indexed by celltype_level4_scanvi with the NEW_COLS columns (strings)."""
    key_raw = pd.read_csv(EFS_KEY_CSV if efs_copy else C.ANNOTATION_KEY_CSV)
    missing = [c for c in [LEVEL4_KEY] + NEW_COLS if c not in key_raw.columns]
    if missing:
        raise ValueError(f"Missing columns in key file: {missing}")
    key = key_raw.drop_duplicates(subset=[LEVEL4_KEY]).copy()
    key[LEVEL4_KEY] = key[LEVEL4_KEY].astype(str)
    return key.set_index(LEVEL4_KEY)[NEW_COLS].apply(lambda s: s.astype(str))


def add_celltype_key(adata, efs_copy: bool = False) -> pd.DataFrame:
    """Store the key in adata.uns['celltype_annotation_key'] (as the notebooks did) and return it."""
    key = celltype_key(efs_copy)
    adata.uns["celltype_annotation_key"] = key
    return key


def map_levels(level4: pd.Series, efs_copy: bool = False) -> pd.DataFrame:
    """Curated / figure labels for a series of level-4 scANVI labels."""
    return celltype_key(efs_copy).reindex(level4.astype(str).to_numpy()).set_index(level4.index)
