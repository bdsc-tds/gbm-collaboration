"""Shared configuration for the Watson et al. (MOSAIC GBM) data release.

Every export script imports from here so that cohort filters, paths and
labels are defined exactly once.
"""

from pathlib import Path

from gbmcollaboration.paper import config as P

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
EFS = P.EFS

OUT_DIR = EFS / "exports" / "mosaic_gbm_watson2027"
WORK_DIR = OUT_DIR / "_work"  # caches, not part of the release
# Full-cohort build (all 273 patients): what the paper analysed and what 99_sanity_checks
# compares against. Not released; 07_release_filter.py writes the released subset to OUT_DIR.
BUILD_DIR = WORK_DIR / "full"

# Inputs produced by the paper's analyses (frozen spot/sample annotations)
# (shared with the paper scripts: gbmcollaboration/paper/config.py)
TIL_ZONE_CSV = P.TIL_ZONES_CSV
ECM_NEIGHBORHOOD_CSV = P.ECM_ANNOTATIONS_CSV
CNV_MATRIX_CSV = P.CNV_MATRIX_CSV
CNV_BINS_CSV = P.GENOMIC_WINDOWS_CSV
CNV_CLUSTERS_CSV = P.CNV_CLUSTERS_CSV
MUTATION_MATRIX_CSV = P.MUTATION_MATRIX_CSV
NICHENET_AXES = [P.NICHENET_DIR / z / "axes_selected.csv" for z in ("TILplus", "BZ1", "BZ2")]
ECM_GENES_CSV = EFS / "supplementary_tables/4a_ECM_genes.csv"  # = Supp Table 3a
GLIMMER_MARKERS_XLSX = P.GLIMMER_MARKERS_XLSX
ANNOTATION_KEY_CSV = P.ANNOTATION_KEY_CSV
MANUSCRIPT_TXT = WORK_DIR / "manuscript.txt"  # pdftotext of the manuscript

XENIUM_BASE = P.XENIUM_BASE

# --------------------------------------------------------------------------
# Cohort definitions (as in the paper)
# --------------------------------------------------------------------------
# Loader flags and the ER_G_290a drop (273 patients) are the paper's
LOADER_FLAGS = P.LOADER_FLAGS
DUPLICATE_BLOCK_DROP = P.DUPLICATE_BLOCK_DROP

# Sample-ID prefix -> released centre code
CENTRE_CODES = {"CH": "CH-1", "CR": "CR-2", "ER": "ER-3"}

# Released patients: only deceased patients may be shared (any recorded death in the
# clinical table, including Charité patients excluded from survival for immature follow-up).
RELEASE_DECEASED_ONLY = True

# TIL definition (Methods: "TIL+ spot definition", "TIL-high versus TIL-low")
TIL_THRESHOLD = P.TIL_THRESHOLD
TIL_HIGH_QUANTILE = P.TIL_HIGH_QUANTILE

# Cell-type level used for pseudobulks (72 types in the paper = level 4 curated)
CELLTYPE_COL = "celltype_level4_curated"
CELLTYPE_SOURCE_COL = "celltype_level4_scanvi"

# Highly variable genes
N_HVG = 3000
HVG_FLAVOR = "seurat_v3"

# Xenium (what generated Fig. 6b/6g: figures/xenium_fig_6.zip)
XENIUM_SEGMENTATION = "10x_mm_5um"
XENIUM_PANELS = ["hBrain_v1.2", "hAtlas_v1.1"]
XENIUM_CORRECTION = "split_fully_purified"
XENIUM_SPLIT_LEVEL = "celltype_level2_curated"
XENIUM_TILE_UM = 55.0
XENIUM_MIN_CELLS_PER_TILE = 5
XENIUM_CELLTYPE_COL = "celltype_level4_curated_first_type"
XENIUM_TIL_COL = "celltype_level3_curated_first_type"
XENIUM_TIL_TYPES = ["T_cell", "NK"]
