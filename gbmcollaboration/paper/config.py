"""Paths and constants shared by every script in paper/.

Cohort exclusions are not defined here: they are the lists in gbmcollaboration/
definitions.py, applied by the DataCenter loaders and by gbmcollaboration.apply_exclusions
for frozen tables that do not go through a loader. To remove a patient from every
analysis, add them to PATIENT_EXCLUSION there.
"""

import os
import subprocess
import sys
from pathlib import Path

# rpy2 needs R_HOME; use the R of the active conda environment
if "R_HOME" not in os.environ:
    _r = Path(sys.executable).parent / "R"
    if _r.exists():
        os.environ["R_HOME"] = subprocess.run([str(_r), "RHOME"], capture_output=True, text=True).stdout.strip()

REPO_ROOT = Path(__file__).resolve().parents[2]
EFS = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd")

# Figures, tables and printed numbers go to OUT_ROOT/<script name>/
OUT_ROOT = Path(os.environ.get("PAPER_OUT", EFS / "figures" / "paper"))

# --------------------------------------------------------------------------- cohort
# Loader exclusions used by every analysis of the paper (lists in definitions.py)
LOADER_FLAGS = dict(exclude_idh1_mut=True, exclude_patient_ids=True, exclude_relapse=True)
# ER_G_290 has two baseline blocks; where one block per patient is needed, block a
# (low tumour content, also in BAD_QUALITY_UKER_VISIUM) is dropped, as in the
# cohort description (273 patients).
DUPLICATE_BLOCK_DROP = ["ER_G_290a"]

CENTER_LABELS = {"CH": "CHUV", "ER": "UKER", "CR": "Charité"}

# --------------------------------------------------------------------------- frozen inputs
# Outputs of the upstream steps (analysis/paper/fig*_0_*), frozen on EFS when the paper was made
TIL_ZONES_CSV = EFS / "results/visium/til_zones/til_zone_annotations.csv"
ECM_ANNOTATIONS_CSV = EFS / "Analysis_Repo/Notebooks/ecm_neighborhood_annotations_v4.csv"
CT_ANNOTATIONS_CSV = EFS / "Analysis_Repo/Notebooks/celltype_neighborhood_annotations_k10.csv"
CNV_CLUSTERS_CSV = EFS / "results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv"
CNV_MATRIX_CSV = EFS / "results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv"
GENOMIC_WINDOWS_CSV = EFS / "results/scRNAseq/infercnv/genomic_windows.csv"
MUTATION_MATRIX_CSV = EFS / "results/WES/mutations/patient_gene_matrix.csv"
NICHENET_DIR = EFS / "results/visium/til_zones/nichenet_analysis"
ZONE_DE_DIR = EFS / "results/visium/til_zones/deg_pseudobulk_edger"
NICHENET_FILTERING_DIR = EFS / "results/visium/til_zones/nichenet_filtering"
XENIUM_BASE = EFS / "data"

# Annotation key: level-4 scANVI labels -> curated / figure labels
ANNOTATION_KEY_CSV = REPO_ROOT / "gbmcollaboration/metadata/SingleCell_Annotation_Key_v3b.csv"
GLIMMER_MARKERS_XLSX = REPO_ROOT / "gbmcollaboration/metadata/Cell type annotation markers Glimmer.xlsx"

# --------------------------------------------------------------------------- definitions
TIL_THRESHOLD = 2.0          # summed cell2location T/NK abundance for a TIL+ spot (>=)
TIL_HIGH_QUANTILE = 0.75     # TIL-high = top quartile of TIL+ spot fraction over all sections
TIL_ZONES = ["TIL+", "BZ1", "BZ2", "TIL-"]
RMST_TAU_YEARS = 3.0
