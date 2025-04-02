"""Hyperparameters and definitions for the GBM MOSAIC data on abstra."""

import os
from gbmcollaboration.utils.env_utils import Path


def get_zarr_directories(main_folder):
    """Get all zarr files in a folder and its subfolders."""

    folder = Path(main_folder)
    zarr_dirs = [
        p for p in folder.iterdir() if p.name.endswith(".zarr")
    ]

    return zarr_dirs


def generate_data_paths(
    clinical_path: Path,
    bulkrna_path: Path,
    wes_path: Path,
    he_path: Path,
    he_features_path: Path,
    load_zarr_func: list[Path],
    single_cell_path: Path,
):
    """Generate a dictionary of data paths."""
    data_paths = {
        "Processed clinical": clinical_path / "tables/processed_clinical_data.csv",
        "Key events": clinical_path / "tables/patient_key_events.csv",
        "Bulk RNA raw counts": bulkrna_path / "raw_counts.tsv",
        "Bulk RNA TPM": bulkrna_path / "tpm_counts.tsv",
        "Bulk RNA normalized counts": bulkrna_path / "deseq2norm_counts.tsv",
        "WES CNV deletion": wes_path / "cohort_level/cnv/binary_deletion_status.csv",
        "WES CNV amplification": wes_path
        / "cohort_level/cnv/binary_duplication_status.csv",
        "WES mutations": wes_path
        / "cohort_level/snv_indel/binary_oncogenic_alteration_status.csv",
        "HE": he_path,
        "HE_features_H0": he_features_path / "features_H0",
        "Spatial Zarr Files": load_zarr_func,
        "Single-cell RNA": single_cell_path / "sc_merged_annotated.zarr",
    }
    return data_paths


DATA_PATH = Path("s3://pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias")

### CHUV DATA
DATA_PATH_CHUV = {
    "CHUV_clinical": DATA_PATH / "d9bbbfeb-da3b-47cb-b31f-badbf030f28d",
    "CHUV_bulkrna": DATA_PATH / "2dc2baed-7851-49bd-a53b-7a9f18ec2a5c",
    "CHUV_wes": DATA_PATH / "cf474777-4c58-45d3-8392-cc5482bbc334",
    "CHUV_he": DATA_PATH / "a3857b67-3f28-40e6-b1f9-4842a2ad49af" / "HE",
    "CHUV_he_features": DATA_PATH
    / "b5bd8840-cc6e-4de2-a78a-22315e9d1734"
    / "image_features",
    "CHUV_spatial": DATA_PATH / "b5bd8840-cc6e-4de2-a78a-22315e9d1734",
    "CHUV_sc_rna": DATA_PATH
    / "34d9d6ca-e701-4f7f-9ea5-3f238a267076"
    / "v3"
    / "dataset_level"
    / "34d9d6ca-e701-4f7f-9ea5-3f238a267076"
    / "merge_samples",
}

### UKER DATA
DATA_PATH_UKER = {
    "UKER_clinical": DATA_PATH / "01af163d-4da2-41ce-b7f3-8135bf640d65",
    "UKER_bulkrna": DATA_PATH / "558a04f3-1334-4596-b692-31525dd88d89",
    "UKER_wes": DATA_PATH / "d239b19d-fe6e-4362-8a1c-35b4182c2511",
    "UKER_he": DATA_PATH / "3205eea0-5076-4294-8876-08cf98474603" / "HE",
    "UKER_he_features": DATA_PATH
    / "40ce7e61-5ff6-467e-9972-0c85705f086c"
    / "image_features",
    "UKER_spatial": DATA_PATH / "40ce7e61-5ff6-467e-9972-0c85705f086c",
    "UKER_sc_rna": DATA_PATH
    / "09ce5d6d-e8de-464a-91c1-c534823200d0"
    / "v3"
    / "dataset_level"
    / "09ce5d6d-e8de-464a-91c1-c534823200d0"
    / "merge_samples",
}

# Check if running in GitHub Actions
IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
# Conditionally set up the data path dictionary
if not IN_GITHUB_ACTIONS:
    load_zarr_chuv_files = get_zarr_directories(
        DATA_PATH_CHUV["CHUV_spatial"] / "without_spotclean/SpatialData"
    )
    load_zarr_uker_files = get_zarr_directories(
        DATA_PATH_UKER["UKER_spatial"] / "without_spotclean/SpatialData"
    )
else:
    load_zarr_chuv_files = []
    load_zarr_uker_files = []

DATA_PATH_CHUV_FILES = generate_data_paths(
    DATA_PATH_CHUV["CHUV_clinical"],
    DATA_PATH_CHUV["CHUV_bulkrna"],
    DATA_PATH_CHUV["CHUV_wes"],
    DATA_PATH_CHUV["CHUV_he"],
    DATA_PATH_CHUV["CHUV_he_features"],
    load_zarr_chuv_files,
    DATA_PATH_CHUV["CHUV_sc_rna"],
)

DATA_PATH_UKER_FILES = generate_data_paths(
    DATA_PATH_UKER["UKER_clinical"],
    DATA_PATH_UKER["UKER_bulkrna"],
    DATA_PATH_UKER["UKER_wes"],
    DATA_PATH_UKER["UKER_he"],
    DATA_PATH_UKER["UKER_he_features"],
    load_zarr_uker_files,
    DATA_PATH_UKER["UKER_sc_rna"],
)

MOSAIC_DATA_DICT_FILENAME = Path(__file__).parent / "third_party" / "gbm_dict.xlsx"

# Clinical preprocessing
COL_DATES = [
    "Event date",
    "Date of sample harvest",
    "Date of sample de-archival",
    "treatment_start_date | Adjuvant:Temozolomide",
    "treatment_start_date | Concurrent:Temozolomide ",
    "treatment_start_date | Maintenance:Temozolomide",
    "treatment_start_date | Post-operative:nan",
    "treatment_start_date | Surgery:PrimaryTumour",
    "treatment_start_date | Systemic therapy:Bevacizumab",
    "treatment_start_date | Systemic therapy:Carboplatin",
    "treatment_start_date | Systemic therapy:Lomustin",
    "treatment_start_date | Systemic therapy:Pembrolizumab",
    "treatment_start_date | Systemic therapy:Temozolomide",
    "Chemotherapy start date | Adjuvant:Temozolomide",
    "Chemotherapy start date | Concurrent:Temozolomide ",
    "Chemotherapy start date | Maintenance:Temozolomide",
    "Chemotherapy start date | Systemic therapy:Carboplatin",
    "Chemotherapy start date | Systemic therapy:Lomustin",
    "Chemotherapy start date | Systemic therapy:Temozolomide",
    "Targeted therapy end date | Systemic therapy:Bevacizumab",
    "Immunotherapy end date | Systemic therapy:Pembrolizumab",
    "Surgery date | Surgery:PrimaryTumour",
    "Targeted therapy start date | Systemic therapy:Bevacizumab",
    "Immunotherapy start date | Systemic therapy:Pembrolizumab",
    "Radiotherapy start date | Post-operative:nan",
    "Chemotherapy end date | Adjuvant:Temozolomide",
    "Chemotherapy end date | Concurrent:Temozolomide ",
    "Chemotherapy end date | Maintenance:Temozolomide",
    "Chemotherapy end date | Systemic therapy:Carboplatin",
    "Chemotherapy end date | Systemic therapy:Lomustin",
    "Chemotherapy end date | Systemic therapy:Temozolomide",
    "Radiotherapy end date | Post-operative:nan",
    "Date of cancer diagnosis",
    "Date of last follow-up before inclusion",
    "Date of death",
    "year_of_birth_datetime",
    "Date of signature",
    "Date of inclusion in MOSAIC",
    "Date of baseline whole blood count.1",
    "Date of liver function test",
    "Date of diagnosis",
]
# Columns to keep for the essential clinical variables
COL_ESSENTIALS = [
    "Administrative gender",
    "Age (years)",
    "MGMT promoter methylation",
    "Vital status at inclusion",
    "Sample collection chronology",
    "Sample source",
]
KEYWORDS_COL_TO_REMOVE = [
    "Breast Cancer",
    "BC -",
    "Bladder",
    "BL -",
    "DL -",
    "Mesothelioma",
    "MS -",
    "NSCLC",
    "NS -",
    "Ovarian",
    "OV -",
]

# HYPERPARAMETERS
AVG_NUM_DAYS_PER_MONTH = 30.437
AVG_NUM_DAYS_PER_YEAR = 365.25

# patient with IDH1 mutation, to potentially exclude
# from tabular data (not present in visium data)
PATIENTS_IDH1_MUT = ["ER_G_272a_vis", "ER_G_273a_vis", "ER_G_292b_vis", "ER_G_292a_vis"]

# samples to exclude from the UKER dataset
BAD_QUALITY_UKER_VISIUM = [
    "ER_G_196b_vis",
    "ER_G_199a_vis",
    "ER_G_200a_vis",
    "ER_G_203a_vis",
    "ER_G_209a_vis",
    "ER_G_210a_vbu",
    "ER_G_265a_vis",
    "ER_G_266a_vis",
    "ER_G_267a_vis",
    "ER_G_267b_vis",
    "ER_G_273a_vis",
    "ER_G_274a_vis",
    "ER_G_276b_vis",
    "ER_G_278a_vis",
    "ER_G_278b_vis",
    "ER_G_281a_vis",
    "ER_G_282a_vis",
    "ER_G_284a_vis",
    "ER_G_286a_vis",
    "ER_G_288a_vis",
    "ER_G_288b_vis",
    "ER_G_290a_vis",
    "ER_G_290b_vis",
    "ER_G_292a_vis",
    "ER_G_308b_vis",
]

# samples to remove from the CHUV dataset
BAD_QUALITY_CHUV_VISIUM = ["CH_G_157a_vis", "CH_G_178a_vis"]
