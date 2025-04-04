"""Hyperparameters and definitions for the GBM MOSAIC data on abstra."""

import os
from pathlib import Path
from typing import Any
from gbmcollaboration.utils import get_s3_path
import boto3
from urllib.parse import urlparse
import subprocess


def get_zarr_directories(bucket_name, prefix=""):
    try:
        command = f"aws s3 ls s3://{bucket_name}/{prefix}"
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)

        # Filter .zarr folders
        zarr_folders = [line.split()[-1].rstrip('/') for line in result.stdout.splitlines() if line.endswith('.zarr/')]

        # Build full paths
        full_zarr_paths = [f"s3://{bucket_name}/{prefix}{zarr_folder}" for zarr_folder in zarr_folders]
        return full_zarr_paths
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        print(e.output)
        return []


def generate_data_path(
    clinical_id: Path,
    bulkrna_id: Path,
    wes_id: Path,
    he_images_id: Path,
    he_features_id: Path,
    load_zarr_func: list[Path],
    single_cell_id: Path,
):
    """Generate a dictionary of data paths."""
    data_paths = {
        "Processed clinical": get_s3_path(clinical_id, "tables/processed_clinical_data.csv"),
        "Key events": get_s3_path(clinical_id, "tables/patient_key_events.csv"),
        "Bulk RNA raw counts": get_s3_path(bulkrna_id, "raw_counts.tsv"),
        "Bulk RNA TPM": get_s3_path(bulkrna_id, "tpm_counts.tsv"),
        "Bulk RNA normalized counts": get_s3_path(bulkrna_id, "deseq2norm_counts.tsv"),
        "WES CNV deletion": get_s3_path(wes_id, "cohort_level/cnv/binary_deletion_status.csv"),
        "WES CNV amplification": get_s3_path(wes_id, "cohort_level/cnv/binary_duplication_status.csv"),
        "WES mutations": get_s3_path(wes_id, "cohort_level/snv_indel/binary_oncogenic_alteration_status.csv"),
        "HE_images": get_s3_path(he_images_id, "HE"),
        "HE_features_H0": get_s3_path(he_features_id, "image_features/features_H0"),
        "Spatial Zarr Files": load_zarr_func,
        "Single-cell RNA": str(get_s3_path(single_cell_id, f"v3/dataset_level/dataset-{single_cell_id}/merge_samples/sc_merged_annotated.zarr")),
    }
    return data_paths

### CHUV DATA
DATA_ID_CHUV = {
    "CHUV_clinical": "d9bbbfeb-da3b-47cb-b31f-badbf030f28d",
    "CHUV_bulkrna": "2dc2baed-7851-49bd-a53b-7a9f18ec2a5c",
    "CHUV_wes": "cf474777-4c58-45d3-8392-cc5482bbc334",
    "CHUV_he_images": "a3857b67-3f28-40e6-b1f9-4842a2ad49af",
    "CHUV_he_features": "b5bd8840-cc6e-4de2-a78a-22315e9d1734",
    "CHUV_sc_rna": "34d9d6ca-e701-4f7f-9ea5-3f238a267076"
}

### UKER DATA
DATA_ID_UKER = {
    "UKER_clinical": "01af163d-4da2-41ce-b7f3-8135bf640d65",
    "UKER_bulkrna": "558a04f3-1334-4596-b692-31525dd88d89",
    "UKER_wes": "d239b19d-fe6e-4362-8a1c-35b4182c2511",
    "UKER_he_images": "3205eea0-5076-4294-8876-08cf98474603",
    "UKER_he_features": "40ce7e61-5ff6-467e-9972-0c85705f086c",
    "UKER_sc_rna": "09ce5d6d-e8de-464a-91c1-c534823200d0"
}

# Check if running in GitHub Actions
IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
# Conditionally set up the data path dictionary
if not IN_GITHUB_ACTIONS:
    load_zarr_chuv_files = get_zarr_directories(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix="b5bd8840-cc6e-4de2-a78a-22315e9d1734/v3/without_spotclean/SpatialData/"
    )
    load_zarr_uker_files = get_zarr_directories(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix="40ce7e61-5ff6-467e-9972-0c85705f086c/v3/without_spotclean/SpatialData/"
    )
else:
    load_zarr_chuv_files = []
    load_zarr_uker_files = []

DATA_PATH_CHUV_FILES = generate_data_path(
    DATA_ID_CHUV["CHUV_clinical"],
    DATA_ID_CHUV["CHUV_bulkrna"],
    DATA_ID_CHUV["CHUV_wes"],
    DATA_ID_CHUV["CHUV_he_images"],
    DATA_ID_CHUV["CHUV_he_features"],
    load_zarr_chuv_files,
    DATA_ID_CHUV["CHUV_sc_rna"],
)

DATA_PATH_UKER_FILES = generate_data_path(
    DATA_ID_UKER["UKER_clinical"],
    DATA_ID_UKER["UKER_bulkrna"],
    DATA_ID_UKER["UKER_wes"],
    DATA_ID_UKER["UKER_he_images"],
    DATA_ID_UKER["UKER_he_features"],
    load_zarr_uker_files,
    DATA_ID_UKER["UKER_sc_rna"],
)

MOSAIC_DATA_DICT_FILENAME = Path(__file__).parent / "metadata" / "gbm_dict.xlsx"

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
