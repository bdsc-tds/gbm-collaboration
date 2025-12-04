"""Hyperparameters and definitions for the GBM MOSAIC data on abstra."""

import os
from pathlib import Path
from typing import Any
from gbmcollaboration.utils import get_s3_path
import boto3
from urllib.parse import urlparse
import subprocess

# choice of version of the processing for visium
FOLDER_VERSION = {
    "sc_rna": "v5",
    "visium": "v7",
    "clinical": "v3",
    "wes": "v2",
    "rnaseq":"v2"
}
# can be v1/v2/v3/live, live is not stable and can change through time


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


def load_wes_sample_level_metrics(bucket_name, search_path="sample_level_metrics.csv", prefix=""):
    try:
        command = f"aws s3 ls s3://{bucket_name}/{prefix}/sample_level/snv_indel/ --recursive"
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        matching_files = []
        for line in result.stdout.splitlines():
            parts = line.strip().split()
            if len(parts) == 4:  # aws s3 ls --recursive gives lines like: date time size path
                filepath = parts[-1]
                if filepath.endswith(search_path):
                    full_path = f"s3://{bucket_name}/{filepath}"
                    matching_files.append(full_path)
        return matching_files

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
    load_tmb_func: list[Path]
):
    """Generate a dictionary of data paths."""
    data_paths = {
        "Processed clinical": get_s3_path(clinical_id, f"{FOLDER_VERSION['clinical']}/tables/processed_clinical_data.csv"),
        "Key events": get_s3_path(clinical_id, f"{FOLDER_VERSION['clinical']}/tables/patient_key_events.csv"),
        "Bulk RNA raw counts": get_s3_path(bulkrna_id, f"{FOLDER_VERSION['rnaseq']}/raw_counts.tsv"),
        "Bulk RNA TPM": get_s3_path(bulkrna_id, f"{FOLDER_VERSION['rnaseq']}/tpm_counts.tsv"),
        "Bulk RNA normalized counts": get_s3_path(bulkrna_id, f"{FOLDER_VERSION['rnaseq']}/deseq2norm_counts.tsv"),
        "WES CNV deletion": get_s3_path(wes_id, f"{FOLDER_VERSION['wes']}/cohort_level/cnv/binary_deletion_status.csv"),
        "WES CNV amplification": get_s3_path(wes_id, f"{FOLDER_VERSION['wes']}/cohort_level/cnv/binary_duplication_status.csv"),
        "WES mutations": get_s3_path(wes_id, f"{FOLDER_VERSION['wes']}/cohort_level/snv_indel/binary_oncogenic_alteration_status.csv"),
        "HE_images": get_s3_path(he_images_id, "HE"),
        "HE_features_H0": get_s3_path(he_features_id, "image_features/features_H0"),
        "Spatial Zarr Files": load_zarr_func,
        "Single-cell RNA": str(get_s3_path(single_cell_id, f"{FOLDER_VERSION['sc_rna']}/dataset_level/dataset-{single_cell_id}/merge_samples/sc_merged_annotated.zarr")),
        "WES TMB": load_tmb_func,
        "Cell2location": get_s3_path(he_features_id, f"{FOLDER_VERSION['visium']}/without_spotclean/cell2loc/AnnData_objects/"),
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

### CR DATA
DATA_ID_CR = {
    "CR_clinical": "9e8255bd-f5fd-4be2-9b14-71f45c90ef77",
    "CR_bulkrna": "34bed287-1bf3-4a70-bff2-cb154061586b",
    "CR_wes": "53cde683-5a68-4340-a817-e92d74084281",
    "CR_he_images": "b0d0f0c7-a904-4331-ad1c-f86ca0527120",
    "CR_he_features": "f0415978-2385-4067-b215-13e699a0f7a7",
    "CR_sc_rna":"ebf3efa8-c8e8-47e1-a56f-68835d75c653"
}




# Check if running in GitHub Actions
IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
# Conditionally set up the data path dictionary
if not IN_GITHUB_ACTIONS:
    load_zarr_chuv_files = get_zarr_directories(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix=f"b5bd8840-cc6e-4de2-a78a-22315e9d1734/{FOLDER_VERSION['visium']}/without_spotclean/SpatialData/"
    )
    load_zarr_uker_files = get_zarr_directories(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix=f"40ce7e61-5ff6-467e-9972-0c85705f086c/{FOLDER_VERSION['visium']}/without_spotclean/SpatialData/"
    )
    load_zarr_cr_files = get_zarr_directories(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix=f"f0415978-2385-4067-b215-13e699a0f7a7/{FOLDER_VERSION['visium']}/without_spotclean/SpatialData/"
    )
    load_tmb_chuv = load_wes_sample_level_metrics(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix=DATA_ID_CHUV["CHUV_wes"]+f"/{FOLDER_VERSION['wes']}"
    )
    load_tmb_uker = load_wes_sample_level_metrics(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix=DATA_ID_UKER["UKER_wes"]+f"/{FOLDER_VERSION['wes']}"
    )
    load_tmb_cr = load_wes_sample_level_metrics(
        bucket_name="pa-bi4jtrdsa151p4d3m-jepjxh9qrwgxphe8hd3r6exicahnoeuw1b-s3alias",
        prefix=DATA_ID_CR["CR_wes"]+f"/{FOLDER_VERSION['wes']}"
    )
else:
    load_zarr_chuv_files = []
    load_zarr_uker_files = []
    load_zarr_cr_files = []
    load_tmb_chuv = []
    load_tmb_uker = []
    load_tmb_cr = []

DATA_PATH_CHUV_FILES = generate_data_path(
    DATA_ID_CHUV["CHUV_clinical"],
    DATA_ID_CHUV["CHUV_bulkrna"],
    DATA_ID_CHUV["CHUV_wes"],
    DATA_ID_CHUV["CHUV_he_images"],
    DATA_ID_CHUV["CHUV_he_features"],
    load_zarr_chuv_files,
    DATA_ID_CHUV["CHUV_sc_rna"],
    load_tmb_chuv
)

DATA_PATH_UKER_FILES = generate_data_path(
    DATA_ID_UKER["UKER_clinical"],
    DATA_ID_UKER["UKER_bulkrna"],
    DATA_ID_UKER["UKER_wes"],
    DATA_ID_UKER["UKER_he_images"],
    DATA_ID_UKER["UKER_he_features"],
    load_zarr_uker_files,
    DATA_ID_UKER["UKER_sc_rna"],
    load_tmb_uker,
)

DATA_PATH_CR_FILES = generate_data_path(
    DATA_ID_CR["CR_clinical"],
    DATA_ID_CR["CR_bulkrna"],
    DATA_ID_CR["CR_wes"],
    DATA_ID_CR["CR_he_images"],
    DATA_ID_CR["CR_he_features"],
    load_zarr_cr_files,
    DATA_ID_CR["CR_sc_rna"],
    load_tmb_cr,
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
PATIENTS_IDH1_MUT = ["ER_G_272a", "ER_G_273a", "ER_G_292b", "ER_G_292a"]

# samples to exclude: pending review/no data/blank from excel sheet
BAD_QUALITY_CHUV_VISIUM = [
    "CH_G_157a",
    "CH_G_188a",
    "CH_G_532a",
    "CH_G_533a",
    "CH_G_534a",
    "CH_G_535a",
    "CH_G_536a",
]
BAD_QUALITY_UKER_VISIUM = [
    "ER_G_196b",
    "ER_G_199a",
    "ER_G_200a",
    "ER_G_203a",
    "ER_G_204a",
    "ER_G_265a",
    "ER_G_266a",
    "ER_G_267b",
    "ER_G_278a",
    "ER_G_290a",
    "ER_G_301a",
    "ER_G_301b",
    "ER_G_308b",
    "ER_G_665b",
    "ER_G_744a",
    "ER_G_745a",
    "ER_G_746a",
    "ER_G_748a",
    "ER_G_749a",
    "ER_G_750a",
    "ER_G_752a",
    "ER_G_753a",
    "ER_G_754a",
    "ER_G_755a",
    "ER_G_756a",
    "ER_G_757a",
    "ER_G_298b",
    "ER_G_645b",
    "ER_G_747a",
    "ER_G_751a",
    "ER_G_759a",
    "ER_G_760a",
    "ER_G_664a",
    "ER_G_675a",
    "ER_G_758a",

    # duplicates: exclude all slides except the one marked passing QC
    "ER_G_209a_vis",
    "ER_G_210a_vbu",
    "ER_G_267a_vis",
    "ER_G_274a_vis",
    "ER_G_281a_vis",
    "ER_G_286a_vis",
    "ER_G_288a_vis",
    "ER_G_635a_vis",
    "ER_G_636a_vis",

    # exclude this one with missing deconvolution columns
    "ER_G_288a_vb2",


    # from Spencer (checking clinical data)
    "ER_G_195b",
    "ER_G_267b",
    "ER_G_274b",
    "ER_G_276b",
    "ER_G_278a",
    "ER_G_282b",
    "ER_G_284a",
    "ER_G_288b",
    "ER_G_290a",
    "ER_G_292b",
    "ER_G_298b",
    "ER_G_301a",
    "ER_G_301b",
    "ER_G_303b",
    "ER_G_308b",
    "ER_G_310b",
    "ER_G_645b",
    "ER_G_665b",
    "ER_G_272a",
    "ER_G_664a",
    "ER_G_675a"
]

BAD_QUALITY_CR_VISIUM = []

# list of samples from the "bad" batch of probes
# can be used to exclude them in the analysis to ensure
# it has no strong effect on the results
PROBE_ISSUE_CHUV = [
    "CH_G_191a",
    "CH_G_199a",
    "CH_G_192a",
    "CH_G_186b",
    "CH_G_190a",
    "CH_G_207a",
    "CH_G_164b",
    "CH_G_198a",
    "CH_G_194a",
    "CH_G_180a",
    "CH_G_193a",
    "CH_G_189c",
    "CH_G_205a",
    "CH_G_206a",
    "CH_G_197a",
    "CH_G_203a",
    "CH_G_194b",
    "CH_G_202a",
    "CH_G_195a",
    "CH_G_179a",
    "CH_G_201a",
    "CH_G_189a",
    "CH_G_189b",
    "CH_G_168b",
    "CH_G_196a",
    "CH_G_204a",
    "CH_G_200a",
    "CH_G_180b",
]

PROBE_ISSUE_UKER = [
    "ER_G_276a",
    "ER_G_288a",
    "ER_G_308b",
    "ER_G_265a",
    "ER_G_310a",
    "ER_G_299a",
    "ER_G_195b",
    "ER_G_210a",
    "ER_G_295a",
    "ER_G_280a",
    "ER_G_267a",
    "ER_G_274b",
    "ER_G_278b",
    "ER_G_297a",
    "ER_G_278a",
    "ER_G_272a",
    "ER_G_303a",
    "ER_G_208a",
    "ER_G_303b",
    "ER_G_270a",
    "ER_G_298a",
    "ER_G_306a",
    "ER_G_308a",
    "ER_G_296a",
    "ER_G_285a",
    "ER_G_273a",
    "ER_G_286a",
    "ER_G_292b",
    "ER_G_196b",
    "ER_G_288b",
    "ER_G_290a",
    "ER_G_276b",
    "ER_G_266a",
    "ER_G_209a",
    "ER_G_274a",
    "ER_G_310b",
    "ER_G_292a",
    "ER_G_290b",
    "ER_G_267b",
    "ER_G_269a",
    "ER_G_281a",
    "ER_G_300a",
    "ER_G_284b",
    "ER_G_282a",
    "ER_G_294a",
    "ER_G_282b",
    "ER_G_284a",
]

# infercnv
PATH_GTF = "/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/infercnv/gencode.v32.annotation.gtf"
