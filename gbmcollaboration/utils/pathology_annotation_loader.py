"""Util function to list, read and merge annotation files for each samples from S3 bucket."""

# Packages
import boto3
import pandas as pd
from io import BytesIO
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

@dataclass(frozen=True)
class DataConfig:
    s3_bucket: str
    s3_prefixes: Tuple[str, ...]
    ecm_cluster_path: str ="/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/Analysis_Repo/Notebooks/ecm_neighborhood_annotations_v4.csv"
    sample_col: str = "sample_name"
    barcode_col: str = "Barcode"
    annot_col: str = "pathology_annotations"
    ecm_col: str = "ECM_neighborhood"
    spot_id_col: str = "spot_id"


@dataclass(frozen=True)
class AnalysisConfig:
    # Filtering
    min_spots_per_annotation: int = 0  # 0 keeps all
    keep_only_shared_annotations: bool = True
    shared_annotations: Tuple[str, ...] = (
        "blood vessels",
        "cortex with tumor infiltration",
        "cortex without visible tumor infiltration",
        "infiltration zone",
        "necrosis",
        "tumor",
    )

    # plotting
    clip_residuals: float = 3.0
    q_threshold: float = 0.05


# Function
def concat_pathology_annotations_from_s3(bucket: str, prefix: str, **read_csv_kwargs) -> pd.DataFrame:
    """
    Concatenate pathology annotation CSV files stored in an S3 bucket, using flexible filename matching.

    This function scans an S3 bucket and prefix for pathology annotation CSV files,
    supporting multiple filename conventions (including legacy and variant formats).

    A file is considered valid if:
        - It ends with '.csv'
        - AND matches one of the known regex patterns
        - OR contains both 'pathology' and 'annot' in its name (case-insensitive)

    For each valid file:
        - The sample name is extracted from the filename using tolerant heuristics
        - The CSV is loaded into a pandas DataFrame
        - Columns named 'Default', 'Necrosis', or 'Untitled' are normalized to 'default'
        - Metadata columns 'sample_name' and 's3_key' are added

    Parameters
    ----------
    bucket : str
        Name of the S3 bucket.

    prefix : str
        Prefix (folder path) under which files are searched.

    **read_csv_kwargs :
        Additional keyword arguments passed directly to pandas.read_csv().

    Returns
    -------
    pandas.DataFrame
        A concatenated DataFrame containing all matched annotation files.

    Raises
    ------
    ValueError
        If no matching pathology annotation files are found.
    """

    s3 = boto3.client("s3")

    # Known filename patterns
    # Accept any separator between tokens (_ or .)
    SEP = r"[._]"

    p_visium = re.compile(
        rf"^Visium{SEP}(.+?){SEP}annot{SEP}patho\.csv$",
        re.IGNORECASE,
    )
    p_simple = re.compile(
        rf"^(.+?){SEP}?pathology{SEP}?annotations?\.csv$",
        re.IGNORECASE,
    )
    p_cloupe = re.compile(
        rf"^(.+?){SEP}?cloupe.*pathology{SEP}?annotations?\.csv$",
        re.IGNORECASE,
    )

    def is_generic_pathology_annotation(name: str) -> bool:
        n = name.lower()
        return n.endswith(".csv") and ("pathology" in n) and ("annot" in n)

    def extract_sample_name(filename: str) -> str:
        """
        Extract sample name from any pathology annotation filename.

        Strategy:
        - Remove extension
        - Remove trailing 'pathology_annotation(s)'
        - Remove trailing technical suffixes (like .cloupe, _cloupe, etc.)
        - Return the leading part
        """

        base = os.path.splitext(os.path.basename(filename))[0]

        base = re.sub(
            r"[._]?pathology[._]?annotations?$",
            "",
            base,
            flags=re.IGNORECASE,
        )
        base = re.sub(
            r"[._]?annot_patho.*$",
            "",
            base,
            flags=re.IGNORECASE,
        )
        base = re.sub(
            r"[._]?cloupe.*$",
            "",
            base,
            flags=re.IGNORECASE,
        )
        base = re.sub(
            r"Visium_",
            "",
            base,
            flags=re.IGNORECASE,
        )

        base = re.sub(r"[._]+$", "", base)

        return base

    def should_take_file(filename: str) -> bool:
        if not filename.lower().endswith(".csv"):
            return False
        if p_visium.match(filename) or p_cloupe.match(filename) or p_simple.match(filename):
            return True
        return is_generic_pathology_annotation(filename)

    dfs = []
    paginator = s3.get_paginator("list_objects_v2")

    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            filename = os.path.basename(key)

            if not should_take_file(filename):
                continue

            sample_name = extract_sample_name(filename)

            resp = s3.get_object(Bucket=bucket, Key=key)
            df = pd.read_csv(BytesIO(resp["Body"].read()), **read_csv_kwargs)

            rename_map = {
                c: "pathology_annotations"
                for c in df.columns
                if c in ["default", "Default", "Necrosis", "Untitled", "Pathology_annotations"]
            }
            if rename_map:
                df.rename(columns=rename_map, inplace=True)

            df["sample_name"] = sample_name
            df["s3_key"] = key

            dfs.append(df)

    if not dfs:
        raise ValueError(
            "No pathology annotation CSV files were found under the specified bucket/prefix."
        )

    return pd.concat(dfs, ignore_index=True)



def load_annotations_from_s3(bucket: str, prefixes: Iterable[str]) -> pd.DataFrame:
    dfs = []
    for p in prefixes:
        dfs.append(concat_pathology_annotations_from_s3(bucket=bucket, prefix=p))
    out = pd.concat(dfs, ignore_index=False)
    return out


def standardize_annotations(df: pd.DataFrame, annot_col: str) -> pd.DataFrame:
    """Lowercase + fix known typos/variants from your script."""
    out = df.copy()

    out[annot_col] = out[annot_col].astype(str).str.lower()

    mapping = {
        "necrose": "necrosis",
        "white matter with tumor" : "infiltration zone",
        "white matter with tumor " : "infiltration zone",
        "infiltration zone with tumor" : "infiltration zone",
        "blood vessel" : "blood vessels",
        "cortex without visible tumor cell infiltration" : "cortex without visible tumor infiltration",
        "cortex" : "cortex without visible tumor infiltration",
        "gray matter (cortex) with tumor infiltration" : "cortex with tumor infiltration",
        "non-tumoral (gray matter)" : "cortex without visible tumor infiltration",
        "tumour (microvascular proliferation)" : "microvascular proliferation",
        "tumor (microvascular prolifeation)" : "microvascular proliferation",
        "tumor (mircovascular proliferation)" : "microvascular proliferation",
        "tumor (microvascular proliferations)" : "microvascular proliferation",
        "tumor (microvacular proliferation)" : "microvascular proliferation",
        "tumor (microvascular proliferation)" : "microvascular proliferation",
        "tumor (with microvascular proliferation)" : "microvascular proliferation",
        "tumor - divergent morphology - partially necrotic?" : "necrosis",
        "tumor (necrotic)" : "necrosis",
        "blood with necrosis" : "necrosis",
        "tumor (gliosarcoma phenotype)" : "tumor",
        "tumor (high cellularity)" : "tumor",
        "tumor - component with divergent morphology" : "tumor",
        "tumor (inflammatory)" : "tumor",
        "tumor (gliosarcoma phenotype)" : "tumor",
        "tumor (lymphocyte-rich)" : "tumor",
        "tumor (morphologically distinct component)" : "tumor",
        "tumor (with cautery artifacts)" : "tumor",
        "blood only" : "hemorrhage",
        "tumor (hemorrhage)" : "hemorrhage",
        "tumor with fibrosis" : "fibrosis",
        "tumor with fibrotic changes" : "fibrosis",
        "fibrous tissue" : "fibrosis",
        "infiltration zone with fibrosis" : "fibrosis",
        "connective tissue" : "other",
        "white matter" : "other",
        "intratumoral lymphocytic infiltrate" : "other"
    }
    out[annot_col] = out[annot_col].replace(mapping)
    return out


def exclude_bad_spots(df: pd.DataFrame, annot_col: str) -> pd.DataFrame:
    """Exclude folds/shifted/no tissue items from your script."""
    bad = {
        "tumor (section fold)",
        "tumor (attention: spots seem shifted compared to the h&e)",
        "no tissue",
        "empty (detachment, maybe tumor in vis)",
        "potentially tumor in vis (fold in he) ",
        "tumor (attention: spots seem shifted compared to h&e ) on h&e the smaller part is superimposed on the bigger one",
        "no tissue (attention: spots seem shifted compared to h&e)",
        "potentially tumor in vis (fold in he)",
    }
    out = df[~df[annot_col].isin(bad)].copy()
    out = out.dropna(subset=[annot_col])
    return out

def filter_to_shared_annotations(df: pd.DataFrame, annot_col: str, shared: Iterable[str]) -> pd.DataFrame:
    shared = set(shared)
    return df[df[annot_col].isin(shared)].copy()


def build_spot_id_and_index(df: pd.DataFrame, sample_col: str, barcode_col: str, spot_id_col: str) -> pd.DataFrame:
    out = df.copy()

    # Your one-off correction
    out[sample_col] = out[sample_col].replace({"CH_G_1789c_vis": "CH_G_189c_vis"})

    out[spot_id_col] = out[sample_col].astype(str) + "_" + out[barcode_col].astype(str)
    out = out.set_index(spot_id_col)
    return out

def add_center_from_sample(df: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    mapping = {"E": "UKER", "C": "CHUV"}
    out = df.copy()
    first = out[sample_col].astype(str).str[0].fillna("")
    out["center"] = first.map(mapping).fillna("Other")
    return out

def load_and_prepare(cfg: DataConfig, acfg: AnalysisConfig) -> pd.DataFrame:
    df_annot = load_annotations_from_s3(cfg.s3_bucket, cfg.s3_prefixes)
    df_annot = build_spot_id_and_index(df_annot, cfg.sample_col, cfg.barcode_col, cfg.spot_id_col)
    df_annot = standardize_annotations(df_annot, cfg.annot_col)
    df_annot = exclude_bad_spots(df_annot, cfg.annot_col)

    # filter nan and arachnoid values in annotation
    df_annot = df_annot[~df_annot["pathology_annotations"].isin(["nan", "arachnoid"])]

    if acfg.keep_only_shared_annotations:
        df_annot = filter_to_shared_annotations(df_annot, cfg.annot_col, acfg.shared_annotations)

    return df_annot
