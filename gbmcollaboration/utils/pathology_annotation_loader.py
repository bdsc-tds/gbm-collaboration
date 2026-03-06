"""Util function to list, read and merge annotation files for each samples from S3 bucket."""

# Packages
import boto3
import pandas as pd
from io import BytesIO
import os
import re


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
