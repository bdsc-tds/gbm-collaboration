"""Utilities for loading WES (Whole Exome Sequencing) data from S3."""

import os
import pandas as pd
from typing import Optional, Tuple, Dict, List
import s3fs

from gbmcollaboration.definitions import FOLDER_VERSION
from gbmcollaboration.utils.env_utils import get_manifest


def glob_s3_bucket(dataset_id: str, center_name: str) -> Tuple[Optional[Dict[str, List[str]]], Optional[List[str]]]:
    """Glob files in an S3 bucket for a specific center.
    
    Parameters
    ----------
    dataset_id : str
        Dataset ID from manifest
    center_name : str
        Name of the center (for error messages)
        
    Returns
    -------
    tuple[dict, list] or tuple[None, None]
        (filenames_dict, samples_list) or (None, None) if error
    """
    try:
        # Get bucket and prefix from manifest
        manifest = get_manifest()
        if dataset_id not in manifest.get("datasets", {}):
            return None, None
        
        s3_info = manifest["datasets"][dataset_id]
        if isinstance(s3_info, list) and len(s3_info) == 1:
            s3_info = s3_info[0]
        
        bucket = s3_info["bucket"]
        prefix = s3_info["prefix"]
        
        s3 = s3fs.S3FileSystem(anon=False)
        file_suffixes = {
            "snv_indel_vcf": f"{FOLDER_VERSION['wes']}/sample_level/snv_indel/*_vep_status_annotated_pass_somatic.vcf",
            "cnv_vcf": f"{FOLDER_VERSION['wes']}/sample_level/cnv/*_annotated.vcf",
        }
        filenames = {
            key: [
                "s3://" + fn for fn in s3.glob(
                    f"s3://{bucket}/{prefix}{suffix}"
                )
            ] for key, suffix in file_suffixes.items()
        }
        samples = [
            os.path.basename(fn).rsplit("_vep_status_annotated_pass_somatic.vcf", 1)[0]
            for fn in filenames["snv_indel_vcf"]
        ]
        
        return filenames, samples
    except Exception as e:
        print(f"Error globbing S3 bucket for {center_name}: {e}")
        return None, None


def read_binary_status_all_genes(fn: str) -> pd.DataFrame:
    """Read the binary status CSV file and return ALL genes as a DataFrame.
    
    Parameters
    ----------
    fn : str
        Path to binary status CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with binary status for all genes
    """
    df = pd.read_csv(fn)
    df["sample"] = os.path.basename(fn).rsplit("_binary_status.csv", 1)[0]
    df = df.rename(columns={
        "deletion": "binary_deletion",
        "duplication": "binary_duplication",
        "oncogenic": "binary_oncogenic"
    })
    if "ensembl_gene_id" in df.columns:
        df = df.drop(columns=["ensembl_gene_id"])
    return df


def read_binary_status(fn: str, gene_name: str) -> pd.DataFrame:
    """Read the binary status CSV file and return a DataFrame filtered for a gene.
    
    Parameters
    ----------
    fn : str
        Path to binary status CSV file
    gene_name : str
        Gene name to filter for
        
    Returns
    -------
    pd.DataFrame
        DataFrame with binary status for the specified gene
    """
    df = read_binary_status_all_genes(fn)
    df = df[df['gene_name'] == gene_name]
    return df

