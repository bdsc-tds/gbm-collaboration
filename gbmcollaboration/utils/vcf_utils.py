"""Utilities for parsing and processing VCF files."""

import boto3
import os
import pandas as pd
from typing import Optional, List, Dict


def get_csq_field_index(fn_vcf: str) -> tuple[int, int]:
    """Get CSQ field indices from VCF header.
    
    Parameters
    ----------
    fn_vcf : str
        Path to VCF file (can be S3 path)
        
    Returns
    -------
    tuple[int, int]
        (aa_coord_index, aa_change_index) - indices for protein position and amino acid change
    """
    bucket = fn_vcf.lstrip("s3://").split('/')[0]
    key = '/'.join(fn_vcf.lstrip("s3://").split('/')[1:])
    s3 = boto3.client('s3')
    lines = s3.get_object(Bucket=bucket, Key=key)['Body'].iter_lines()
    csq_header = [
        line.decode('utf-8') for line in lines
        if line.decode('utf-8').startswith('##INFO=<ID=CSQ')
    ][0]
    description = csq_header.split('Description="')[1].split('"')[0]
    fields = description.split('|')
    aa_coord_index = fields.index('Protein_position')
    aa_change_index = fields.index('Amino_acids')
    
    return aa_coord_index, aa_change_index


def extract_variants(csq_field: str, aa_coord_index: int, aa_change_index: int) -> List[Dict[str, str]]:
    """Extracts the variant type for the canonical transcript from a CSQ field.
    
    Parameters
    ----------
    csq_field : str
        CSQ annotation field from VCF INFO column
    aa_coord_index : int
        Index of protein position in CSQ fields
    aa_change_index : int
        Index of amino acid change in CSQ fields
        
    Returns
    -------
    list[dict]
        List of variant dictionaries with gene_name, gene_id, transcript_id, 
        variant_type, amino_acid
    """
    annotations = csq_field.split(',')
    
    variants = []
    for annotation in annotations:
        fields = annotation.split('|')
        if len(fields) < 7:
            continue
        
        aa_coord = fields[aa_coord_index]
        aa_change = fields[aa_change_index]
        aa_orig = aa_change.split("/")[0]
        aa_mut = aa_change.split("/")[1] if '/' in aa_change else ''
        aa_variant = f"{aa_orig}{aa_coord}{aa_mut}"
        
        variants.append({
            "gene_name": fields[3],
            "gene_id": fields[4],
            "transcript_id": fields[6],
            "variant_type": fields[1],
            "amino_acid": aa_variant,
        })
    
    return variants


def extract_cnvs(df: pd.DataFrame, sample: str, gene_coords: Dict[str, any]) -> pd.DataFrame:
    """Extract CNVs from a VCF DataFrame for a given gene.
    
    Parameters
    ----------
    df : pd.DataFrame
        VCF DataFrame with columns: CHROM, POS, INFO, FILTER, and sample column
    sample : str
        Sample name (column name in df)
    gene_coords : dict
        Dictionary with chromosome, start, end keys
        
    Returns
    -------
    pd.DataFrame
        DataFrame with CNVs overlapping the gene
    """
    df['START'] = df['POS'] + 1
    df['END'] = df['INFO'].str.split('END=').str[1].str.split(';').str[0].astype(int)
    if df['INFO'].str.contains('FOLD_CHANGE_LOG').any():
        df['log2_ratio'] = df['INFO']\
            .str.split('FOLD_CHANGE_LOG=').str[1].str.split(';').str[0].astype(float)
    else:
        df['log2_ratio'] = df[sample].str.split(':').str[1].astype(float)
    df['width'] = df['END'] - df['START']

    df['gene_overlap'] = None
    df.loc[
        (df['CHROM'] == gene_coords['chromosome']) & 
        (df['START'] <= gene_coords['start']) & 
        (df['END'] >= gene_coords['end']),
        'gene_overlap'
    ] = 'full'
    df.loc[
        (df['CHROM'] == gene_coords['chromosome']) &
        ((
            (df['START'] >= gene_coords['start']) &
            (df['START'] <= gene_coords['end'])
        ) | (
            (df['END'] >= gene_coords['start']) &
            (df['END'] <= gene_coords['end'])
        )),
        'gene_overlap'
    ] = 'partial'

    gene_cnvs = df.loc[
        (df['gene_overlap'].notnull()) & (df['FILTER'].isin(['PASS', '.'])),
        ['ALT', 'FILTER', 'log2_ratio', 'width', 'gene_overlap']
    ]

    return gene_cnvs


def read_vcf_all_genes(fn_vcf: str, vcf_type: str) -> pd.DataFrame:
    """Read a VCF file and return ALL variants (all genes) as a DataFrame.
    
    Parameters
    ----------
    fn_vcf : str
        Path to VCF file (can be S3 path)
    vcf_type : str
        Type of VCF: "snv_indel" or "cnv"
        
    Returns
    -------
    pd.DataFrame
        DataFrame with all variants from the VCF file
    """
    sample = '_'.join(os.path.basename(fn_vcf).split('_')[:3])
    df_vcf = pd.read_csv(fn_vcf, sep='\t', comment='#', header=None)
    df_vcf.columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample
    ]

    if vcf_type == "snv_indel":
        # Extract all variants from CSQ field (all genes)
        aa_coord_index, aa_change_index = get_csq_field_index(fn_vcf)
        
        # Create list to store all variant records with VCF info
        all_variants = []
        for idx, row in df_vcf.iterrows():
            variants = extract_variants(row["INFO"], aa_coord_index, aa_change_index)
            for variant in variants:
                variant["sample"] = sample
                variant["CHROM"] = row["CHROM"]
                variant["POS"] = row["POS"]
                variant["REF"] = row["REF"]
                variant["ALT"] = row["ALT"]
                variant["FILTER"] = row["FILTER"]
                variant["QUAL"] = row["QUAL"]
                all_variants.append(variant)
        
        if not all_variants:
            return pd.DataFrame()
        
        df = pd.DataFrame(all_variants)
        df = df[["sample", "variant_type", "amino_acid", "gene_name", "gene_id",
                 "transcript_id", "CHROM", "POS", "REF", "ALT", "FILTER", "QUAL"]]
        return df

    elif vcf_type == "cnv":
        # Extract all CNVs (we'll need to match to genes later)
        df = df_vcf.copy()
        df['START'] = df['POS'] + 1
        df['END'] = df['INFO'].str.split('END=').str[1].str.split(';').str[0].astype(int)
        if df['INFO'].str.contains('FOLD_CHANGE_LOG').any():
            df['log2_ratio'] = df['INFO']\
                .str.split('FOLD_CHANGE_LOG=').str[1].str.split(';').str[0].astype(float)
        else:
            df['log2_ratio'] = df[sample].str.split(':').str[1].astype(float)
        df['width'] = df['END'] - df['START']
        df["sample"] = sample
        df = df[["sample", "CHROM", "START", "END", "ALT", "FILTER", "log2_ratio", "width"]]
        # Only keep PASS variants
        df = df[df['FILTER'].isin(['PASS', '.'])]
        return df
    
    return pd.DataFrame()


def extract_gene_variants(
    df: pd.DataFrame, 
    gene_name: str, 
    transcript_id: Optional[str] = None, 
    canonical_transcript_only: bool = False
) -> pd.DataFrame:
    """Extract variants for a specific gene from a DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with variant data
    gene_name : str
        Gene symbol to filter for
    transcript_id : str, optional
        Specific transcript ID to filter for
    canonical_transcript_only : bool
        If True, automatically get canonical transcript ID
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with variants for the specified gene
    """
    if df.empty:
        return df
    samples = df["sample"].unique()
    df_gene = df[df["gene_name"] == gene_name]
    if df_gene.empty:
        return df_gene
    if canonical_transcript_only:
        from gbmcollaboration.utils.ensembl_utils import get_canonical_transcript
        gene_id = df_gene["gene_id"].iloc[0]
        transcript_id = get_canonical_transcript(gene_id).split('.')[0]
    if transcript_id:
        df_gene = df_gene[df_gene["transcript_id"] == transcript_id]

    return df_gene

