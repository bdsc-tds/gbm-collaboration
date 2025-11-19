"""Utilities for parsing and processing mutation data."""

import re
import pandas as pd
import requests
from typing import List, Optional, Dict

from gbmcollaboration.utils.ensembl_utils import get_gene_coordinates


def _get_promoter_genomic_positions(gene_name: str, mutations_list: List[str], df_all_snv_indel: pd.DataFrame) -> Dict[str, int]:
    """
    Get genomic positions for promoter mutations using TSS-relative positions.
    
    Promoter mutation names like "C228T" refer to positions relative to transcription start site (TSS):
    - C228T = -124 bp upstream of TSS
    - C250T = -146 bp upstream of TSS
    
    This function:
    1. Gets the canonical transcript and TSS from Ensembl
    2. Calculates genomic positions based on TSS-relative positions
    3. Verifies variants exist at those positions matching REF/ALT
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "TERT")
    mutations_list : list[str]
        List of mutation names (e.g., ["C228T", "C250T"])
    df_all_snv_indel : pd.DataFrame
        All SNV/indel variants to search in
        
    Returns
    -------
    dict[str, int]
        Mapping of mutation name to genomic position
    """    
    positions = {}
    
    # Map of known promoter mutations to their TSS-relative positions
    # Format: mutation_name -> bp upstream of TSS
    promoter_tss_offsets = {
        'C228T': -124,  # TERT promoter mutation
        'C250T': -146,  # TERT promoter mutation
    }
    
    # Get gene coordinates
    gene_coords = get_gene_coordinates(gene_name)
    if gene_coords is None:
        return positions
    
    chrom_str = gene_coords['chromosome']  # e.g., "chr5"
    
    # Get canonical transcript to find TSS
    try:
        gene_id = gene_coords['ensembl_id']
        url = f"https://rest.ensembl.org/lookup/id/{gene_id}?content-type=application/json"
        response = requests.get(url)
        if response.status_code != 200:
            return positions
        
        gene_data = response.json()
        canonical_transcript_id = gene_data.get('canonical_transcript')
        
        if canonical_transcript_id:
            # Get transcript details
            transcript_url = f"https://rest.ensembl.org/lookup/id/{canonical_transcript_id}?content-type=application/json"
            transcript_response = requests.get(transcript_url)
            if transcript_response.status_code == 200:
                transcript_data = transcript_response.json()
                strand = transcript_data.get('strand', 1)
                tss = transcript_data.get('start') if strand == 1 else transcript_data.get('end')
                
                # For each mutation, calculate genomic position from TSS
                for mutation in mutations_list:
                    if mutation not in promoter_tss_offsets:
                        continue
                    
                    tss_offset = promoter_tss_offsets[mutation]
                    ref_nuc = mutation[0]  # First character (REF)
                    alt_nuc = mutation[-1]  # Last character (ALT)
                    
                    # Calculate genomic position
                    # For forward strand: upstream = TSS - offset (offset is negative)
                    # For reverse strand: upstream = TSS + offset (offset is negative, so we add)
                    if strand == 1:
                        expected_pos = tss + tss_offset  # tss_offset is negative, so this subtracts
                    else:
                        expected_pos = tss - tss_offset  # For reverse strand, upstream is + direction
                    
                    # Verify variant exists at this position with matching REF/ALT
                    matching_variants = df_all_snv_indel[
                        (df_all_snv_indel['gene_name'] == gene_name) &
                        (df_all_snv_indel['CHROM'] == chrom_str) &
                        (df_all_snv_indel['POS'] == expected_pos) &
                        (df_all_snv_indel['REF'] == ref_nuc) &
                        (df_all_snv_indel['ALT'] == alt_nuc)
                    ]
                    
                    if not matching_variants.empty:
                        positions[mutation] = expected_pos
                    else:
                        # Try searching within a small window (±5bp) in case of coordinate differences
                        window_variants = df_all_snv_indel[
                            (df_all_snv_indel['gene_name'] == gene_name) &
                            (df_all_snv_indel['CHROM'] == chrom_str) &
                            (df_all_snv_indel['POS'] >= expected_pos - 5) &
                            (df_all_snv_indel['POS'] <= expected_pos + 5) &
                            (df_all_snv_indel['REF'] == ref_nuc) &
                            (df_all_snv_indel['ALT'] == alt_nuc)
                        ]
                        if not window_variants.empty:
                            # Use the closest position
                            closest_pos = window_variants['POS'].iloc[0]
                            positions[mutation] = closest_pos
                            print(f"  Note: {gene_name} {mutation} found at {closest_pos} (expected {expected_pos}, TSS={tss})")
                        
    except Exception as e:
        print(f"Warning: Could not get promoter positions for {gene_name} from Ensembl: {e}")
    
    return positions


def parse_mutations(mut_str: str) -> List[str]:
    """Parse mutations string like "['C228T', 'C250T']" into a list.
    
    Parameters
    ----------
    mut_str : str
        String representation of mutation list
        
    Returns
    -------
    list[str]
        List of mutation strings
    """
    if pd.isna(mut_str) or mut_str == '':
        return []
    # Remove brackets and quotes, split by comma
    mut_str = mut_str.strip().strip('[]').strip("'\"")
    if not mut_str:
        return []
    # Split by comma and clean up
    mutations = [m.strip().strip("'\"") for m in mut_str.split(',')]
    return [m for m in mutations if m]


def clean_gene_name(gene_str: str) -> Optional[str]:
    """Remove comments from gene names (everything after #).
    
    Parameters
    ----------
    gene_str : str
        Gene name string, possibly with comments
        
    Returns
    -------
    str or None
        Cleaned gene name without comments
    """
    if pd.isna(gene_str):
        return None
    return str(gene_str).split('#')[0].strip()


def create_comprehensive_patient_gene_matrix(
    genes_mut_df: pd.DataFrame,
    genes_amp: List[str],
    genes_lof: List[str],
    df_all_snv_indel: pd.DataFrame,
    df_all_cnv: pd.DataFrame
) -> pd.DataFrame:
    """
    Create a comprehensive patient X gene matrix with separate columns for mutations, 
    amplifications, and loss-of-function alterations.
    
    For genes that can have multiple alteration types (e.g., EGFR can have both mutation 
    and amplification), separate columns are created: EGFR_mut, EGFR_amp, etc.
    
    Parameters
    ----------
    genes_mut_df : pd.DataFrame
        DataFrame with columns: gene, mutations
    genes_amp : list
        List of gene names that can be amplified
    genes_lof : list
        List of gene names that can have loss-of-function
    df_all_snv_indel : pd.DataFrame
        All SNV/indel variants with columns: sample, gene_name, amino_acid
    df_all_cnv : pd.DataFrame
        All CNV variants with columns: sample, CHROM, START, END, ALT, log2_ratio
    
    Returns
    -------
    pd.DataFrame
        Patient X gene matrix with columns for each alteration type
    """
    from gbmcollaboration.utils.ensembl_utils import get_gene_coordinates
    
    # Get all unique patients
    all_patients_snv = set(df_all_snv_indel['sample'].unique()) if not df_all_snv_indel.empty else set()
    all_patients_cnv = set(df_all_cnv['sample'].unique()) if not df_all_cnv.empty else set()
    all_patients = sorted(all_patients_snv | all_patients_cnv)
    
    print(f"Found {len(all_patients)} patients")
    
    # Initialize result dictionary
    result_dict = {}
    
    # 1. Process mutations (specific SNV/indel mutations)
    print("\nProcessing mutations...")
    for _, row in genes_mut_df.iterrows():
        gene = row['gene']
        mutations_list = parse_mutations(row['mutations'])
        if not mutations_list:
            continue
        
        # Check if mutations are nucleotide-level (promoter) or amino acid-level
        # If 'location' column exists and is 'promoter', treat as nucleotide mutations
        is_promoter = False
        if 'location' in genes_mut_df.columns:
            location_val = str(row.get('location', '')).lower().strip()
            is_promoter = location_val == 'promoter'
            if is_promoter:
                print(f"  Detected {gene} as promoter mutations")
        
        col_name = f"{gene}_mut"
        result_dict[col_name] = pd.Series(False, index=all_patients)
        
        if df_all_snv_indel.empty:
            continue
        
        # Filter variants for this gene
        gene_variants = df_all_snv_indel[df_all_snv_indel['gene_name'] == gene].copy()
        
        if gene_variants.empty:
            continue
        
        patients_with_mutation = set()
        for mutation in mutations_list:
            if is_promoter:
                # Promoter mutations: identify by checking position relative to TSS using transcript_id from variant
                ref_nuc = mutation[0]  # First character
                alt_nuc = mutation[-1]  # Last character
                
                # Filter variants matching REF/ALT
                matching_variants = gene_variants[
                    (gene_variants['REF'] == ref_nuc) &
                    (gene_variants['ALT'] == alt_nuc)
                ].copy()
                
                # For each matching variant, check if it's the specific promoter mutation
                for _, variant_row in matching_variants.iterrows():
                    identified_mutation = _identify_promoter_mutation(
                        variant_row['POS'],
                        variant_row['transcript_id'],
                        variant_row['REF'],
                        variant_row['ALT']
                    )
                    if identified_mutation == mutation:
                        patients_with_mutation.add(variant_row['sample'])
            else:
                # Coding mutations: check amino acid changes
                aa_mask = gene_variants['amino_acid'].str.contains(mutation, na=False, case=False)
                patients_with_mutation.update(gene_variants[aa_mask]['sample'].unique())
                
                # Also try nucleotide check as fallback (for cases like C228T format in coding region)
                if len(mutation) >= 3 and mutation[0].isalpha() and mutation[-1].isalpha():
                    try:
                        pos_match = re.search(r'(\d+)', mutation)
                        if pos_match:
                            pos_str = pos_match.group(1)
                            ref_nuc = mutation[0]
                            alt_nuc = mutation[-1]
                            nuc_mask = (
                                (gene_variants['REF'] == ref_nuc) &
                                (gene_variants['ALT'] == alt_nuc) &
                                (gene_variants['POS'].astype(str).str.contains(pos_str, na=False))
                            )
                            patients_with_mutation.update(gene_variants[nuc_mask]['sample'].unique())
                    except:
                        pass
        
        result_dict[col_name].loc[list(patients_with_mutation)] = True
        print(f"  {gene}: {len(patients_with_mutation)} patients with mutations")
    
    # 2. Process amplifications (CNV duplications)
    print("\nProcessing amplifications...")
    for gene in genes_amp:
        col_name = f"{gene}_amp"
        result_dict[col_name] = pd.Series(False, index=all_patients)
        
        if df_all_cnv.empty:
            continue
        
        # Get gene coordinates
        gene_coords = get_gene_coordinates(gene)
        if gene_coords is None:
            print(f"  {gene}: Could not get coordinates, skipping")
            continue
        
        # Filter CNVs that overlap this gene and are duplications
        # Threshold: log2_ratio > 0.3 indicates amplification (typically ~1.2x+ copy number)
        # Reference: VisCap methodology (PMC4940431) and IDH-mutant astrocytomas study (PMC6660955)
        gene_cnvs = df_all_cnv[
            (df_all_cnv['CHROM'] == gene_coords['chromosome']) &
            (df_all_cnv['START'] <= gene_coords['end']) &
            (df_all_cnv['END'] >= gene_coords['start']) &
            (df_all_cnv['ALT'].str.contains('DUP', na=False, case=False)) &
            (df_all_cnv['log2_ratio'] > 0.3)  # Amplification threshold: > 0.3 (~1.2x+ copy number)
        ]
        
        patients_with_amp = gene_cnvs['sample'].unique()
        result_dict[col_name].loc[list(patients_with_amp)] = True
        print(f"  {gene}: {len(patients_with_amp)} patients with amplification")
    
    # 3. Process loss-of-function (CNV deletions or truncating mutations)
    print("\nProcessing loss-of-function...")
    for gene in genes_lof:
        col_name = f"{gene}_lof"
        result_dict[col_name] = pd.Series(False, index=all_patients)
        
        patients_with_lof = set()
        
        # Check for CNV deletions
        if not df_all_cnv.empty:
            gene_coords = get_gene_coordinates(gene)
            if gene_coords:
                gene_cnvs = df_all_cnv[
                    (df_all_cnv['CHROM'] == gene_coords['chromosome']) &
                    (df_all_cnv['START'] <= gene_coords['end']) &
                    (df_all_cnv['END'] >= gene_coords['start']) &
                    (df_all_cnv['ALT'].str.contains('DEL', na=False, case=False)) &
                    (df_all_cnv['log2_ratio'] < -0.3)  # Deletion threshold: < -0.3 (~0.7x copy number)
                ]
                patients_with_lof.update(gene_cnvs['sample'].unique())
        
        # Check for truncating mutations (stop_gained, frameshift)
        if not df_all_snv_indel.empty:
            gene_variants = df_all_snv_indel[df_all_snv_indel['gene_name'] == gene].copy()
            truncating_mask = gene_variants['variant_type'].str.contains(
                'stop_gained|frameshift', na=False, case=False, regex=True
            )
            patients_with_lof.update(gene_variants[truncating_mask]['sample'].unique())
        
        result_dict[col_name].loc[list(patients_with_lof)] = True
        print(f"  {gene}: {len(patients_with_lof)} patients with loss-of-function")
    
    # Create final DataFrame
    result_df = pd.DataFrame(result_dict)
    
    # Add summary columns for genes with multiple alteration types
    print("\nCreating summary columns...")
    genes_mut_list = genes_mut_df['gene'].tolist()
    all_genes = set(genes_mut_list + genes_amp + genes_lof)
    for gene in all_genes:
        cols_to_check = [col for col in result_df.columns if col.startswith(f"{gene}_")]
        if len(cols_to_check) > 1:
            # Create combined column: gene has ANY alteration
            result_df[gene] = result_df[cols_to_check].any(axis=1)
            print(f"  {gene}: Created summary column (combines {len(cols_to_check)} alteration types)")
    
    return result_df

