"""Utilities for parsing and processing mutation data."""

import re
import pandas as pd
import requests
from typing import List, Optional, Dict, Tuple

from gbmcollaboration.utils.ensembl_utils import get_gene_coordinates

# Hardcoded TERT promoter mutation positions (single source of truth)
# Format: mutation_name -> (chromosome, position, ref, alt)
# Note: TERT is on reverse strand, so C>T mutations appear as G>A in VCF
TERT_PROMOTER_POSITIONS = {
    'C228T': ('chr5', 1295113, 'C', 'T'),  # c.−124C>T (forward strand) = G>A (reverse strand/VCF)
    'C250T': ('chr5', 1295135, 'C', 'T'),  # c.−146C>T (forward strand) = G>A (reverse strand/VCF)
}

# Reverse complement mapping for reverse strand genes
REVERSE_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# TSS-relative positions for fallback (bp upstream of TSS)
TERT_PROMOTER_TSS_OFFSETS = {
    'C228T': -124,  # TERT promoter mutation
    'C250T': -146,  # TERT promoter mutation
}


def _identify_promoter_mutation(pos: int, chrom: str, ref: str, alt: str, transcript_id: Optional[str] = None) -> Optional[str]:
    """
    Identify if a variant is a known promoter mutation.
    
    Uses two methods:
    1. Hardcoded genomic positions (primary method):
       - C228T (c.−124C>T): chr5:1295113
       - C250T (c.−146C>T): chr5:1295135
    2. TSS-relative positions (fallback method):
       - C228T = -124 bp upstream of TSS
       - C250T = -146 bp upstream of TSS
    
    Parameters
    ----------
    pos : int
        Genomic position of the variant
    chrom : str
        Chromosome (e.g., "chr5" or "5")
    ref : str
        Reference allele (e.g., "C")
    alt : str
        Alternate allele (e.g., "T")
    transcript_id : str, optional
        Ensembl transcript ID (e.g., "ENST00000310581") for TSS-relative fallback
    
    Returns
    -------
    str or None
        Mutation name (e.g., "C228T", "C250T") or None if not a known promoter mutation
    """
    # Check C>T mutations (forward strand) or G>A mutations (reverse strand/VCF)
    # TERT is on reverse strand, so C>T appears as G>A in VCF
    is_forward_strand = (ref == 'C' and alt == 'T')
    is_reverse_strand = (ref == 'G' and alt == 'A')
    
    if not (is_forward_strand or is_reverse_strand):
        return None
    
    # Method 1: Check hardcoded genomic positions (primary method)
    for mutation_name, (chrom_str, expected_pos, ref_nuc, alt_nuc) in TERT_PROMOTER_POSITIONS.items():
        # Check if this matches forward strand (C>T) or reverse strand (G>A)
        matches_forward = (ref == ref_nuc and alt == alt_nuc)
        matches_reverse = (ref == REVERSE_COMPLEMENT.get(ref_nuc) and alt == REVERSE_COMPLEMENT.get(alt_nuc))
        
        if matches_forward or matches_reverse:
            # Check if position matches (allow both chr5 and 5 format)
            if pos == expected_pos:
                if chrom == chrom_str or chrom == chrom_str.lstrip('chr') or \
                   (chrom.startswith('chr') and chrom[3:] == chrom_str.lstrip('chr')) or \
                   (not chrom.startswith('chr') and f'chr{chrom}' == chrom_str):
                    return mutation_name
    
    # Method 2: Fallback to TSS-relative positions if transcript_id is available
    if transcript_id:
        try:
            
            # Get transcript details from Ensembl
            transcript_url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?content-type=application/json"
            response = requests.get(transcript_url)
            if response.status_code == 200:
                transcript_data = response.json()
                strand = transcript_data.get('strand', 1)
                tss = transcript_data.get('start') if strand == 1 else transcript_data.get('end')
                
                # Calculate offset from TSS
                # For forward strand: upstream = decreasing coordinates, so offset = pos - tss (negative)
                # For reverse strand: upstream = increasing coordinates, so offset = tss - pos (negative)
                if strand == 1:
                    offset = pos - tss  # Negative if upstream
                else:
                    offset = tss - pos  # Negative if upstream (TERT is reverse strand)
                
                # Check if offset matches a known promoter mutation (allow ±2bp tolerance)
                for mutation_name, known_offset in TERT_PROMOTER_TSS_OFFSETS.items():
                    if abs(offset - known_offset) <= 2:
                        return mutation_name
        except Exception:
            # Silently fail - will return None
            pass
    
    return None


def _get_promoter_genomic_positions(gene_name: str, mutations_list: List[str]) -> Dict[str, Tuple[str, int]]:
    """
    Get genomic positions for promoter mutations.
    
    Uses two methods:
    1. Hardcoded genomic positions (primary method): Uses TERT_PROMOTER_POSITIONS
    2. TSS-relative positions (fallback method): Uses TERT_PROMOTER_TSS_OFFSETS
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "TERT")
    mutations_list : list[str]
        List of mutation names (e.g., ["C228T", "C250T"])
        
    Returns
    -------
    dict[str, tuple[str, int]]
        Mapping of mutation name to (chromosome, position) tuple
    """    
    positions = {}
    
    # For each mutation, first try hardcoded positions
    for mutation in mutations_list:
        if mutation in TERT_PROMOTER_POSITIONS:
            chrom_str, expected_pos, _, _ = TERT_PROMOTER_POSITIONS[mutation]
            positions[mutation] = (chrom_str, expected_pos)
            continue
        
        # Method 2: Fallback to TSS-relative calculation
        if mutation in TERT_PROMOTER_TSS_OFFSETS:
            try:
                # Get gene coordinates
                gene_coords = get_gene_coordinates(gene_name)
                if gene_coords is None:
                    continue
                
                chrom_str = gene_coords['chromosome']  # e.g., "chr5"
                
                # Get canonical transcript to find TSS
                gene_id = gene_coords['ensembl_id']
                url = f"https://rest.ensembl.org/lookup/id/{gene_id}?content-type=application/json"
                response = requests.get(url)
                if response.status_code != 200:
                    continue
                
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
                        
                        tss_offset = TERT_PROMOTER_TSS_OFFSETS[mutation]
                        
                        # Calculate genomic position
                        # For forward strand: upstream = TSS - offset (offset is negative)
                        # For reverse strand: upstream = TSS + offset (offset is negative, so we add)
                        if strand == 1:
                            expected_pos = tss + tss_offset  # tss_offset is negative, so this subtracts
                        else:
                            expected_pos = tss - tss_offset  # For reverse strand, upstream is + direction
                        
                        positions[mutation] = (chrom_str, expected_pos)
            except Exception as e:
                print(f"Warning: Could not get promoter positions for {gene_name} {mutation} from Ensembl: {e}")
    
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
        
        # Use display name with "p" suffix for promoter genes, original name for computation
        display_gene = f"{gene}p" if is_promoter else gene
        col_name = f"{display_gene}_mut"
        result_dict[col_name] = pd.Series(False, index=all_patients)
        
        if df_all_snv_indel.empty:
            continue
        
        # Filter variants for this gene
        gene_variants = df_all_snv_indel[df_all_snv_indel['gene_name'] == gene].copy()
        
        if gene_variants.empty:
            continue
        
        patients_with_mutation = set()
        if is_promoter:
            # Promoter mutations: use hardcoded genomic positions
            mutation_positions = _get_promoter_genomic_positions(gene, mutations_list)
            
            # Filter variants directly at the known positions (single filtering step)
            for mutation in mutations_list:
                if mutation in mutation_positions:
                    chrom_str, expected_pos = mutation_positions[mutation]
                    _, _, ref_nuc, alt_nuc = TERT_PROMOTER_POSITIONS[mutation]
                    
                    # TERT is on reverse strand, so C>T mutations appear as G>A in VCF
                    # Check for both forward strand (C>T) and reverse strand (G>A) representations
                    ref_rev = REVERSE_COMPLEMENT.get(ref_nuc)
                    alt_rev = REVERSE_COMPLEMENT.get(alt_nuc)
                    
                    # Filter variants at the exact position (check both strand orientations)
                    matching_variants = df_all_snv_indel[
                        (
                            (df_all_snv_indel['CHROM'] == chrom_str) |
                            (df_all_snv_indel['CHROM'] == chrom_str.lstrip('chr'))
                        ) &
                        (df_all_snv_indel['POS'] == expected_pos) &
                        (
                            # Forward strand: C>T
                            ((df_all_snv_indel['REF'] == ref_nuc) & (df_all_snv_indel['ALT'] == alt_nuc)) |
                            # Reverse strand: G>A (reverse complement)
                            ((ref_rev is not None) & (alt_rev is not None) & 
                             (df_all_snv_indel['REF'] == ref_rev) & (df_all_snv_indel['ALT'] == alt_rev))
                        )
                    ]
                    patients_with_mutation.update(matching_variants['sample'].unique())
                else:
                    # Fallback: use _identify_promoter_mutation if position not found
                    if mutation in TERT_PROMOTER_POSITIONS:
                        _, _, ref_nuc, alt_nuc = TERT_PROMOTER_POSITIONS[mutation]
                    else:
                        ref_nuc = mutation[0]
                        alt_nuc = mutation[-1]
                    
                    matching_variants = gene_variants[
                        (gene_variants['REF'] == ref_nuc) &
                        (gene_variants['ALT'] == alt_nuc)
                    ].copy()
                    
                    for _, variant_row in matching_variants.iterrows():
                        transcript_id = variant_row.get('transcript_id') if 'transcript_id' in variant_row else None
                        identified_mutation = _identify_promoter_mutation(
                            variant_row['POS'],
                            variant_row['CHROM'],
                            variant_row['REF'],
                            variant_row['ALT'],
                            transcript_id=transcript_id
                        )
                        if identified_mutation == mutation:
                            patients_with_mutation.add(variant_row['sample'])
        else:
            # Coding mutations: check amino acid changes
            for mutation in mutations_list:
                aa_mask = gene_variants['amino_acid'].str.contains(mutation, na=False, case=False)
                patients_with_mutation.update(gene_variants[aa_mask]['sample'].unique())
                
                # Also try nucleotide check as fallback (for cases like C228T format in coding region)
                # Note: If mutation is in cDNA notation for reverse strand genes, VCF will have reverse complement
                if len(mutation) >= 3 and mutation[0].isalpha() and mutation[-1].isalpha():
                    try:
                        pos_match = re.search(r'(\d+)', mutation)
                        if pos_match:
                            pos_str = pos_match.group(1)
                            ref_nuc = mutation[0]
                            alt_nuc = mutation[-1]
                            
                            # Check if gene is on reverse strand (for strand-aware nucleotide matching)
                            gene_coords = get_gene_coordinates(gene)
                            is_reverse_strand = False
                            if gene_coords:
                                try:
                                    gene_id = gene_coords['ensembl_id']
                                    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?content-type=application/json"
                                    response = requests.get(url, timeout=5)
                                    if response.status_code == 200:
                                        gene_data = response.json()
                                        strand = gene_data.get('strand', 1)
                                        is_reverse_strand = (strand == -1)
                                except:
                                    pass  # Fall back to forward strand assumption
                            
                            # Match forward strand (as reported) or reverse complement (as in VCF)
                            ref_rev = REVERSE_COMPLEMENT.get(ref_nuc)
                            alt_rev = REVERSE_COMPLEMENT.get(alt_nuc)
                            
                            if is_reverse_strand and ref_rev and alt_rev:
                                # Gene is on reverse strand: check reverse complement in VCF
                                nuc_mask = (
                                    (gene_variants['REF'] == ref_rev) &
                                    (gene_variants['ALT'] == alt_rev) &
                                    (gene_variants['POS'].astype(str).str.contains(pos_str, na=False))
                                )
                            else:
                                # Gene is on forward strand or strand unknown: check as reported
                                nuc_mask = (
                                    (gene_variants['REF'] == ref_nuc) &
                                    (gene_variants['ALT'] == alt_nuc) &
                                    (gene_variants['POS'].astype(str).str.contains(pos_str, na=False))
                                )
                            
                            patients_with_mutation.update(gene_variants[nuc_mask]['sample'].unique())
                    except:
                        pass
        
        result_dict[col_name].loc[list(patients_with_mutation)] = True
        print(f"  {display_gene}: {len(patients_with_mutation)} patients with mutations")
    
    # 2. Process amplifications (CNV duplications)
    print("\nProcessing amplifications...")
    for gene in genes_amp:
        # Check if this gene is a promoter gene in mutations list
        is_promoter = False
        if 'location' in genes_mut_df.columns:
            gene_rows = genes_mut_df[genes_mut_df['gene'] == gene]
            if not gene_rows.empty:
                location_val = str(gene_rows.iloc[0].get('location', '')).lower().strip()
                is_promoter = location_val == 'promoter'
        
        display_gene = f"{gene}p" if is_promoter else gene
        col_name = f"{display_gene}_amp"
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
        # Check if this gene is a promoter gene in mutations list
        is_promoter = False
        if 'location' in genes_mut_df.columns:
            gene_rows = genes_mut_df[genes_mut_df['gene'] == gene]
            if not gene_rows.empty:
                location_val = str(gene_rows.iloc[0].get('location', '')).lower().strip()
                is_promoter = location_val == 'promoter'
        
        display_gene = f"{gene}p" if is_promoter else gene
        col_name = f"{display_gene}_lof"
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
    # Build mapping of original gene names to display names
    gene_to_display = {}
    for _, row in genes_mut_df.iterrows():
        gene = row['gene']
        is_promoter = False
        if 'location' in genes_mut_df.columns:
            location_val = str(row.get('location', '')).lower().strip()
            is_promoter = location_val == 'promoter'
        display_gene = f"{gene}p" if is_promoter else gene
        gene_to_display[gene] = display_gene
    
    # Also check genes in amp/lof lists
    for gene in genes_amp + genes_lof:
        if gene not in gene_to_display:
            # Check if it's in mutations list
            if 'location' in genes_mut_df.columns:
                gene_rows = genes_mut_df[genes_mut_df['gene'] == gene]
                if not gene_rows.empty:
                    location_val = str(gene_rows.iloc[0].get('location', '')).lower().strip()
                    is_promoter = location_val == 'promoter'
                    display_gene = f"{gene}p" if is_promoter else gene
                    gene_to_display[gene] = display_gene
                else:
                    gene_to_display[gene] = gene
            else:
                gene_to_display[gene] = gene
    
    # Create summary columns using display names
    all_genes = set(genes_mut_df['gene'].tolist() + genes_amp + genes_lof)
    for gene in all_genes:
        display_gene = gene_to_display.get(gene, gene)
        cols_to_check = [col for col in result_df.columns if col.startswith(f"{display_gene}_")]
        if len(cols_to_check) > 1:
            # Create combined column: gene has ANY alteration
            result_df[display_gene] = result_df[cols_to_check].any(axis=1)
            print(f"  {display_gene}: Created summary column (combines {len(cols_to_check)} alteration types)")
    
    return result_df

