# %%
import pandas as pd
import joblib
import numpy as np
from tqdm import tqdm
from pathlib import Path
import matplotlib.pyplot as plt

from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter
from gbmcollaboration.definitions import DATA_ID_CHUV, DATA_ID_UKER, DATA_ID_CR
from gbmcollaboration.utils.vcf_utils import read_vcf_all_genes
from gbmcollaboration.utils.mutation_utils import (
    parse_mutations,
    clean_gene_name,
    create_comprehensive_patient_gene_matrix,
)
from gbmcollaboration.utils.wes_io_utils import glob_s3_bucket
from gbmcollaboration.viz.mutations import create_oncoplot

# %%
RESULTS_DIR = Path('/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/WES/mutations')
RESULTS_DIR.mkdir(parents=True,exist_ok=True)

# %%
# Load driver gene lists
genes_mut_df = pd.read_csv("genes_mut.csv")
genes_amp_df = pd.read_csv("genes_amp.csv", header=None, names=['gene'])
genes_lof_df = pd.read_csv("genes_lof.csv", header=None, names=['gene'])

# Clean gene names and parse mutations
genes_amp_df['gene'] = genes_amp_df['gene'].apply(clean_gene_name)
genes_lof_df['gene'] = genes_lof_df['gene'].apply(clean_gene_name)
genes_amp = genes_amp_df['gene'].dropna().tolist()
genes_lof = genes_lof_df['gene'].dropna().tolist()

genes_mut_df['mutations_parsed'] = genes_mut_df['mutations'].apply(parse_mutations)
genes_mut = genes_mut_df['gene'].tolist()

print(f"Loaded gene lists:")
print(f"  Mutations: {len(genes_mut)} genes")
print(f"  Amplifications: {len(genes_amp)} genes")
print(f"  Loss-of-function: {len(genes_lof)} genes")
print(f"\nGenes with mutations: {', '.join(genes_mut)}")
print(f"\nGenes with amplifications: {', '.join(genes_amp)}")
print(f"\nGenes with loss-of-function: {', '.join(genes_lof)}")

# %%
# Load ALL variants from ALL genes into comprehensive tables
print("=" * 80)
print("STEP 1: Loading ALL variants from ALL genes into comprehensive tables")
print("=" * 80)

all_snv_indel_all_genes = []
all_cnv_all_genes = []

centers = {
    "mosaic_chuv": (CHUVDataCenter, DATA_ID_CHUV["CHUV_wes"]),
    "mosaic_uker": (UKERDataCenter, DATA_ID_UKER["UKER_wes"]),
    "mosaic_cr": (CRDataCenter, DATA_ID_CR["CR_wes"]),
}

for center_name, (center_obj, wes_id) in centers.items():
    print(f"\nLoading all variants from {center_name}...")

    # Get VCF file paths
    fns, samples = glob_s3_bucket(wes_id, center_name)
    if fns is None or not samples:
        print(f"  No WES data found for {center_name}")
        continue

    print(f"  Found {len(samples)} samples")

    # Load ALL SNV/indel variants (all genes)
    if fns.get("snv_indel_vcf"):
        print(f"  Loading SNV/indel variants from {len(fns['snv_indel_vcf'])} VCF files...")
        df_snv_indel = pd.concat(
            joblib.Parallel(n_jobs=-1)(
                joblib.delayed(read_vcf_all_genes)(fn_vcf, "snv_indel")
                for fn_vcf in tqdm(fns["snv_indel_vcf"], desc=f"    {center_name} SNV/indel")
            ),
            ignore_index=True
        )
        df_snv_indel["center"] = center_name
        all_snv_indel_all_genes.append(df_snv_indel)
        print(f"    Loaded {len(df_snv_indel)} variant records")

    # Load ALL CNV variants (all genes)
    if fns.get("cnv_vcf"):
        print(f"  Loading CNV variants from {len(fns['cnv_vcf'])} VCF files...")
        df_cnv = pd.concat(
            joblib.Parallel(n_jobs=-1)(
                joblib.delayed(read_vcf_all_genes)(fn_vcf, "cnv")
                for fn_vcf in tqdm(fns["cnv_vcf"], desc=f"    {center_name} CNV")
            ),
            ignore_index=True
        )
        df_cnv["center"] = center_name
        all_cnv_all_genes.append(df_cnv)
        print(f"    Loaded {len(df_cnv)} CNV records")

# Combine all centers
print("\nCombining all centers...")
df_all_snv_indel = pd.concat(all_snv_indel_all_genes, ignore_index=True) if all_snv_indel_all_genes else pd.DataFrame()
df_all_cnv = pd.concat(all_cnv_all_genes, ignore_index=True) if all_cnv_all_genes else pd.DataFrame()

print(f"\nSummary:")
print(f"  SNV/indel variants: {len(df_all_snv_indel)} records")
print(f"  CNV variants: {len(df_all_cnv)} records")

# %%
# Filter TERT promoter mutations using hardcoded genomic positions
print("\n" + "=" * 80)
print("STEP 2: Filtering TERT promoter mutations (hardcoded positions)")
print("=" * 80)

# Hardcoded TERT promoter mutation positions
TERT_PROMOTER_POSITIONS = {
    'C228T': {'chrom': 'chr5', 'pos': 1295113, 'ref': 'C', 'alt': 'T'},  # c.−124C>T
    'C250T': {'chrom': 'chr5', 'pos': 1295135, 'ref': 'C', 'alt': 'T'},  # c.−146C>T
}

if not df_all_snv_indel.empty:
    print("\nFiltering TERT promoter mutations using hardcoded positions...")
    tert_promoter_variants = []
    
    for mutation_name, pos_info in TERT_PROMOTER_POSITIONS.items():
        # Filter variants matching hardcoded position and REF/ALT
        matching = df_all_snv_indel[
            (
                (df_all_snv_indel['CHROM'] == pos_info['chrom']) |
                (df_all_snv_indel['CHROM'] == pos_info['chrom'].lstrip('chr'))
            ) &
            (df_all_snv_indel['POS'] == pos_info['pos']) &
            (df_all_snv_indel['REF'] == pos_info['ref']) &
            (df_all_snv_indel['ALT'] == pos_info['alt'])
        ]
        
        if not matching.empty:
            matching = matching.copy()
            matching['tert_mutation'] = mutation_name
            tert_promoter_variants.append(matching)
            print(f"  Found {len(matching)} variants for {mutation_name} at {pos_info['chrom']}:{pos_info['pos']}")
            print(f"    Patients: {matching['sample'].nunique()}")
        else:
            print(f"  No variants found for {mutation_name} at {pos_info['chrom']}:{pos_info['pos']}")
    
    if tert_promoter_variants:
        df_tert_promoter = pd.concat(tert_promoter_variants, ignore_index=True)
        print(f"\nTotal TERT promoter variants found: {len(df_tert_promoter)}")
        print(f"Unique patients with TERT promoter mutations: {df_tert_promoter['sample'].nunique()}")
    else:
        print("\nNo TERT promoter mutations found using hardcoded positions.")
        df_tert_promoter = pd.DataFrame()
else:
    print("No SNV/indel data available.")
    df_tert_promoter = pd.DataFrame()

# %%
# Create comprehensive patient X gene matrix
print("\n" + "=" * 80)
print("STEP 3: Creating comprehensive patient X gene alteration matrix")
print("=" * 80)

patient_gene_matrix = create_comprehensive_patient_gene_matrix(
    genes_mut_df, genes_amp, genes_lof, df_all_snv_indel, df_all_cnv
)

print(f"\nPatient X gene matrix shape: {patient_gene_matrix.shape}")
print(f"\nColumns: {list(patient_gene_matrix.columns)}")
print(f"\nFirst few rows:")
print(patient_gene_matrix.head(10))
print(f"\nAlteration frequency per column:")
print(patient_gene_matrix.sum().sort_values(ascending=False))

# %%
# Create oncoplot-style heatmap
create_oncoplot(patient_gene_matrix, figsize=(20, max(8, len(patient_gene_matrix.columns) * 0.3)))
