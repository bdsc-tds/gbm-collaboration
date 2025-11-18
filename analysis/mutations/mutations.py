# %%
import pandas as pd
import joblib
import numpy as np
from tqdm import tqdm
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
# NOTE: How to identify mutation types from VCF CSQ field:
# The VCF CSQ (Consequence) field contains variant consequences from VEP (Variant Effect Predictor).
# Common values in fields[1] (Consequence field) include:
# - "missense_variant" -> maps to "missense"
# - "stop_gained" or "stop_lost" -> maps to "truncation" 
# - "splice_acceptor_variant" or "splice_donor_variant" -> maps to "splice"
# - "frameshift_variant" -> maps to "frameshift"
# - "inframe_deletion" or "inframe_insertion" -> maps to "inframe_deletion" or "inframe_insertion"
# - "synonymous_variant" -> maps to "synonymous"
# - "splice_region_variant" -> can be considered "splice"
#
# The extract_variants() function extracts fields[1] as "variant_type" which contains these values.
# You can filter variants by checking if variant_type contains the mutation type you're looking for.
# For example: df[df['variant_type'].str.contains('missense_variant', na=False)]

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
# FILTERING FUNCTIONS - Use these to filter the comprehensive tables
print("\n" + "=" * 80)
print("STEP 2: Filtering functions available")
print("=" * 80)
print("""
You can now filter the comprehensive tables using:

1. Filter by driver genes:
   df_filtered = df_all_snv_indel[df_all_snv_indel['gene_name'].isin(driver_genes)]

2. Filter by mutation type (from variant_type column):
   df_missense = df_all_snv_indel[df_all_snv_indel['variant_type'].str.contains('missense_variant', na=False)]
   df_truncation = df_all_snv_indel[df_all_snv_indel['variant_type'].str.contains('stop_gained', na=False)]
   df_splice = df_all_snv_indel[df_all_snv_indel['variant_type'].str.contains('splice', na=False, case=False)]
   df_frameshift = df_all_snv_indel[df_all_snv_indel['variant_type'].str.contains('frameshift', na=False)]

3. Filter by specific mutations (from amino_acid column):
   df_braf_v600e = df_all_snv_indel[
       (df_all_snv_indel['gene_name'] == 'BRAF') & 
       (df_all_snv_indel['amino_acid'].str.contains('V600E', na=False))
   ]

4. Filter CNVs by gene overlap (requires gene coordinates):
   # First get gene coordinates, then filter CNVs that overlap
""")

# %%
# Create comprehensive patient X gene matrix
print("\n" + "=" * 80)
print("STEP 3: Creating comprehensive patient X gene alteration matrix")
print("=" * 80)
print("""
This matrix includes:
- {gene}_mut: Specific driver mutations (SNV/indel)
- {gene}_amp: Gene amplifications (CNV duplications)
- {gene}_lof: Loss-of-function alterations (CNV deletions or truncating mutations)
- {gene}: Summary column (True if patient has ANY alteration type for that gene)

For example, EGFR can have both EGFR_mut and EGFR_amp columns if a patient has
both a mutation and an amplification.
""")

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
# Plot number of patients with True per column
if not patient_gene_matrix.empty:
    # Calculate frequencies
    frequencies = patient_gene_matrix.sum().sort_values(ascending=True)
    
    # Create plot
    plt.figure(figsize=(12, max(8, len(frequencies) * 0.3)))
    frequencies.plot(kind='barh')
    plt.xlabel('Number of Patients')
    plt.ylabel('Alteration Type')
    plt.title('Number of Patients with Driver Alterations')
    plt.tight_layout()
    plt.show()
    
    print(f"\nTotal alterations found: {frequencies.sum()}")
    print(f"Columns with at least one patient: {(frequencies > 0).sum()} / {len(frequencies)}")

# %%
# Create oncoplot-style heatmap
create_oncoplot(patient_gene_matrix, figsize=(20, max(8, len(patient_gene_matrix.columns) * 0.3)))
