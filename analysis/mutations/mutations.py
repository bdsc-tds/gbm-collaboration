# %%
# Autoreload extension for IPython/Jupyter (safe to run as script)
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('load_ext', 'autoreload')
        ipython.run_line_magic('autoreload', '2')
except (ImportError, NameError):
    pass  # Not in IPython environment, skip magic commands

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
# Load all variants from WES data
print("=" * 80)
print("Loading all variants from all genes")
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
# Create patient X gene alteration matrix
# Columns: {gene}_mut (specific mutations), {gene}_amp (amplifications), 
#          {gene}_lof (loss-of-function), {gene} (summary: any alteration)
print("\n" + "=" * 80)
print("Creating patient X gene alteration matrix")
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
# Create oncoplot
create_oncoplot(patient_gene_matrix, figsize=(15, max(6, len(patient_gene_matrix.columns) * 0.2)))
