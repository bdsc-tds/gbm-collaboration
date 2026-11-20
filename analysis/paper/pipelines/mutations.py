"""Pipeline: WES patient x alteration matrix (frozen as results/WES/mutations/patient_gene_matrix.csv).

Output: boolean matrix, 299 samples x 34 columns (config.MUTATION_MATRIX_CSV), read by
06_oncoplot.py (Fig 3a) and the data release:
- <gene>_mut: listed SNV / indel (genes_mut.csv; TERT promoter C228T / C250T as TERTp_mut);
- <gene>_amp: CNVkit DUP call overlapping the gene with log2 ratio > 0.3 (genes_amp.csv);
- <gene>_lof: CNVkit DEL call with log2 ratio < -0.3, or a stop-gained / frameshift
  SNV/indel (genes_lof.csv); PTEN_lof is the "PTEN LOH" of the Results;
- EGFR, PTEN: any alteration of the gene (redundant, dropped by 06_oncoplot.py).
Logic in gbmcollaboration.utils.mutation_utils.create_comprehensive_patient_gene_matrix and
vcf_utils.read_vcf_all_genes.

Inputs: SNV/indel and CNV VCFs of the three centres on S3 (WES datasets DATA_ID_*["*_wes"]),
the three gene lists next to this file.

Runtime: tens of minutes (parallel VCF reads from S3), a few GB RAM. Not re-run for the paper.

Known issues: no sample exclusion here (relapse, IDH-mutant and excluded patients are in
the matrix; 06_oncoplot keeps the clinical cohort).

Changes: the gene lists are read next to this file (the source read them from the working
directory, so it only ran from its own folder); the matrix is written to the output folder
(`output.setup`), not over the frozen file; dropped the IPython autoreload cell and the
oncoplot-style heatmap (`create_oncoplot`, superseded by Fig 3a).
"""

# %%
# Import modules
import pandas as pd
import joblib
from tqdm import tqdm
from pathlib import Path

from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter
from gbmcollaboration.definitions import DATA_ID_CHUV, DATA_ID_UKER, DATA_ID_CR
from gbmcollaboration.utils.vcf_utils import read_vcf_all_genes
from gbmcollaboration.utils.mutation_utils import (
    parse_mutations,
    clean_gene_name,
    create_comprehensive_patient_gene_matrix,
)
from gbmcollaboration.utils.wes_io_utils import glob_s3_bucket
from gbmcollaboration.paper import output

# %%
RESULTS_DIR = output.setup("pipelines_mutations")
GENE_LIST_DIR = Path(__file__).parent  # genes_mut.csv, genes_amp.csv, genes_lof.csv

# %%
# Load driver gene lists
genes_mut_df = pd.read_csv(GENE_LIST_DIR / "genes_mut.csv")
genes_amp_df = pd.read_csv(GENE_LIST_DIR / "genes_amp.csv", header=None, names=['gene'])
genes_lof_df = pd.read_csv(GENE_LIST_DIR / "genes_lof.csv", header=None, names=['gene'])

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
# Create comprehensive patient X gene matrix
print("\n" + "=" * 80)
print("STEP 2: Creating comprehensive patient X gene alteration matrix")
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
# Save patient gene matrix to CSV
output_file = RESULTS_DIR / "patient_gene_matrix.csv"
patient_gene_matrix.to_csv(output_file)
print(f"\nSaved patient X gene matrix to {output_file}")
