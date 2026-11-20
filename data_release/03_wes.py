"""WES exports (somatic, summarised only; no reads/VCFs/germline/HLA).

1. cnv_log2ratio_sample_x_bin.csv : mean CNVkit log2 copy ratio per sample per genomic
   bin (autosomes; the 200-gene / step-10 windows shared with the scRNA inferCNV
   analysis), as used for Fig 3b clustering.
2. cnv_bins.csv                   : bin_id, chromosome, start, end (GENCODE v32, GRCh38).
3. cnv_clusters.csv               : WES CNV Leiden cluster (0-3) per sample (Fig 3b-d).
4. oncogenic_alterations_sample_x_gene.csv : 0/1 per sample for oncogenic SNV/indel
   (_mut, incl. TERT promoter C228T/C250T), amplification (_amp) and loss of function
   (_lof) as plotted in the Fig 3a oncoplot.

All restricted to the 273-patient cohort.
"""

import pandas as pd

import config as C
from common import block_of, out


def main():
    clin = pd.read_csv(out("clinical", "clinical_binned.csv"), index_col="sample_id")
    cohort = clin.index

    # --- CNV bins
    bins = pd.read_csv(C.CNV_BINS_CSV, index_col=0).rename(columns={"bin": "bin_id"})
    cnv = pd.read_csv(C.CNV_MATRIX_CSV, index_col=0)
    cnv.index = block_of(cnv.index)
    cnv.columns = cnv.columns.astype(int)
    bins = bins.set_index("bin_id").loc[cnv.columns]
    # The saved matrix still holds the 49 chrX bins (and lacks chr21's single bin 1274
    # and chrY); infercnv_wes_analyses.py removes chrX/Y before clustering, so the
    # paper's Fig 3b is autosomal. Release the same autosomal set.
    auto = ~bins.chromosome.isin(["chrX", "chrY"])
    print(f"Dropping {int((~auto).sum())} sex-chromosome bins; chr21 present: {'chr21' in set(bins.chromosome)}")
    bins, cnv = bins[auto], cnv.loc[:, auto.to_numpy()]
    n_all = len(cnv)
    cnv = cnv.loc[cnv.index.intersection(cohort)].sort_index()
    cnv.index.name = "sample_id"
    cnv.to_csv(out("wes", "cnv_log2ratio_sample_x_bin.csv"), float_format="%.5g")
    bins.reset_index().to_csv(out("wes", "cnv_bins.csv"), index=False)
    print(f"CNV: {n_all} samples in matrix -> {len(cnv)} in cohort; {cnv.shape[1]} bins")

    # --- CNV clusters
    cl = pd.read_csv(C.CNV_CLUSTERS_CSV, index_col=0)["cnv_leiden"]
    cl.index = block_of(cl.index)
    cl = cl.loc[cl.index.intersection(cohort)].rename("cnv_cluster")
    cl.index.name = "sample_id"
    cl.sort_index().to_csv(out("wes", "cnv_clusters.csv"))
    print("CNV clusters:", len(cl), cl.value_counts().sort_index().to_dict())

    # --- Oncogenic alterations (Fig 3a)
    mut = pd.read_csv(C.MUTATION_MATRIX_CSV, index_col=0)
    mut = mut.drop(columns=["EGFR", "PTEN"])  # redundant aggregate columns, dropped for the oncoplot too
    mut.index = block_of(mut.index)
    n_all = len(mut)
    mut = mut.loc[mut.index.intersection(cohort)].astype(int).sort_index()
    mut.index.name = "sample_id"
    mut.to_csv(out("wes", "oncogenic_alterations_sample_x_gene.csv"))
    print(f"Mutations: {n_all} samples in matrix -> {len(mut)} in cohort; {mut.shape[1]} alteration columns")
    print((mut.mean() * 100).round(1).sort_values(ascending=False).head(8).to_dict())


if __name__ == "__main__":
    main()
