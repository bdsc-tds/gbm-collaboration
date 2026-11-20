# fig3b_0_cnv_wes_matrix.py

Origin: `analysis/infercnv/generate_cnv_wes_matrix.py` (WES CNV matrix, frozen as
results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv).

Changes:
- Header docstring (output, inputs, runtime, comparison with the frozen matrix); change list removed.
- The matrix is written to `output.setup("pipelines_generate_cnv_wes_matrix")`, not over the frozen
  file. `RESULTS_DIR` / `PATH_BINS_SINGLECELL` as in the original.
- Dropped Part 2 (exploration: PCA, outlier QC, single-sample plots).

Package helpers: `output` only.

Verification: not re-run (73 min of S3 reads, needs AWS credentials; frozen output). `py_compile` OK.
