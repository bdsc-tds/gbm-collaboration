# fig3b_0_infercnv_sc.py

Origin: `analysis/infercnv/infercnv_sc_analyses.py` (single-cell inferCNV; frozen outputs in
results/scRNAseq/infercnv/).

Changes:
- Header docstring (outputs, readers, runtime, comparison with the frozen outputs); change list removed.
- `RESULTS_DIR = output.setup("pipelines_infercnv_sc_analyses")` instead of the frozen EFS folder
  (Part 2 reloads Part 1's object from it).
- Part 2 subset on `obs["tumor_lineage"]`, a column the file never creates (KeyError): it is
  set from `lineage_status == "tumor-lineage"` before the subset (on the frozen obs this gives
  exactly the 782,542 frozen cells).
- `display(dominant_df)` (notebook-only) -> `print`; unused imports removed.
- Loading: the original `AllDataCenters.load_singlecell(...)` (full object, needed for inferCNV).

Package helpers: `output` only.

Verification: not re-run (hours, well over 60 GB; outputs frozen). `py_compile` OK.
