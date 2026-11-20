# fig3c_suppfig3b_wes_profile_os_rmst.ipynb

Origin: `analysis/survival/wes_profile_os_rmst.ipynb` (Spencer Watson): Fig 3c, Supp Fig 3b.

Changes:
- Clinical rows from the notebook's own `load_clinical_data()` (per-centre `load_tabular` with
  the loader exclusions, SURVIVAL_EXCLUSION dropped on CR; 261 rows), now called (the call was
  commented out in the original, which relied on an earlier session state).
- `COVARIATE_COLS`: "p53" replaces "Site code" (Spencer's revision of the model); p53 coded
  Negative 0 / Positive 1 in the same `mapping` dict. `load_and_merge_cluster_data` and
  `run_rmst_and_km_cluster_vs_rest` take the covariate list as an argument (default
  `COVARIATE_COLS`) so the published model runs with the same code.
- `plot_global_km_curves` gave `kmplot` OS in years, unsigned, with `convert_days_to_months=True`:
  kmplot cast the durations to whole years and read the sign as the event. It now passes the
  signed OS in days. The global log-rank test is unchanged.
- KM confidence bands are rasterized in the PDF (`rasterize_ci_bands`), so they keep their
  transparency (they rendered as solid blocks in the manuscript PDF).
- Removed the per-cluster Cox models and the categorical Cox forest (not in the paper).
- Outputs: `OUT = output.setup("fig3c_suppfig3b_wes_profile_os_rmst")`; the data saved with the KM
  figures are the KM curves; forest plot saved with `output.save` (panel Supp Fig 3b); Fig 3c
  panel tag on cluster 2 vs rest.
- New cells: the published model, a table of both models next to the published numbers, versions.
  `row` (model-table helper) is in the functions cell.
- Removed decision codes (D9) and review comments; header without the change list.

Package helpers: `ids.patient_of` (printed patient count of the joined rows); `output`.

Verification: re-executed (1 min). All printed numbers identical to the previous run: 261
clinical rows, 253 joined (252 patients); model age+MGMT+sex+p53 n = 173, cluster 2 vs rest
+5.08 mo [1.52, 8.64], p = 0.0052; published model n = 247, +3.35 mo [0.28, 6.43], p = 0.0326
(= published). All KM / forest PNGs and summary CSVs byte-identical to
`figures/paper/10_wes_profile_os_rmst/` (fig3c_cluster2_models.csv differs only in the model label).
