# fig2de_celltype_modules.ipynb

Origin: `pipelines/NicheNet_v3.ipynb` (Spencer Watson; originally `analysis/survival/NicheNet_v3.ipynb`),
cells 162 (Fig 2d), 172 (grid search adjusted for MGMT + age) and 174 (Fig 2e forest), numbered as
in the EFS copy `Final_Notebooks/NicheNet_v3.ipynb`. Cells 164–170 and 173 (exploratory) not copied.

Changes:
- Inputs with the regular loaders: clinical rows `AllDataCenters.load_clinical(exclude_idh1_mut=False)`
  (= the notebook's cell-3 `load_processed_clinical` concatenation, 313 rows, index reset as
  `ignore_index=True` did); snRNA-seq obs `AllDataCenters.load_singlecell(obs_only=True)` (no
  exclusions, as cell 5 loaded it; the matrices are not used); curated labels
  `annotation.map_levels(efs_copy=True)` (the EFS key cell 13 read).
- All functions (`is_excluded`, `build_proportion_matrix`, `clean_corr`, `compute_linkage`,
  `plot_corr_clustermap_with_clusters`, `build_cluster_dict`, `print_clusters`,
  `get_clusters_at_threshold`, `compute_score`, `stratify`, `run_rmst_mgmt_age`,
  `apply_stratification`, `run_rmst_module`, `run_modules`, `fmt_p`, `get_vals`, `y_frac`) are in
  one cell after the config; cell 172's identical redefinitions, the unused `contains_tcell` and
  the `if not in dir()` guards removed. Imports and the three cells' config constants at the top.
- Cell 174's per-module loop is `run_modules(tau, clin_df)` (body unchanged) so it can be rerun at
  τ = 3 years for the †.
- Outputs through `output.save` / `output.save_table`, `OUT = output.setup("fig2de_celltype_modules")`,
  panels Fig 2d, Fig 2d modules, Fig 2e, Fig 2e dagger.
- Added: check against the saved output of cell 174 (published Fig 2e), count of fitted patients
  whose missing MGMT / sex / p53 the `(x == level)` coding counts as the reference level, the †
  at τ = 3 years, and the corrected versions for the response letter (survival cohort with
  complete-case covariates; modules recomputed on the 266 samples of the paper cohort, loaded
  with `AllDataCenters.load_singlecell(exclude_idh1_mut=True, exclude_patient_ids=True,
  exclude_relapse=True, obs_only=True)`).
- Header: "Changes" list moved here; decision codes (D8) and source-cell references removed from
  section titles and comments.

Package helpers used:
- `annotation.map_levels(efs_copy=True)` / `annotation.NEW_COLS`: the annotation-key cell (13) of
  NicheNet_v3, shared with the other notebooks that read the EFS key.
- `cohort.clinical()`, `cohort.survival_cohort()`, `cohort.code_covariates()`,
  `cohort.n_coded_as_reference()`: the 273 / 260 cohort definitions and covariate coding used for
  the printed checks and the corrected version (not in the original notebook; no regular loader
  defines the paper cohorts).
- `paper.survival` (import only: R_HOME, survival + survRM2), `paper.output`.

Verification (re-run 2026-11-20, 11 min): all 11 CSV outputs byte-identical to the previous run
(`figures/paper/05_celltype_modules/`); printed output identical apart from loader log lines.
Fig 2e matches the saved cell-174 output for all 7 modules (e.g. module 9 +2.02 mo, p = 0.0008);
273 fitted patients; † at τ = 3 y: modules 9, 8, 4, 2; corrected (b): 266 samples, 23 modules,
15 / 25 identical.
