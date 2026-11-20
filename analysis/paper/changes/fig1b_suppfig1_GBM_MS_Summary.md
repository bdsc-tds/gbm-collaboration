# fig1b_suppfig1_GBM_MS_Summary.ipynb

Origin: `analysis/survival/GBM_MS_Summary.ipynb` (Spencer Watson): Fig 1b, Supp Fig 1a-c and the
cohort numbers of the text.

Changes:
- Kept only the cells that produce paper output (original cells 4, 6, 44, 53, 54); removed the
  snRNA-seq and Visium loading, normalisation/UMAP, exports, ECM gene checks, TIL abundance,
  modality counts (Table 1a file), marker genes, exploratory clinical covariate cells and the
  exhausted T-cell cell.
- Clinical loading: `AllDataCenters.load_clinical(exclude_idh1_mut=False)` (the 313 rows of the
  original `load_processed_clinical` cell); the recurrence-kept count uses
  `AllDataCenters.load_clinical(exclude_idh1_mut=True, exclude_patient_ids=True, exclude_relapse=False)`.
  `paper_cohort` (273) and `surv` (260) from `cohort.clinical()` / `cohort.survival_cohort()`.
- New cell "Cohort numbers in the text"; new cell re-running the RMST model on the survival cohort.
- Original cell 44 used `clin_dedup`, defined only in cell 53 (run later): the deduplication lines
  are in the "Range in overall survival" cell.
- The PH-diagnostics cell (Supp Fig 1a-b) is a function `ph_diagnostics(clin_table, save_dir,
  suffix, panels)` in the functions cell of the notebook (run twice: N = 287 and the survival
  cohort). It was in `gbmcollaboration/paper/clinical_survival.py`; it is used only by this
  notebook, so it moved back and the module was deleted.
- Functions at the top: the RMST CONFIG block of cell 53 moved into the Parameters cell (the
  functions use TAU and the plot constants as defaults), and `ph_diagnostics`, `run_rmst`,
  `km_rmst_plot`, `km_center_plot`, `get_res`, `fmt_p`, `y_to_frac` in one functions cell after it.
  Cells 53 and 54 each defined `get_res` (same result for the default `key="adj"`): one definition
  (cell 53's) is kept.
- Cell 53 drew a forest plot that cell 54 redrew under the same file name: the drawing part of cell
  53 removed (`forest_df` and the summary table kept).
- R packages (survival, survRM2) loaded in the imports cell.
- Outputs: `OUT = output.setup("fig1b_suppfig1_GBM_MS_Summary")`, panels Fig 1b (age, mgmt, p53,
  sex, center), Supp Fig 1a, 1b, 1c.
- Removed review comments and decision codes (D10, D12, D13); header without the change list.

Package helpers: `cohort.clinical()` / `cohort.survival_cohort()` (the 273/260 cohort definition
shared by all survival files: ER_G_290a drop + SURVIVAL_EXCLUSION flag, not done by the loader);
`output`.

Verification: re-executed (1 min 30 s). All printed numbers identical to the previous run
(313 rows; 273 = UKER 135 / Charite 73 / CHUV 65; 260; 72/152/73; N = 287, Schoenfeld p MGMT
0.2031, age 0.0147, sex 0.8606; survival cohort 0.3839 / 0.0586 / 0.7660; RMST MGMT +7.1, age
-5.4, female +3.1, p53 +3.7; survival cohort p53 +2.9, p = 0.0847). All 11 PNGs byte-identical to
`figures/paper/02_GBM_MS_Summary/`. Only the order of the loader's thread-interleaved log lines differs.
