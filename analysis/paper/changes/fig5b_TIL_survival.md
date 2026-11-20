# fig5b_TIL_survival.ipynb

Origin: Spencer Watson, TIL_Niche_Comparisons (EFS `Analysis_Repo/Notebooks/`), cell "RMST Survival
Analysis", as reproduced in `analysis/scratch/review_2026-11/source_cells_spencer_2026-11.ipynb` cells 16–17.

Changes (port and this pass):
- Spot table from the frozen TIL zone file (`til.til_zones()`) in place of `adata.obs` of the full
  Visium object (the cell only reads `sample` and `TIL_positive`).
- Clinical table: `AllDataCenters.load_clinical(exclude_idh1_mut=False)` (313 rows), the table the
  cell's `clin_out` was; this pass replaces the removed `cohort.clinical_as_ran()`.
- RMST through `rmst()` (survRM2 via rpy2, unadjusted + adjusted, results in months); this pass
  moves `rmst` and `to_r` from `gbmcollaboration/paper/survival.py` into the notebook (used only
  here) and removes `rmst` from the module.
- Figure settings in the Parameters cell; KM plot code as in the source cell; output folder
  `OUT = output.setup("fig5b_TIL_survival")`, panel `Fig 5b`.
- Added: printed checks of the join (sections without a clinical match, missing MGMT, sections
  without TIL+ spots); TIL+ fraction vs age and by MGMT; the sensitivity table (other covariate sets,
  without probe-batch blocks, corrected cohort).
- Bug fixed in an added check (this pass): "Fitted patients outside the paper cohort" compared
  9-character join keys with block IDs and the full cohort, so it always printed 0. It now lists
  the fitted patients outside the 260-patient survival cohort: 10 Charité patients (the same 10 as
  `describe_fit` in fig4h).
- Header: change narration and decision codes removed; comments neutral.

Package helpers used:
- `til.til_zones()`, `til.til_groups()`: the frozen TIL zone table and the TIL-high / TIL-low
  sections, shared with fig5c, fig5de and data_release.
- `cohort.code_covariates()`, `cohort.survival_cohort()`: covariate coding (as reference / complete
  cases) and the 260-patient survival cohort, shared by the survival files.
- `ids.parse_ids`, `paper.config` (R_HOME), `paper.output`.

Verification (re-run 2026-11-20, 58 s): all printed numbers identical to the previous run (n = 207,
54 / 153, +3.93 mo [0.5, 7.4], p = 0.0251; ρ = −0.22; sensitivity table) except the corrected
check line above. fig5b CSVs and PNG byte-identical to `figures/paper/17_TIL_survival/` (PDF differs
only in its creation date).
