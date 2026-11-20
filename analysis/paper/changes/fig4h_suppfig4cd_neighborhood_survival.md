# fig4h_suppfig4cd_neighborhood_survival.ipynb

Origin: Spencer Watson, notebook Updated_figures (local, not in the repository): the RMST sweep under
"Updated ECM Neighborhood Plot" and the tau = 1.5 years check after it, as reproduced in
`analysis/scratch/review_2026-11/source_cells_spencer_2026-11.ipynb` cells 6–9.

Changes (port and this pass):
- The sweep cell is split into its Config block (Config cell), its functions (`clinical_frame`,
  `apply_split`, `run_rmst`, `fmt_p`, `run_sweep`, `forest_plot`, `dagger`, `describe_fit`; one
  functions cell after the Config) and the calls, so the same sweep runs on the ECM neighborhoods
  (Fig 4h) and on the cell-type neighborhoods (Supp Fig 4c–d). Analysis code as in the source cells.
- This pass: the functions were in `gbmcollaboration/paper/neighborhood_rmst.py`; they are used only
  by this notebook, so they are back at its top and the module is deleted.
- Clinical table: `AllDataCenters.load_clinical(exclude_idh1_mut=False)` in place of the local CSV
  export of the processed clinical tables (`clinical_frame`).
- Neighborhood labels from the frozen per-spot annotation CSVs (`config.ECM_ANNOTATIONS_CSV`,
  `config.CT_ANNOTATIONS_CSV`).
- Added: `describe_fit` (sections without a subject, fitted patients outside the 260-patient
  survival cohort, covariates coded as the reference level), the ECM_3 / ECM_1 print, and the
  CT_3 / CT_0 mutual adjustment table.
- Outputs: `OUT = output.setup("fig4h_suppfig4cd_neighborhood_survival")`, panels `Fig 4h`, `Supp Fig 4c`.
- Header: change narration and decision codes removed.

Package helpers used:
- `cohort.survival_cohort()` in `describe_fit`: the paper's 260-patient survival cohort definition
  (loader exclusions + the 13 Charité patients), shared by all survival files.
- `paper.config` (paths, R_HOME), `paper.output`.

Verification (re-run 2026-11-20, 56 s): all printed output identical to the previous run (ECM_3
adjusted −3.50 mo, p = 0.0114; ECM_1 0.2595; ECM_9 0.4743; † ECM_3 −1.7 mo p = 0.0129; CT_3 −3.91,
p = 0.0057; CT_0 +3.30, p = 0.0257; mutual adjustment table). fig4h / suppfig4c CSVs and PNGs and
the mutual-adjustment CSV byte-identical to `figures/paper/15_neighborhood_survival/`.
