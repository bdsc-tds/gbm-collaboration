# table1bc_cohort_description.py

Origin: `analysis/cohort_desc/cohort_description.py` (Table 1b-c).

Changes:
- Header docstring (panels, origin, comparison with the published table).
- Output: `OUT = output.setup("table1bc_cohort_description")`, table written with
  `output.save_table(..., panel="Table 1b-c")` instead of
  `descriptive_stats/processed/demographic_table_<DATE>_<PROJECT>.csv` (`PROJECT`, `DATE` removed).
- Loading is the original per-centre `DataCenter.load_tabular(exclude_idh1_mut, exclude_patient_ids,
  exclude_relapse)` + Baseline filter + ER_G_290a drop (`_load_center_df`, as written by the author);
  the excluded-features list stays in this file.

Package helpers: `paper.output` only (output folder and panel copy).

Verification: re-run (1 min 10 s). `demographic_table.csv` is byte-identical to the published
`descriptive_stats/processed/demographic_table_2026-06-10_Glimmer.csv` (centre n: UKER 135,
CHUV 65, CR 73).
