# fig4f_pathologist_annotation_analyses.py

Origin: `analysis/hne/pathologist_annotation_analyses.py` (Owkin).

Changes:
- The annotation loading (DataConfig, AnalysisConfig, S3 reading, standardisation, exclusions,
  spot IDs, centre) comes from `gbmcollaboration.utils.pathology_annotation_loader`, of which the
  script held a copy; `load_and_prepare` merges its output with the ECM table as before.
- ECM table path from `gbmcollaboration.paper.config.ECM_ANNOTATIONS_CSV`.
- Output folder `OUT = output.setup("fig4f_pathologist_annotation_analyses")` instead of
  `results/hne/ecm_pathology_figures/`; Fig 4f also written to the panels folder (`panel="Fig 4f"`).
- Removed the figures not in the paper (agreed with the user): the three-panel centre
  replication (`fig_three_panel_centers`), the dot plot (`fig_dotplot`) and the sparse top-effects
  heatmap (`fig_sparse_top_effects`), with their calls; `savefig` import no longer needed.
- Added at the end: the Fig 4f table is read back and compared with the published CSV
  (`results/hne/ecm_pathology_figures/Fig_residual_heatmap.csv`).
- Per-file "Changes" list removed from the docstring (this file).

Package helpers used:
- `gbmcollaboration.utils.pathology_annotation_loader`: existing package module (committed,
  not paper-specific) holding the annotation loading the script duplicated.
- `paper.config` (paths), `paper.output` (output folder, panels).

Verification (re-run 2026-11-20, ~1 min): Chi2 = 77872.63 (dof 81), Cramér's V = 0.152;
batch-adjusted LR = 130.33 (df 81), p = 4.2e-4; Jaccard UKER vs CHUV = 0.781.
Fig 4f: same rows / columns / order as the published table, max |difference| = 0.
