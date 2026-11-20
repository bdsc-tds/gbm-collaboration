# fig3d_correlation_cnv_wes_ecm.py

Origin: `analysis/multimodal_correlation/correlation_cnv_wes_ecm.py`, Parts 1-2 (Fig 3d, PERMANOVA
numbers of the Results/Methods).

Changes:
- Single-cell data: `ad.AnnData(obs=AllDataCenters.load_singlecell(exclude_idh1_mut=True,
  exclude_patient_ids=True, exclude_relapse=True, obs_only=True))` (the original loaded the whole
  expression cohort; only obs is used). The original `adata.X = layers[...]` lines are gone with X.
- Annotation key read from the original EFS `SingleCell_Annotation_Key_v3.csv` as in the source
  (no package helper); sample IDs `str[:9]` as in the source.
- PERMANOVA: scikit-bio is not in the environment; `permanova()` (numpy/scipy Bray-Curtis
  PERMANOVA, same pseudo-F, seeded permutations) is defined at the top of this file (was
  `gbmcollaboration/paper/wes.py`, used only here). `stars` moved to the top too.
- New: level-3 PERMANOVA next to level 4 and the published values (`permanova_cnv_cluster.csv`);
  check of the Kruskal-Wallis table against the published `figures/figure3/fig3e.csv`.
- Outputs: `RESULTS_DIR = output.setup("fig3d_correlation_cnv_wes_ecm")`, Fig 3d panel
  (`fig3d_celltype_composition_cna_clusters`; the source saved it as figures/figure3/fig3e.pdf).
- Dropped (not in the paper): ECM Kruskal-Wallis and mean-proportion prints, the stacked barplot,
  the AC/OPC ratio test (KeyError on the level-4 table), unused parameters, Part 3 (pseudobulk
  DESeq2 and its dot plot).
- Header without the change list.

Package helpers: `output` only.

Verification: re-run (1 min 05 s). `fig3d_celltype_composition_cna_clusters.png/.csv`,
`kruskal_wallis_test_cellular_landscape_cluster.csv` and `permanova_cnv_cluster.csv`
byte-identical to `figures/paper/09_correlation_cnv_wes_ecm/`. PERMANOVA ECM F = 0.582 p = 0.868
(published 0.58 / 0.86), level 4 F = 4.02, level 3 F = 4.56 (published 4.56), p = 0.001; KW
table = published fig3e.csv (71/71 types, max |dH| 4e-15, same order).
