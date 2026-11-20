# fig4ab_suppfig4abg_Cell_Neighborhoods.ipynb

Origin: `analysis/survival/Cell_Neighborhoods.ipynb` (Spencer Watson).

Changes:
- Loading cells (2–12) replaced: `neighborhoods.load_spots()` (spots and coordinates with the
  loader exclusions; no cell uses the expression matrix, so it is not loaded),
  `AllDataCenters.load_cell2location(visium.sections())` merged with `visium.add_c2l_to_obs`,
  `annotation.add_celltype_key(efs_copy=True)` (the EFS key the notebook read). The clinical table
  (cell 4) was not used and is not loaded.
- Clustering: cell 36 (and the Hellinger representation it needs, first part of cell 14) run only
  with `RUN_CLUSTERING = True` and write to `OUT`; otherwise the frozen annotation and centroid files
  are read (`config.CT_ANNOTATIONS_CSV`).
- Removed the first strategy (Leiden clusters then kNN, superseded by cell 36): Leiden/UMAP part of
  cell 14, cells 16–25 (cell 22 also wrote `celltype_neighborhood_annotations_k10.csv`, with labels
  not sorted by size), 29 (earlier version of cell 38), 31 (earlier version of cell 46), 33
  (cell-type composition bars), 34 (exploratory); cell 1 (versions → `output.versions`).
- Cells 38 and 40 no longer reload the ECM file (cell 27 does).
- Cells 42 and 44: the cell2location → curated-type block is `neighborhoods.c2l_curated_counts`;
  the composition / specificity / marker-score matrices are `marker_scores` (same code, now defined
  in the notebook: `neighborhoods.marker_scores` removed from the package, used only here).
- Functions (`make_color_map`, `global_categories`, `color_array`, `joint_freq_pct`,
  `marker_scores`) in one cell after the parameters.
- Figures go to `OUT = output.setup("fig4ab_suppfig4abg_Cell_Neighborhoods")` with `output.save`
  (Fig 4a, 4b, Supp Fig 4a, 4b, 4g also to the panels folder); cell 46's matrix is compared with the
  one of the last run instead of overwriting it in `Analysis_Repo/Notebooks`.
- Fig 4a: the notebook's subset dotplot (cell 44: curated names, niches 0 and 6 left out) is not
  the published panel and no longer carries the panel tag. Added a reconstruction cell (`panel="Fig 4a"`):
  cell2location means with the uncurated scANVI names (names containing doublet / other / unknown
  removed), `dotplot_composite_score(abund_per_nh, included_niches=["0","1","2","3","4","5","7","8","9"],
  top_n_per_niche=5)`. It gives the published 36 cell types plus Mural_LEPR and Plasma_nos. A second
  version with `DROP = ["Mural_LEPR", "Plasma_nos"]` (`panel="Fig 4a published set"`, labelled as
  inferred) gives exactly the published 36, same order except near-ties (Neuron_neurogliaform with
  CT_0 instead of CT_9; Fibro_CAF-like / Mural_perivasc_fibro, Pericyte / Macro_LYVE1,
  Neuron_inh_PVALB / Neuron_exc_CDH24 swapped). Inference: the published panel was drawn from a
  cell2location run without those two factors (not on S3 / EFS).
- Per-cell imports moved to the top. Header: "Changes" list moved here.

Package helpers used:
- `gbmcollaboration.viz.cell2loc_plotting.dotplot_composite_score` (existing package plot function).
- `neighborhoods.load_spots()`: spot coordinates via the regular `AllDataCenters.load_spatial(load_gex=False)`;
  `load_spatial_adata(genes=...)` returns no coordinates.
- `neighborhoods.c2l_curated_counts`: shared with fig4cdg.
- `visium.sections()`, `visium.add_c2l_to_obs`, `annotation.add_celltype_key(efs_copy=True)`,
  `paper.config`, `paper.output`.

Verification (re-run 2026-11-20, 10 min): all 6 CSV outputs (Fig 4a, 4b, Supp Fig 4a, 4b, 4g)
byte-identical to the previous run (`figures/paper/12_Cell_Neighborhoods/`). Global NMI 0.140,
per-sample median 0.237; Supp Fig 4b matrix vs the notebook's last run: 228 samples, same order,
max |Δ| = 1.1e-16. Only printed difference: `T_other` is now before `total_cells` in `adata.obs`
(column order of `load_cell2location`), with no effect on any output.
Re-run with the Fig 4a reconstruction (10:00, 10 min): the 6 CSVs above unchanged; Fig 4a 38 types
(extra Mural_LEPR, Plasma_nos, none missing); Fig 4a published set: 36 types, same set as published.
