# fig4e_suppfig4e_colocalization_plots.ipynb

Origin: `analysis/cell2loc/colocalization_plots.ipynb` (Owkin, j-bac).

Changes:
- Paths from `gbmcollaboration.paper.config`; figures and tables go to
  `OUT = output.setup("fig4e_suppfig4e_colocalization_plots")` with `output.save` /
  `output.save_table` (Fig 4e, Supp Fig 4e, Supp Tables 3b–c also to the panels folder) instead of
  EFS `figures/visium/` and `supplementary_tables/`; the published versions are read back and compared.
- Loading (cell 7): `neighborhoods.load_spots()` gives the q05 abundances that
  `load_spatial_adata(load_deconvolution=True)` put in `adata.obs`, without the expression matrix;
  cell 9's gene means per neighborhood come from `neighborhoods.pseudobulk_sums` (same means).
- Fig 4e (cell 16): run with the defaults of `dotplot_composite_score` (the published version;
  the cell last ran with `center_zero=True, cmap='bwr'`).
- Bug fixed: `col_prefix` kept the "ECM_" of the previous annotation for `pathology_annotations`
  (cells 12–24); now empty for every annotation but the ECM neighborhoods.
- Removed: cell 2 (FileLink), 8 (duplicate of 7), 12 (clustermap dotplots), 14 (earlier version of
  16), 18 (`plot_neighborhood_marker_dotplot_from_agg` not defined: NameError), 20 and 22
  (pathology × neighborhood dotplots, not paper panels; Fig 4f is the residual heatmap of
  fig4f_pathologist_annotation_analyses.py), unused imports.
- Header: "Changes" list moved here.

Package helpers used:
- `neighborhoods.load_spots()`: the q05 deconvolution table of the full `load_spatial_adata`, read
  with the regular `AllDataCenters.load_spatial(load_gex=False)` (the full object with expression
  is > 100 GB; `load_spatial_adata(genes=...)` does not return the deconvolution).
- `neighborhoods.pseudobulk_sums`: gene sums per section × ECM neighborhood with each centre's
  regular `load_spatial`, a few sections at a time (Supp Table 3b needs every gene; shared with fig4cdg).
- `gbmcollaboration.utils.pathology_annotation_loader` (existing package module, used by the original).
- `paper.config`, `paper.output`.

Verification (re-run 2026-11-20 under the heavy lock, 13 min): all 28 CSV outputs byte-identical
to the previous run (`figures/paper/13_colocalization_plots/`). Against the published files:
Supp Table 3b (10 × 18,082, same genes) max |difference| 1.8e-15; Supp Table 3c (12 tables)
≤ 1e-16; Fig 4e plotted abundances and Supp Fig 4e bar data (228 × 10) max |difference| 0.
