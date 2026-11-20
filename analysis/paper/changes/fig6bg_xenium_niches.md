# fig6bg_xenium_niches.ipynb

Origin: `analysis/xenium/niches.ipynb` (j-bac), cells 0, 2, 4, 6, 28, 31 and 33 (moved with `git mv`).

Changes (port and this pass):
- Both panels run in one notebook: `panels` in the parameters cell replaces `panel`; cells 4, 6 and
  28 are wrapped in functions (`load_panel`, `prepare_panel`, `tile_panel`) run per panel, and the
  plotting cells loop over the panels.
- Removed cells 7–27, 29 and 34–48 (TIL pseudobulk DE, spatial TIL clusters, smoothed ULM, spatial
  regression, plots): exploratory, no paper panel; unused helpers of cell 0 (`run_ulm_smooth`,
  `spearman_cols`, `plot_correlation_heatmap`, `create_non_overlapping_grid`, `plot_grid`) and
  unused imports removed.
- This pass: `pseudobulk_lognorm`, `assign_grid_tiles`, `spearman_tils` were in
  `gbmcollaboration/paper/xenium.py`; used only here, so they are back in the notebook (verbatim from
  cell 0) and the module is deleted. All functions (these three, `load_panel`, `prepare_panel`,
  `tile_panel`, `compare`) are in one functions cell after the parameters; the emptied "Read and
  preprocess data" cells are gone.
- Fixed `from gbmcollaboration.viz import dotplot` (imports the module) →
  `from gbmcollaboration.viz.dotplot import dotplot`.
- Cell 4: disabled joblib cache and the `gex_only=False` branch removed; only the
  `split_fully_purified` counts are concatenated (raw samples still read for coordinates and labels).
- Cell 6: cell-level `normalize_total`/`log1p`, `spatial_neighbors`, `find_spatial_clusters`,
  `TIL_sum_cat` (used only by removed cells) removed.
- Paths: `config.EFS` / `config.XENIUM_BASE`; outputs `OUT = output.setup("fig6bg_xenium_niches")`,
  `OUT/<panel>/` with `output.save`, panels `Fig 6b hBrain`, `Fig 6b hAtlas`, `Fig 6g hBrain`,
  `Fig 6g hAtlas` (file stems kept).
- Added: comparison with the published CSVs (EFS `figures/xenium_fig_6.zip`); package versions.
- Header: change list moved here; edit narration removed from two comments.

Package helpers used: `paper.config` (paths), `paper.output`; `gbmcollaboration.utils.xenium_io_utils`
and `gbmcollaboration.viz.dotplot` as in the original.

Verification (re-run 2026-11-20; 205 min wall time, nearly all waiting for the heavy-run lock): same
cell, tile and category counts as before (hBrain 1,497,063 cells, 85,809 tiles; hAtlas 1,400,518 /
81,525); the check against the published CSVs is exact for both panels (Fig 6b 42 / 48 rows, Fig 6g
408 / 408, rho, p, p_adj identical). Fig 6b CSVs and all PNGs byte-identical to
`figures/paper/21_xenium_niches/`; the Fig 6g CSVs hold the same rows and values in a different row
order.
