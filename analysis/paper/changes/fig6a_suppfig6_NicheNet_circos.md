# fig6a_suppfig6_NicheNet_circos.ipynb

Origin: Spencer Watson, the circos cells of `fig6a_0_NicheNet_runs.ipynb` (`%%R` cells 80 TIL+,
116 BZ1, 152 BZ2; identical apart from the zone).

Changes (port and this pass):
- The R code of the cell kept as written, wrapped in one R function `circos_zone(AXES_SELECTED_CSV,
  zone_dir, OUT_PNG)` run through rpy2 (in place of the `%%R` cell magic), called once per zone.
  Only the head changed: input CSV, title zone label and output PNG are arguments in place of the
  `ZONE` → folder path lines; the function returns the plotted axes.
- This pass: the R code and the Python wrapper `circos_zone()` were in
  `gbmcollaboration/paper/nichenet.py`; used only here, so they are in the notebook's functions
  cell (with `save_panel`, `compare_with_published`) and the module is deleted. The R function is
  defined once when the cell runs (the module's `_defined` flag is gone).
- Figures go to `OUT = output.setup("fig6a_suppfig6_NicheNet_circos")` (the notebook wrote them next
  to `axes_selected.csv`), named by panel; plotted axes saved as `<panel>.csv`; PNG, PDF and CSV also
  in the panels folder (Fig6a, SuppFig6a, SuppFig6b).
- Added: comparison with the published PNGs, the TIL+ catalog check, the printout of the axes named
  in the Results.
- R circlize, dplyr and readr were added to the `gbm_main` environment (conda-forge).
- Header: change list moved here.

Package helpers used: `paper.config` (paths, R_HOME), `paper.output`.

Verification (re-run 2026-11-20, 31 s): printed output identical (46 / 54 / 44 edges, same ranges,
6 TIL+ axis IDs missing from the catalog, same Results axes); the three PNGs byte-identical to the
notebook's published PNGs.
