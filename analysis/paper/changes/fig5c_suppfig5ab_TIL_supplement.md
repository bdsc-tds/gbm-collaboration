# fig5c_suppfig5ab_TIL_supplement.ipynb

Origin: Spencer Watson, `TIL_supplement.ipynb` (dropped in the repository root, November 2026; not
kept), all cells in his order.

Changes (port and this pass):
- Setup cell split into one imports cell and one Config cell. His helpers `load_til_zones`,
  `til_groups`, `zone_groups`, `VISIUM_EXCLUSIONS` / `visium_section_excluded` are in
  `gbmcollaboration.paper.til` (`til_zones()`, `til_groups()`, `zones_with_groups()`; same rules,
  exclusions through `gbmcollaboration.apply_exclusions`); `parse_ids` in `paper.ids`;
  `scanpy_style` in `paper.output`. `GROUPS` is re-indexed by `visium_id` with his `sample` and
  `TIL_percent` columns so his cells run as written.
- Fig 5c cell2location abundances: `AllDataCenters.load_cell2location(<228 sections of the zone
  file>)` (this pass; previously `visium.load_c2l`, removed) in place of his S3 loader
  (`load_sample`, `load_center`); `COL_PREFIX` is the literal column prefix.
- Annotation key: `annotation.celltype_key(efs_copy=True)` (the EFS copy he read).
- Functions at the top (this pass): `save` (from the Config cell) and `scale_size` (from the
  Fig 5c plotting block) are in one functions cell after the Config.
- Outputs: `OUT = output.setup("fig5c_suppfig5ab_TIL_supplement")`, panels Supp Fig 5a, 5b, Fig 5c;
  `SAVE_FIGURES` flag dropped; versions cell → `output.versions`.
- Added: the check of Fig 5c against `Analysis_Repo/Notebooks/TIL_zone_celltype_binary_bubble.csv`.
- Header: change list moved here; references to `fig5b_TIL_survival.ipynb` updated.

Package helpers used:
- `til.til_groups()`, `til.zones_with_groups()`: TIL-high / TIL-low assignment and the zone table,
  shared with fig5b, fig5de and data_release.
- `annotation.celltype_key(efs_copy=True)`: annotation key (EFS copy), shared with other files.
- `ids.parse_ids`, `output.scanpy_style`, `paper.config`, `paper.output`.

Verification (re-run 2026-11-20, 1 min 22 s): all printed output identical to the previous run
(cut-off 0.030350, 57 / 171, 57 / 148 in zones; ρ 0.600, OR 14.57; Fig 5c 276 tests, same 20 cell
types; check vs published table unchanged). All CSVs and PNGs byte-identical to
`figures/paper/16_TIL_supplement/` (PDFs differ only in creation date).
