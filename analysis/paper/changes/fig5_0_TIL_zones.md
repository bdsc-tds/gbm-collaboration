# fig5_0_TIL_zones.ipynb

Origin: Spencer Watson, EFS `Analysis_Repo/Notebooks/TIL_Niche_Comparisons_v2_clean.ipynb`, cells 21
(TIL+ definition; cell 19 is an identical copy), 22 (zones) and 24 (write the table), with their
markdown.

Changes (port and this pass):
- Spots and image coordinates from `AllDataCenters.load_spatial(load_gex=False)` (same spots and
  coordinates as `load_spatial_adata`, without the > 100 GB dense expression these cells do not
  use); `adata` holds only `obs` and the annotation key.
- Cell2location abundances: `AllDataCenters.load_cell2location(<sections>)` (this pass; previously
  the removed `visium.load_c2l`), joined onto obs with `visium.add_c2l_to_obs`.
- Annotation key: `annotation.add_celltype_key(adata, efs_copy=True)` (the EFS copy the notebook read).
- This pass: `get_coords` (defined inside the zone cell) moved to a functions cell after the
  Parameters.
- Cell 24 writes to `OUT = output.setup("fig5_0_TIL_zones")` unless `RUN_PIPELINE = True` (then to the
  frozen `config.TIL_ZONES_CSV`); added the comparison with the frozen file.
- Header: change list moved here.

Package helpers used:
- `visium.add_c2l_to_obs`: the (sample, barcode) join of the abundances onto obs, shared with
  fig4ab and fig4cdg.
- `annotation.add_celltype_key(efs_copy=True)`: annotation key into `uns`, shared with other files.
- `paper.config` (loader flags, paths), `paper.output`.

Verification (re-run 2026-11-20, 8 min 40 s): 870,077 spots, 228 sections, zones TIL− 726,937 /
BZ1 53,474 / TIL+ 45,504 / BZ2 44,162; identical to the frozen file (same spot set, 0 spots differ
in sample, TIL_positive, TIL_zone). Output identical to the previous run apart from the output path.
Supp Table 4a (supp-tables pass, 2026-11-20): the zone table is also saved with
`output.save_table(OUT, "til_zone_annotations", df_til, panel="Supp Table 4a til_zone_annotations", index=False)`
(panels folder `SuppTable4a_til_zone_annotations.csv`); the last cell checks it is byte-identical to the published
`supplementary_tables/5a_til_zone_annotations.csv`: True (re-run 8 min 25 s, zones unchanged, 0 spots differ).
