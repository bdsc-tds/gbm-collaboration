# fig2bc_suppfig2_marker_dotplots.ipynb

Origin: `analysis/annotation/marker_dotplots.ipynb` (Owkin, j-bac).

Changes:
- Loading as in the original: `AllDataCenters.load_singlecell(exclude_idh1_mut=True,
  exclude_patient_ids=True, exclude_relapse=True)`, counts layer into `.X`, `normalize_total` +
  `log1p` (the `snrna.load_singlecell` / `snrna.normalize_log1p` wrappers of the previous pass are
  gone).
- Level mapping loop: `snrna.add_curated_levels(ad, df_ct_map.set_index('celltype_level4_scanvi'))`
  (same loop as the notebook's, shared with fig2a).
- Paths from `paper.config` (marker workbook); `OUT = output.setup("fig2bc_suppfig2_marker_dotplots")`;
  figures through `output.save` with the panel id (`PANELS` dict: Fig 2b, 2c, Supp Fig 2c/f/i).
- Unused imports (seaborn, lifelines) and the commented-out key path removed; each figure is closed
  and the subset deleted after saving.
- `sc.pl.dotplot(..., use_raw=True)` made explicit and `assert ad.raw is not None` after the
  cell-type filter: the published tables (EFS `figures/scRNAseq/annotation/`) were computed from
  `adata.raw` (the integer counts of each centre's zarr `raw/X`), which `sc.pl.dotplot` uses by
  default (`use_raw=None`) when `.raw` is present. The `normalize_total` / `log1p` of
  `layers['unnormalized_counts']` in the load cell therefore does not reach the dot plots (kept as
  the authors wrote it). The earlier non-matching run came from the `snrna.load_singlecell` wrapper
  of the previous pass, which dropped `.raw`.
- Header cell: "Changes" list moved here.

Package helpers used:
- `snrna.add_curated_levels`: the level-1..4 mapping loop, used by fig2a and fig2bc.
- `paper.config` (paths), `paper.output`.

Verification (re-run 2026-11-20 under the heavy lock, 12 min): all 9 dot-plot tables equal the
published ones in EFS `figures/scRNAseq/annotation/` (same rows and order, mean expression and
fraction expressing max |difference| = 0), checked with scratchpad/dotplot/compare.py. The queued
confirmation run of the other agent (scratchpad/dotplot/full_run.py) was cancelled, as this run
covers it.
