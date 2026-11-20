# fig5fg_TIL_zone_GSEA.ipynb

Origin: Spencer Watson, EFS `Analysis_Repo/Notebooks/TIL_Niche_Comparisons_v2_clean.ipynb`, cell 63
("GSEA lollipop plots v4") with its markdown (cell 62).

Changes (port and this pass):
- `GSEA_OUTDIR` (set by cell 61 in the source) is the frozen folder, read only; figures go to
  `OUT = output.setup("fig5fg_TIL_zone_GSEA")` (same `fig.savefig` call, 200 dpi) and to the panels
  folder (Fig 5f, Fig 5g TIL+ vs BZ1, Fig 5g BZ1 vs BZ2), plotted rows as CSV.
- Imports moved to the imports cell; `savefig.transparent` reset to False after `output.setup`
  (the source saved on white).
- This pass: the cell's CONFIG block and its constants (`COMP_LABELS`, `ZONE_COLORS`, `BLOCK_GAP`)
  moved to the Parameters cell, its helpers (`get_library`, `is_excluded`, `_has_force_pattern`,
  `select_top`, `draw_lollipop`) to a functions cell after it; the rest of the cell as written.
- Added: frozen table compared with Supp Table 5b; PNGs compared with those of the source run.
- Header: change list moved here, stale reference to `paper/12_til_zone_gsea.py` removed.

Package helpers used: `paper.config` (paths), `paper.output`.

Verification (re-run 2026-11-20, 23 s): printed output identical; the three PNGs byte-identical to
the source run's PNGs in `GSEA_OUTDIR`; frozen table identical to Supp Table 5b.
