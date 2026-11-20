# suppfig3ef_analyses_wes_sc_infercnv.py

Origin: `analysis/infercnv/analyses_wes_sc_infercnv.py` (Supp Fig 3e-f, clone 4 test).

Changes:
- Only the obs of `infercnv_adata_malignant_only.h5ad` is read: `ad.read_h5ad(path, backed="r")`,
  obs copied, file closed (the source read the whole 16 GB object; was the package helper
  `wes.read_h5ad_obs`, now replaced by plain anndata backed mode).
- The WES label is added with a join on orig.ident, which keeps the cell index (the source's
  `obs.merge` replaced it with 0..n-1; same values). Prints the cells without a WES cluster.
- `frac_to_area` defined at the top.
- Figure saved with `output.save(..., panel="Supp Fig 3e-f")` instead of results/WES/infercnv/;
  `OUT = output.setup("suppfig3ef_analyses_wes_sc_infercnv")`.
- Dropped (not in the paper): the preview dot plot (redrawn as panel e), the unused
  `df_zoom_frac` table and `add_group_medians` helper, the clone 4 cell-state heatmap, and the
  cluster-2-vs-others cell-type composition test (its conclusion: "don't see significant
  differences").
- Header without the change list.

Package helpers: `output` only.

Verification: re-run (38 s). `WES_scCNV_integration_dot_facets.csv` byte-identical to
`figures/paper/08_analyses_wes_sc_infercnv/`; PNG identical except in panel f (the stripplot
jitter is unseeded). 782,542 cells; clone 4 test: n = 61 vs 205, U = 10592, p = 1.92e-16,
median difference 0.115.
