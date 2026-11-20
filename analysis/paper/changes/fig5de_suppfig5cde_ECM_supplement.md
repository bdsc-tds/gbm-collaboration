# fig5de_suppfig5cde_ECM_supplement.ipynb

Origin: Spencer Watson, `ECM_supplement.ipynb` (dropped in the repository root, November 2026; not
kept), all cells in his order.

Changes (port and this pass):
- Setup cell split into one imports cell and one Config cell. His helpers `load_til_zones`,
  `til_groups`, `zone_groups`, `VISIUM_EXCLUSIONS` / `visium_section_excluded` are in
  `gbmcollaboration.paper.til` (`til_groups()`: `TIL_group_zones` is his `zone_groups`;
  `zones_with_groups()`); `scanpy_style` in `paper.output`; unused `parse_ids`, `ids_to_set`,
  `first_block_per_patient` dropped. Same groups (57 / 171; 57 / 148 in zone comparisons).
- Expression: `AllDataCenters.load_spatial_adata(genes=table_genes, list_sample_id=visium.sections(),
  normalization="log_norm_count")` (this pass; previously the removed `visium.load_expression`) in
  place of the dense spots × all-genes frame of the whole cohort (> 100 GB). Same loader exclusions,
  layer and gene space; only the panel genes and the gene-table genes are kept, so `VAR_NAMES` is the
  gene space restricted to them.
- `SUMMARY_LIST_119` and `CHECK_GENES` moved from the gene-table cell to the Config cell (the
  expression cell needs them); the Supp Table 4a gene list is read there too.
- Functions at the top (this pass): `save`, `make_bubble`, `star`, `paired_panel`, `per_sample_means`,
  `test_per_zone_category`, `box_pair`, `stars` and `code_cells` moved from their cells into one
  functions cell after the Config. `make_bubble(..., top_n=TOP_N)` defaulted to the cell's `TOP_N`;
  `top_n` is now passed at the call (same value, 15).
- Fig 5d: `MIN_ZONE_SPOTS = 0` as he ran it; added cell printing the same tests on section-zones with
  ≥ 20 spots (the filter of `Final_Notebooks/focus_niches_significance.csv`).
- Added the Supp Table 4a gene-use table (for each gene, the analysis files whose code uses it).
  This pass: the notebook's own name in `SKIP` updated, and `analysis/scratch/` (unused code, now
  under `analysis/`) left out of the search.
- Outputs: `OUT = output.setup("fig5de_suppfig5cde_ECM_supplement")`, panels Supp Fig 5c, Fig 5d,
  Fig 5e, Supp Fig 5d, 5e, Supp Table 4a; `SAVE_FIGURES` dropped; versions → `output.versions`.
- Added after each panel: comparison with the table of the published run.
- Header: change list moved here.

Package helpers used:
- `til.til_groups()`, `til.zones_with_groups()`: TIL group assignment and zone table, shared with
  fig5b, fig5c and data_release.
- `visium.sections()`: the 228 Visium sections of the paper (list_sample_id of the gene-subset loader).
- `output.scanpy_style`, `paper.config`, `paper.output`.

Verification (re-run 2026-11-20; 205 min wall time, almost all of it waiting for the heavy-run lock):
all printed numbers identical to the previous run (Supp Fig 5c, Fig 5d both filters, Fig 5e n = 192,
177 / 192, p = 1.66e-30, Supp Fig 5d–e, matrisome gene table, comparisons with the published tables).
All panel CSVs and PNGs byte-identical to `figures/paper/18_ECM_supplement/` (PDFs differ only in
creation date). Only `supp_table_4a_gene_use.csv` changed: 388 code units searched (319 before), since
new paper files exist; the collagen genes and IGFBP2 now also list `supptables_4a_ecm_til_zone_tables.ipynb`
(the cell-72 matrisome list of TIL_Niche_Comparisons_v2_clean). CTSD, IGFBP7, LGALS3 are still used
by no analysis.
