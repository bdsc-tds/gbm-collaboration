# fig4cdg_suppfig4f_ECM_niches_v4.ipynb

Origin: `analysis/survival/ECM_niches_v4.ipynb` (Spencer Watson).

Changes:
- Loading cells (2–12) replaced by the regular loaders:
  `AllDataCenters.load_spatial_adata(genes=ecm_genes_brain, list_sample_id=visium.sections(), normalization="log_norm_count")`
  (matrisome genes only; the full dense matrix is > 100 GB; cast to float64),
  `AllDataCenters.load_cell2location(visium.sections())` merged into `adata.obs` with
  `visium.add_c2l_to_obs` (the "Merge C2L data" cell), `annotation.add_celltype_key(efs_copy=True)`
  (the EFS key the notebook read), `AllDataCenters.load_clinical(exclude_idh1_mut=False)` for the
  clinical cell (all 313 rows, as `load_processed_clinical` concatenated them).
  Coordinates (`imagerow`, `imagecol`) come from `neighborhoods.load_spots()`; the notebook prints
  that the gene-subset load has the same 870,077 spots in the same order (the order of the frozen
  annotation file, which the adjacency permutations depend on) and indexes it by that order.
- Cell 14 split: the gene list comes before loading; the Supp Table 3a gene table (with the
  category of each gene from the comments of the list) is written there and compared with the
  published table (which also lists the 6 genes commented out of the list).
- Clustering cells 16, 18, 24, 40 run only with `RUN_CLUSTERING = True` and write to `OUT`;
  otherwise cell 42 reads `config.ECM_ANNOTATIONS_CSV`. Cell 20 (reload of the Leiden clusters) removed.
- Fig 4c: `sample_id` set to the published section CH_G_196a_vis (last run: ER_G_666a_vis).
- Cell 46 (full transcriptome): pseudobulk means from `neighborhoods.pseudobulk_sums` (per
  section × neighborhood sums, a few sections at a time) instead of the full matrix. The last run
  tested `imagerow` / `imagecol` as genes; they are left out (BH over 18,082 instead of 18,084
  tests; max change 1e-4, no signature gene changes).
- Cell 30: the cell2location → curated-type block is `neighborhoods.c2l_curated_counts`.
- Cell 58: `compute_contacts` counts edges with `count_contacts` (`np.add.at`) instead of a Python
  loop over 2.5 M edges × 1000 permutations; same label shuffles, same counts.
- Functions (`s3key`, `draw_heatmap`, `compute_contacts`, `count_contacts`) in one cell after the
  parameters. `count_contacts` (formerly `neighborhoods.contact_counts`, renamed: the notebook uses `contact_counts` as a variable) moved from `gbmcollaboration.paper.neighborhoods` into the notebook
  (used only here).
- Figures and tables go to `OUT = output.setup("fig4cdg_suppfig4f_ECM_niches_v4")`; paper panels
  also to the panels folder (Supp Table 3a, Fig 4c, Fig 4d, Supp Fig 4f, Fig 4g). Comparisons with
  the outputs of the last run (`Analysis_Repo/Notebooks`, `Final_Notebooks`) are printed.
- Per-cell imports moved to the top; versions cell → `output.versions`.
- Removed: 1 (versions), 2–12 (loading, see above), 20, 22 (UMAP), 26 (one-sample plot),
  28 (per-sample barplot; Supp Fig 4e is drawn in fig4e_suppfig4e_colocalization_plots), 34–38
  (survival, now fig4h_suppfig4cd_neighborhood_survival), 49 (composition-marker heatmap, not in the
  paper), 51–54 and 66 (TIL zones, now fig5de_suppfig5cde_ECM_supplement), 60 (summary table),
  62–63 (sample lookups).
- Header: "Changes" list moved here.

Package helpers used:
- `neighborhoods.load_spots()`: spot coordinates, read with the regular
  `AllDataCenters.load_spatial(load_gex=False)`; the gene-subset path of `load_spatial_adata`
  returns neither coordinates nor deconvolution.
- `neighborhoods.pseudobulk_sums`: full-transcriptome sums per section × neighborhood with each
  centre's regular `load_spatial`, a few sections at a time (all genes of all spots do not fit in
  memory; shared with fig4e for Supp Table 3b).
- `neighborhoods.c2l_curated_counts`: the cell2location → curated-type block (shared with fig4ab).
- `visium.sections()` (the 228 sections), `visium.add_c2l_to_obs` (merge cell),
  `annotation.add_celltype_key(efs_copy=True)` (annotation-key cell), `paper.config`, `paper.output`.

Verification (re-run 2026-11-20 under the heavy lock, 65 min; a first run failed at Fig 4g because
the moved function `contact_counts` clashed with the notebook variable of the same name, hence the
rename to `count_contacts`): all 11 CSV outputs byte-identical to the previous run
(`figures/paper/11_ECM_niches_v4/`). The gene-subset load has the same 870,077 spots in the same
order as `load_spots`. Supp Table 3a: 117 genes (published 123: also COL16A1, CTSD, IGFBP2, IGFBP7,
LGALS3, PLOD3, commented out of the list). Fig 4d: same 40 genes and order as the published panel.
Against the notebook's last run: ECM-restricted signatures 1170 / 1170 rows, full transcriptome
180,820 / 180,840 rows (the 20 missing are imagerow / imagecol), max |Δ logFC| 4e-16; Fig 4g
contacts, enrichment and permutation p values Δ = 0; OS Spearman max |Δ rho| 8e-17; 313 clinical rows.
