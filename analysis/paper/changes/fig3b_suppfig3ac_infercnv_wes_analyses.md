# fig3b_suppfig3ac_infercnv_wes_analyses.py

Origin: `analysis/infercnv/infercnv_wes_analyses.py` (Fig 3b, Supp Fig 3a, Supp Fig 3c, CNV scores).

Changes:
- Sample IDs: only the obs of `infercnv_adata.h5ad`, read with `ad.read_h5ad(path, backed="r")`
  (12 s). The source read the whole object (1.35 M cells, 25 GB) and densified X into a
  per-sample pseudobulk used only to list the samples (hundreds of GB of RAM). The pseudobulk was
  sorted by sample (groupby), so the samples are taken in sorted order, as there.
- `RUN_CLUSTERING = False` (default): Leiden is not re-run, `cnv_leiden` is read from the frozen
  `results/WES/infercnv/leiden_cluster_wes_baseline_samples.csv` (`CNV_CLUSTERS_PATH`); PCA and the
  neighbour graph are still computed (UMAP). With True, the new table goes to the output folder.
- Outputs to `output.setup("fig3b_suppfig3ac_infercnv_wes_analyses")` and the panels folder
  (Fig 3b, Supp Fig 3a, Supp Fig 3c) instead of the scanpy figure folder and hard-coded EFS paths;
  the Tukey figure (saved by the source as figures/figure3/fig3c.pdf) is Supp Fig 3c; strict-loci
  and gene tables written to the output folder and compared with the frozen results/WES/ files.
- `sc = ax.scatter(...)`, which shadowed scanpy, is `sc_plot`; the Tukey rcParams applied in an
  `rc_context`.
- pyranges is not in the environment: `read_gtf_genes` / `genes_overlapping_segments` (pandas, same
  rule as `PyRanges.join`; identical gene table) are defined in this file's functions cell (were in
  `gbmcollaboration/paper/wes.py`, used only here; the module was deleted).
- Functions cell after the paths: `col_to_bin_id`, `chr_rank` (defined twice in the source, once
  here), `read_gtf_genes`, `genes_overlapping_segments`.
- Paths: original literals (`RESULTS_DIR`, matrix, bins, GTF) instead of `config.*`.
- New: printed checks against the published values, the ER_G_673a check cell.
- Header without the change list.

Package helpers: `output` only.

Verification: re-run (7 min 08 s). All outputs byte-identical to `figures/paper/07_infercnv_wes_analyses/`
(heatmap, UMAP, Supp Fig 3c PNG/CSV, strict loci, segment genes, label check, CNV scores).
CNV scores 0.170/0.172/0.197/0.183 (published 0.170/0.172/0.196/0.183; 0.196 with ER_G_673a in
cluster 2); same 3 segments (chr10, chr19, chr20) and bins as published, same gene lists
(1044/2940/1417 genes).
