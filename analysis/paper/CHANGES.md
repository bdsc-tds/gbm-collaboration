# Changes to the paper analyses

For each file: where it comes from, what was changed and why, which shared helpers it uses, and how its outputs were checked. Files are in the order of the README table; `git log --follow -p <file>` shows the full diff.

## table1a_sample_number_per_modality.ipynb

Origin: `analysis/cohort_desc/sample_number_per_modality.ipynb` (Supp Table 1a).

Changes:
- `load_wes` loaded the clinical table instead of the WES table, so the published "WES" and
  "WES_surv" rows (274 / 261) counted clinical rows. The WES table is loaded
  (`load_clinical_and_wes`, user's fix, plus the per-centre summary cell).
- Visium backup sections (`_vbu`, `_vb2`, `_bu2`, `_bu3`) were cut to `CH_G_185a_vb` etc. and matched
  no block; sections are mapped to their block with `ids.block_of`, blocks to patients with
  `ids.patient_of` (was `.str[:-1]`).
- Modality-combination counts counted per patient; the as-published values are recomputed next to
  them with the two bugs (column "as published").
- snRNA sample names from `AllDataCenters.load_singlecell(exclude_idh1_mut=True,
  exclude_patient_ids=True, exclude_relapse=True, obs_only=True)` (was the removed
  `snrna.load_singlecell_obs`; the original loaded the full object).
- Functions (`load_clinical_and_wes`, `get_center`, `modality_counts`) in the first code cell.
- Tables written with `output.save_table` (`OUT = output.setup("table1a_sample_number_per_modality")`,
  panel "Table 1a").
- Header without the change list.

Package helpers: `ids.block_of` / `ids.patient_of` (parse section -> block -> patient IDs incl.
the Visium backup suffixes, which the original `.str` slicing got wrong; shared with other
files); `output`.

Verification: re-executed (about 10 min, mostly `AllDataCenters.load_spatial`). Outputs identical
to the previous run: snRNA samples 266; per-centre table (Total 273 / 250 / 228 / 273 / 265 / 263
/ 204); Table 1a "as published" 228/266/274/261/201/207/197/266/254/201/192 (= published Supp
Table 1a), "corrected" 228/265/263/250/222/218/208/256/244/213/204.

## table1bc_cohort_description.py

Origin: `analysis/cohort_desc/cohort_description.py` (Table 1b-c).

Changes:
- Header docstring (panels, origin, comparison with the published table).
- Output: `OUT = output.setup("table1bc_cohort_description")`, table written with
  `output.save_table(..., panel="Table 1b-c")` instead of
  `descriptive_stats/processed/demographic_table_<DATE>_<PROJECT>.csv` (`PROJECT`, `DATE` removed).
- Loading is the original per-centre `DataCenter.load_tabular(exclude_idh1_mut, exclude_patient_ids,
  exclude_relapse)` + Baseline filter + ER_G_290a drop (`_load_center_df`, as written by the author);
  the excluded-features list stays in this file.

Package helpers: `paper.output` only (output folder and panel copy).

Verification: re-run (1 min 10 s). `demographic_table.csv` is byte-identical to the published
`descriptive_stats/processed/demographic_table_2026-06-10_Glimmer.csv` (centre n: UKER 135,
CHUV 65, CR 73).

## supptables_1bc_pathology_clinical_variables.ipynb

Origin: 1b — annotation loading of analysis/cell2loc/colocalization_plots.ipynb (now fig4e_suppfig4e_colocalization_plots.ipynb): `pathology_annotation_loader.load_and_prepare(keep_only_shared_annotations=False)`, distinct labels in order of appearance. 1c — `analysis/scratch/survival/survival_patient_features.ipynb`, `load_clinical()` + the commented `df_clin_processed.columns.to_frame().to_csv('1c_clinical_variables.csv', index=False)` (published file without header line).
Data: `AllDataCenters.load_clinical(exclude_idh1_mut=True)` (centre column moved first, as the original concat produced) + the four os/pfs day columns of the original function; S3 pathology annotations through the existing package loader.
Outputs: SuppTable1b.csv, SuppTable1c.csv (panels folder).
Verification: 1b identical (11 labels, same order; published file lacks the final newline); 1c byte-identical (255 variables). Runtime 45 s.

## fig1b_suppfig1_GBM_MS_Summary.ipynb

Origin: `analysis/survival/GBM_MS_Summary.ipynb` (Spencer Watson): Fig 1b, Supp Fig 1a-c and the
cohort numbers of the text.

Changes:
- Kept only the cells that produce paper output (original cells 4, 6, 44, 53, 54); removed the
  snRNA-seq and Visium loading, normalisation/UMAP, exports, ECM gene checks, TIL abundance,
  modality counts (Table 1a file), marker genes, exploratory clinical covariate cells and the
  exhausted T-cell cell.
- Clinical loading: `AllDataCenters.load_clinical(exclude_idh1_mut=False)` (the 313 rows of the
  original `load_processed_clinical` cell); the recurrence-kept count uses
  `AllDataCenters.load_clinical(exclude_idh1_mut=True, exclude_patient_ids=True, exclude_relapse=False)`.
  `paper_cohort` (273) and `surv` (260) from `cohort.clinical()` / `cohort.survival_cohort()`.
- New cell "Cohort numbers in the text"; new cell re-running the RMST model on the survival cohort.
- Original cell 44 used `clin_dedup`, defined only in cell 53 (run later): the deduplication lines
  are in the "Range in overall survival" cell.
- The PH-diagnostics cell (Supp Fig 1a-b) is a function `ph_diagnostics(clin_table, save_dir,
  suffix, panels)` in the functions cell of the notebook (run twice: N = 287 and the survival
  cohort). It was in `gbmcollaboration/paper/clinical_survival.py`; it is used only by this
  notebook, so it moved back and the module was deleted.
- Functions at the top: the RMST CONFIG block of cell 53 moved into the Parameters cell (the
  functions use TAU and the plot constants as defaults), and `ph_diagnostics`, `run_rmst`,
  `km_rmst_plot`, `km_center_plot`, `get_res`, `fmt_p`, `y_to_frac` in one functions cell after it.
  Cells 53 and 54 each defined `get_res` (same result for the default `key="adj"`): one definition
  (cell 53's) is kept.
- Cell 53 drew a forest plot that cell 54 redrew under the same file name: the drawing part of cell
  53 removed (`forest_df` and the summary table kept).
- R packages (survival, survRM2) loaded in the imports cell.
- Outputs: `OUT = output.setup("fig1b_suppfig1_GBM_MS_Summary")`, panels Fig 1b (age, mgmt, p53,
  sex, center), Supp Fig 1a, 1b, 1c.
- Removed review comments and decision codes (D10, D12, D13); header without the change list.

Package helpers: `cohort.clinical()` / `cohort.survival_cohort()` (the 273/260 cohort definition
shared by all survival files: ER_G_290a drop + SURVIVAL_EXCLUSION flag, not done by the loader);
`output`.

Verification: re-executed (1 min 30 s). All printed numbers identical to the previous run
(313 rows; 273 = UKER 135 / Charite 73 / CHUV 65; 260; 72/152/73; N = 287, Schoenfeld p MGMT
0.2031, age 0.0147, sex 0.8606; survival cohort 0.3839 / 0.0586 / 0.7660; RMST MGMT +7.1, age
-5.4, female +3.1, p53 +3.7; survival cohort p53 +2.9, p = 0.0847). All 11 PNGs byte-identical to
`figures/paper/02_GBM_MS_Summary/`. Only the order of the loader's thread-interleaved log lines differs.

## fig2a_suppfig2_figure2.py

Origin: `figures/figure2.py` (Owkin, Valerie Ducret).

Changes:
- Loading as in the original: `AllDataCenters.load_singlecell(exclude_idh1_mut=True,
  exclude_patient_ids=True, exclude_relapse=True)` (full object, `.X` = log1p counts as loaded).
  The `snrna.load_singlecell` wrapper and the `UMAP_INPUT` switch of the previous pass are removed.
- `compute_pca_umap` and `downsample_stratified` are back in the script (verbatim from
  figure2.py; used only here), defined at the top with the other functions. The `uns`
  pca / neighbors / umap reset loop of the original is restored.
- Annotation key: the repo key (`annotation.celltype_key()`, `SingleCell_Annotation_Key_v3b.csv`)
  instead of EFS `resources/SingleCell_Annotation_Key_v3.csv`. Level-1 and level-4 labels are
  identical; at level 2 Mono_Macro is Mono_Mac, at level 3 DC1 / DC2 are cDC1 / cDC2; v3b maps
  Doublet_TNK_nos to level-1 Doublet (v3 left it unmapped, hence the extra name in the doublet
  filter, kept).
- Bug: figure2.py reassigned `adata_filt = compute_pca_umap(adata_filt)`, so the lineage subsets
  were taken from the scaled 2,000-HVG matrix and normalised + logged again (NaN; the script fails
  at highly_variable_genes). Fig 2a is computed on a copy; the lineage subsets come from
  `adata_filt` as loaded.
- Bug: the filter `"malignant_other_1"` (lower case) matched nothing. Malignant_other_1 is level-1
  Other and already removed by the Other filter; the corrected filter is kept and prints that it
  removes 0 nuclei.
- Fig 2a saves level 4 only (the panel) plus the UMAP coordinates; the three copy-pasted
  stacked-bar blocks are one function (`stacked_bar`), the lineage blocks one function
  (`lineage_panels`); bars are drawn before each lineage UMAP (they only use obs).
- Outputs: `OUT = output.setup("fig2a_suppfig2_figure2")`, `output.save` (pdf + png + plotted data
  csv) with `panel=` tags (Fig 2a, Supp Fig 2a/b/d/e/g/h per level); unused imports removed.
- Docstring: "Changes" list moved here.

Package helpers used:
- `snrna.add_curated_levels`: the level-1..4 mapping loop, shared with fig2bc.
- `annotation.celltype_key()`: the repo annotation key (indexed by level-4 scANVI label).
- `paper.output`.

Verification (re-run 2026-11-20 under the heavy lock, 120 min, ~150–205 GB RSS): all 13 CSVs
(1,350,338 nuclei, 266 samples loaded; 1,123,420 annotated nuclei, 265 samples, 70 level-4 types;
Fig 2a UMAP coordinates, Supp Fig 2b/e/h stacked-bar tables, levels 1–4) byte-identical to the
previous run (`figures/paper/03_figure2/`); same 63 files. Fig 2a done after 44 min; tumour
nuclei 747,089 → 375,148 after downsampling; Malignant_other_1 filter removes 0 nuclei.

## fig2bc_suppfig2_marker_dotplots.ipynb

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

## fig2de_celltype_modules.ipynb

Origin: `fig6a_0_NicheNet_runs.ipynb (then NicheNet_v3.ipynb)` (Spencer Watson; originally `analysis/survival/NicheNet_v3.ipynb`),
cells 162 (Fig 2d), 172 (grid search adjusted for MGMT + age) and 174 (Fig 2e forest), numbered as
in the EFS copy `Final_Notebooks/NicheNet_v3.ipynb`. Cells 164–170 and 173 (exploratory) not copied.

Changes:
- Inputs with the regular loaders: clinical rows `AllDataCenters.load_clinical(exclude_idh1_mut=False)`
  (= the notebook's cell-3 `load_processed_clinical` concatenation, 313 rows, index reset as
  `ignore_index=True` did); snRNA-seq obs `AllDataCenters.load_singlecell(obs_only=True)` (no
  exclusions, as cell 5 loaded it; the matrices are not used); curated labels
  `annotation.map_levels(efs_copy=True)` (the EFS key cell 13 read).
- All functions (`is_excluded`, `build_proportion_matrix`, `clean_corr`, `compute_linkage`,
  `plot_corr_clustermap_with_clusters`, `build_cluster_dict`, `print_clusters`,
  `get_clusters_at_threshold`, `compute_score`, `stratify`, `run_rmst_mgmt_age`,
  `apply_stratification`, `run_rmst_module`, `run_modules`, `fmt_p`, `get_vals`, `y_frac`) are in
  one cell after the config; cell 172's identical redefinitions, the unused `contains_tcell` and
  the `if not in dir()` guards removed. Imports and the three cells' config constants at the top.
- Cell 174's per-module loop is `run_modules(tau, clin_df)` (body unchanged) so it can be rerun at
  τ = 3 years for the †.
- Outputs through `output.save` / `output.save_table`, `OUT = output.setup("fig2de_celltype_modules")`,
  panels Fig 2d, Fig 2d modules, Fig 2e, Fig 2e dagger.
- Added: check against the saved output of cell 174 (published Fig 2e), count of fitted patients
  whose missing MGMT / sex / p53 the `(x == level)` coding counts as the reference level, the †
  at τ = 3 years, and the corrected versions for the response letter (survival cohort with
  complete-case covariates; modules recomputed on the 266 samples of the paper cohort, loaded
  with `AllDataCenters.load_singlecell(exclude_idh1_mut=True, exclude_patient_ids=True,
  exclude_relapse=True, obs_only=True)`).
- Header: "Changes" list moved here; decision codes (D8) and source-cell references removed from
  section titles and comments.

Package helpers used:
- `annotation.map_levels(efs_copy=True)` / `annotation.NEW_COLS`: the annotation-key cell (13) of
  NicheNet_v3, shared with the other notebooks that read the EFS key.
- `cohort.clinical()`, `cohort.survival_cohort()`, `cohort.code_covariates()`,
  `cohort.n_coded_as_reference()`: the 273 / 260 cohort definitions and covariate coding used for
  the printed checks and the corrected version (not in the original notebook; no regular loader
  defines the paper cohorts).
- `paper.survival` (import only: R_HOME, survival + survRM2), `paper.output`.

Verification (re-run 2026-11-20, 11 min): all 11 CSV outputs byte-identical to the previous run
(`figures/paper/05_celltype_modules/`); printed output identical apart from loader log lines.
Fig 2e matches the saved cell-174 output for all 7 modules (e.g. module 9 +2.02 mo, p = 0.0008);
273 fitted patients; † at τ = 3 y: modules 9, 8, 4, 2; corrected (b): 266 samples, 23 modules,
15 / 25 identical.

## fig3a_0_mutation_calling.py

Origin: `analysis/mutations/mutations.py` (WES patient x alteration matrix, frozen as
results/WES/mutations/patient_gene_matrix.csv, read by fig3a_oncoplot.py).

Changes:
- Header docstring (output, logic, inputs, runtime); change list removed.
- The gene lists are read next to this file (`GENE_LIST_DIR = Path(__file__).parent`; the source
  read them from the working directory, so it only ran from its own folder).
- Matrix written to `output.setup("pipelines_mutations")`, not over the frozen file.
- Dropped the IPython autoreload cell and the oncoplot-style heatmap (`create_oncoplot`,
  superseded by Fig 3a); unused imports removed.

Package helpers: `output` only.

Verification: re-run (17 min 38 s, S3 VCF reads). Not identical to the frozen matrix:
296 vs 299 samples (ER_G_001ult, ER_G_002ult, ER_G_003ult absent from the current S3 listing)
and MET_amp differs for 36 samples; all other 33 columns identical on the 296 common samples.
The run printed "Could not get coordinates for PIK3CA / TP53" (the gene coordinates come from an
online lookup inside mutation_utils), so the gene-coordinate step is not stable over time. The
code edits (paths, header) do not touch the logic; Fig 3a keeps reading the frozen matrix.

## fig3a_oncoplot.py

Origin: `analysis/mutations/mutations_updated_oncoplot.py` (Fig 3a and the alteration frequencies
of the text).

Changes (on the original file; plot functions in the file as written by the author):
- `_maybe_save` writes .pdf and .png with `output.save` (and to the panels folder when `panel` is
  given); `plot_oncoplot` takes `panel` (`"Fig 3a"` from `run_all_visualizations`).
- Oncoplot x label prints the number of plotted samples (was the hard-coded "Samples (n=261)",
  the number of clinical rows; 251 blocks are plotted).
- `RESULTS_DIR = output.setup("fig3a_oncoplot")` instead of figures/oncoplot.
- Prints the number of patients kept (`ids.patient_of`).
- New: alteration frequencies of the Results with three denominators (`TEXT_PERCENT` in Config;
  `gene_any` / `freq_table` with the other functions at the top), tables
  `alteration_frequencies_by_gene/_by_class.csv`; package versions.
- The functions are no longer in `gbmcollaboration/paper/wes.py` (used only here; the module
  was deleted). The co-occurrence heatmap and the original per-centre clinical loading loop
  (`load_tabular` + SURVIVAL_EXCLUSION on CR) are as in the original.

Package helpers: `ids.patient_of` (block -> patient for the printed patient count; shared ID
parser); `output`.

Verification: re-run (54 s). `oncoplot_top12.png`, `driver_frequencies_grouped.png` and both
frequency tables byte-identical to `figures/paper/06_oncoplot/`; 299 matrix samples, 261 clinical,
251 plotted blocks = 250 patients; EGFR 52.2, CDK6 25.9, PTEN 20.7, TERTp 19.1, CDK4 17.9,
NTRK2 16.3 % (text 52/26/21/19/18/16). `driver_cooccurrence_heatmap` is written again (not a panel).

## fig3b_0_infercnv_sc.py

Origin: `analysis/infercnv/infercnv_sc_analyses.py` (single-cell inferCNV; frozen outputs in
results/scRNAseq/infercnv/).

Changes:
- Header docstring (outputs, readers, runtime, comparison with the frozen outputs); change list removed.
- `RESULTS_DIR = output.setup("pipelines_infercnv_sc_analyses")` instead of the frozen EFS folder
  (Part 2 reloads Part 1's object from it).
- Part 2 subset on `obs["tumor_lineage"]`, a column the file never creates (KeyError): it is
  set from `lineage_status == "tumor-lineage"` before the subset (on the frozen obs this gives
  exactly the 782,542 frozen cells).
- `display(dominant_df)` (notebook-only) -> `print`; unused imports removed.
- Loading: the original `AllDataCenters.load_singlecell(...)` (full object, needed for inferCNV).

Package helpers: `output` only.

Verification: not re-run (hours, well over 60 GB; outputs frozen). `py_compile` OK.

## fig3b_0_cnv_wes_matrix.py

Origin: `analysis/infercnv/generate_cnv_wes_matrix.py` (WES CNV matrix, frozen as
results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv).

Changes:
- Header docstring (output, inputs, runtime, comparison with the frozen matrix); change list removed.
- The matrix is written to `output.setup("pipelines_generate_cnv_wes_matrix")`, not over the frozen
  file. `RESULTS_DIR` / `PATH_BINS_SINGLECELL` as in the original.
- Dropped Part 2 (exploration: PCA, outlier QC, single-sample plots).

Package helpers: `output` only.

Verification: not re-run (73 min of S3 reads, needs AWS credentials; frozen output). `py_compile` OK.

## fig3b_suppfig3ac_infercnv_wes_analyses.py

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

## fig3c_suppfig3b_wes_profile_os_rmst.ipynb

Origin: `analysis/survival/wes_profile_os_rmst.ipynb` (Spencer Watson): Fig 3c, Supp Fig 3b.

Changes:
- Clinical rows from the notebook's own `load_clinical_data()` (per-centre `load_tabular` with
  the loader exclusions, SURVIVAL_EXCLUSION dropped on CR; 261 rows), now called (the call was
  commented out in the original, which relied on an earlier session state).
- `COVARIATE_COLS`: "p53" replaces "Site code" (Spencer's revision of the model); p53 coded
  Negative 0 / Positive 1 in the same `mapping` dict. `load_and_merge_cluster_data` and
  `run_rmst_and_km_cluster_vs_rest` take the covariate list as an argument (default
  `COVARIATE_COLS`) so the published model runs with the same code.
- `plot_global_km_curves` gave `kmplot` OS in years, unsigned, with `convert_days_to_months=True`:
  kmplot cast the durations to whole years and read the sign as the event. It now passes the
  signed OS in days. The global log-rank test is unchanged.
- KM confidence bands are rasterized in the PDF (`rasterize_ci_bands`), so they keep their
  transparency (they rendered as solid blocks in the manuscript PDF).
- Removed the per-cluster Cox models and the categorical Cox forest (not in the paper).
- Outputs: `OUT = output.setup("fig3c_suppfig3b_wes_profile_os_rmst")`; the data saved with the KM
  figures are the KM curves; forest plot saved with `output.save` (panel Supp Fig 3b); Fig 3c
  panel tag on cluster 2 vs rest.
- New cells: the published model, a table of both models next to the published numbers, versions.
  `row` (model-table helper) is in the functions cell.
- Removed decision codes (D9) and review comments; header without the change list.

Package helpers: `ids.patient_of` (printed patient count of the joined rows); `output`.

Verification: re-executed (1 min). All printed numbers identical to the previous run: 261
clinical rows, 253 joined (252 patients); model age+MGMT+sex+p53 n = 173, cluster 2 vs rest
+5.08 mo [1.52, 8.64], p = 0.0052; published model n = 247, +3.35 mo [0.28, 6.43], p = 0.0326
(= published). All KM / forest PNGs and summary CSVs byte-identical to
`figures/paper/10_wes_profile_os_rmst/` (fig3c_cluster2_models.csv differs only in the model label).

## fig3d_correlation_cnv_wes_ecm.py

Origin: `analysis/multimodal_correlation/correlation_cnv_wes_ecm.py`, Parts 1-2 (Fig 3d, PERMANOVA
numbers of the Results/Methods).

Changes:
- Single-cell data: `ad.AnnData(obs=AllDataCenters.load_singlecell(exclude_idh1_mut=True,
  exclude_patient_ids=True, exclude_relapse=True, obs_only=True))` (the original loaded the whole
  expression cohort; only obs is used). The original `adata.X = layers[...]` lines are gone with X.
- Annotation key read from the original EFS `SingleCell_Annotation_Key_v3.csv` as in the source
  (no package helper); sample IDs `str[:9]` as in the source.
- PERMANOVA: scikit-bio is not in the environment; `permanova()` (numpy/scipy Bray-Curtis
  PERMANOVA, same pseudo-F, seeded permutations) is defined at the top of this file (was
  `gbmcollaboration/paper/wes.py`, used only here). `stars` moved to the top too.
- New: level-3 PERMANOVA next to level 4 and the published values (`permanova_cnv_cluster.csv`);
  check of the Kruskal-Wallis table against the published `figures/figure3/fig3e.csv`.
- Outputs: `RESULTS_DIR = output.setup("fig3d_correlation_cnv_wes_ecm")`, Fig 3d panel
  (`fig3d_celltype_composition_cna_clusters`; the source saved it as figures/figure3/fig3e.pdf).
- Dropped (not in the paper): ECM Kruskal-Wallis and mean-proportion prints, the stacked barplot,
  the AC/OPC ratio test (KeyError on the level-4 table), unused parameters, Part 3 (pseudobulk
  DESeq2 and its dot plot).
- Header without the change list.

Package helpers: `output` only.

Verification: re-run (1 min 05 s). `fig3d_celltype_composition_cna_clusters.png/.csv`,
`kruskal_wallis_test_cellular_landscape_cluster.csv` and `permanova_cnv_cluster.csv`
byte-identical to `figures/paper/09_correlation_cnv_wes_ecm/`. PERMANOVA ECM F = 0.582 p = 0.868
(published 0.58 / 0.86), level 4 F = 4.02, level 3 F = 4.56 (published 4.56), p = 0.001; KW
table = published fig3e.csv (71/71 types, max |dH| 4e-15, same order).

## suppfig3d_wes_cluster_covariate_rmst.ipynb

No original code found. Model = Fig 1b / Supp Fig 1c clinical RMST of fig1b_suppfig1_GBM_MS_Summary (one row per subject, τ = 3 y, each covariate adjusted for the others + centre dummies), on the 264 CNV-clustered subjects, without and with a cluster-2 indicator; plus all 287 subjects (cluster 2 = 0 without WES).
`rmst_adjusted` is a compact version of fig1b's `run_rmst` adjusted branch (arm missing -> dropped; covariates missing -> reference level, as fig1b).
Verification: all-subject model reproduces Fig 1b (+7.1 / −5.4 / +3.1 / +3.7 months). With cluster 2: MGMT +6.7 -> +6.3, age −4.9 -> −4.9, sex +3.5 -> +3.2 (all p < 0.011, unchanged significance); p53 +3.0 (p 0.082) -> +3.7 (p 0.033).

## suppfig3ef_analyses_wes_sc_infercnv.py

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

## fig4ab_suppfig4abg_Cell_Neighborhoods.ipynb

Origin: `analysis/survival/Cell_Neighborhoods.ipynb` (Spencer Watson).

Changes:
- Loading cells (2–12) replaced: `neighborhoods.load_spots()` (spots and coordinates with the
  loader exclusions; no cell uses the expression matrix, so it is not loaded),
  `AllDataCenters.load_cell2location(visium.sections())` merged with `visium.add_c2l_to_obs`,
  `annotation.add_celltype_key(efs_copy=True)` (the EFS key the notebook read). The clinical table
  (cell 4) was not used and is not loaded.
- Clustering: cell 36 (and the Hellinger representation it needs, first part of cell 14) run only
  with `RUN_CLUSTERING = True` and write to `OUT`; otherwise the frozen annotation and centroid files
  are read (`config.CT_ANNOTATIONS_CSV`).
- Removed the first strategy (Leiden clusters then kNN, superseded by cell 36): Leiden/UMAP part of
  cell 14, cells 16–25 (cell 22 also wrote `celltype_neighborhood_annotations_k10.csv`, with labels
  not sorted by size), 29 (earlier version of cell 38), 31 (earlier version of cell 46), 33
  (cell-type composition bars), 34 (exploratory); cell 1 (versions → `output.versions`).
- Cells 38 and 40 no longer reload the ECM file (cell 27 does).
- Cells 42 and 44: the cell2location → curated-type block is `neighborhoods.c2l_curated_counts`;
  the composition / specificity / marker-score matrices are `marker_scores` (same code, now defined
  in the notebook: `neighborhoods.marker_scores` removed from the package, used only here).
- Functions (`make_color_map`, `global_categories`, `color_array`, `joint_freq_pct`,
  `marker_scores`) in one cell after the parameters.
- Figures go to `OUT = output.setup("fig4ab_suppfig4abg_Cell_Neighborhoods")` with `output.save`
  (Fig 4a, 4b, Supp Fig 4a, 4b, 4g also to the panels folder); cell 46's matrix is compared with the
  one of the last run instead of overwriting it in `Analysis_Repo/Notebooks`.
- Fig 4a: the notebook's subset dotplot (cell 44: curated names, niches 0 and 6 left out) is not
  the published panel and no longer carries the panel tag. Added a reconstruction cell (`panel="Fig 4a"`):
  cell2location means with the uncurated scANVI names (names containing doublet / other / unknown
  removed), `dotplot_composite_score(abund_per_nh, included_niches=["0","1","2","3","4","5","7","8","9"],
  top_n_per_niche=5)`. It gives the published 36 cell types plus Mural_LEPR and Plasma_nos. A second
  version with `DROP = ["Mural_LEPR", "Plasma_nos"]` (`panel="Fig 4a published set"`, labelled as
  inferred) gives exactly the published 36, same order except near-ties (Neuron_neurogliaform with
  CT_0 instead of CT_9; Fibro_CAF-like / Mural_perivasc_fibro, Pericyte / Macro_LYVE1,
  Neuron_inh_PVALB / Neuron_exc_CDH24 swapped). Inference: the published panel was drawn from a
  cell2location run without those two factors (not on S3 / EFS).
- Per-cell imports moved to the top. Header: "Changes" list moved here.

Package helpers used:
- `gbmcollaboration.viz.cell2loc_plotting.dotplot_composite_score` (existing package plot function).
- `neighborhoods.load_spots()`: spot coordinates via the regular `AllDataCenters.load_spatial(load_gex=False)`;
  `load_spatial_adata(genes=...)` returns no coordinates.
- `neighborhoods.c2l_curated_counts`: shared with fig4cdg.
- `visium.sections()`, `visium.add_c2l_to_obs`, `annotation.add_celltype_key(efs_copy=True)`,
  `paper.config`, `paper.output`.

Verification (re-run 2026-11-20, 10 min): all 6 CSV outputs (Fig 4a, 4b, Supp Fig 4a, 4b, 4g)
byte-identical to the previous run (`figures/paper/12_Cell_Neighborhoods/`). Global NMI 0.140,
per-sample median 0.237; Supp Fig 4b matrix vs the notebook's last run: 228 samples, same order,
max |Δ| = 1.1e-16. Only printed difference: `T_other` is now before `total_cells` in `adata.obs`
(column order of `load_cell2location`), with no effect on any output.
Re-run with the Fig 4a reconstruction (10:00, 10 min): the 6 CSVs above unchanged; Fig 4a 38 types
(extra Mural_LEPR, Plasma_nos, none missing); Fig 4a published set: 36 types, same set as published.

## fig4cdg_suppfig4f_ECM_niches_v4.ipynb

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

## fig4e_suppfig4e_colocalization_plots.ipynb

Origin: `analysis/cell2loc/colocalization_plots.ipynb` (Owkin, j-bac).

Changes:
- Paths from `gbmcollaboration.paper.config`; figures and tables go to
  `OUT = output.setup("fig4e_suppfig4e_colocalization_plots")` with `output.save` /
  `output.save_table` (Fig 4e, Supp Fig 4e, Supp Tables 3b–c also to the panels folder) instead of
  EFS `figures/visium/` and `supplementary_tables/`; the published versions are read back and compared.
- Loading (cell 7): `neighborhoods.load_spots()` gives the q05 abundances that
  `load_spatial_adata(load_deconvolution=True)` put in `adata.obs`, without the expression matrix;
  cell 9's gene means per neighborhood come from `neighborhoods.pseudobulk_sums` (same means).
- Fig 4e (cell 16): run with the defaults of `dotplot_composite_score` (the published version;
  the cell last ran with `center_zero=True, cmap='bwr'`).
- Bug fixed: `col_prefix` kept the "ECM_" of the previous annotation for `pathology_annotations`
  (cells 12–24); now empty for every annotation but the ECM neighborhoods.
- Removed: cell 2 (FileLink), 8 (duplicate of 7), 12 (clustermap dotplots), 14 (earlier version of
  16), 18 (`plot_neighborhood_marker_dotplot_from_agg` not defined: NameError), 20 and 22
  (pathology × neighborhood dotplots, not paper panels; Fig 4f is the residual heatmap of
  fig4f_pathologist_annotation_analyses.py), unused imports.
- Header: "Changes" list moved here.

Package helpers used:
- `neighborhoods.load_spots()`: the q05 deconvolution table of the full `load_spatial_adata`, read
  with the regular `AllDataCenters.load_spatial(load_gex=False)` (the full object with expression
  is > 100 GB; `load_spatial_adata(genes=...)` does not return the deconvolution).
- `neighborhoods.pseudobulk_sums`: gene sums per section × ECM neighborhood with each centre's
  regular `load_spatial`, a few sections at a time (Supp Table 3b needs every gene; shared with fig4cdg).
- `gbmcollaboration.utils.pathology_annotation_loader` (existing package module, used by the original).
- `paper.config`, `paper.output`.

Verification (re-run 2026-11-20 under the heavy lock, 13 min): all 28 CSV outputs byte-identical
to the previous run (`figures/paper/13_colocalization_plots/`). Against the published files:
Supp Table 3b (10 × 18,082, same genes) max |difference| 1.8e-15; Supp Table 3c (12 tables)
≤ 1e-16; Fig 4e plotted abundances and Supp Fig 4e bar data (228 × 10) max |difference| 0.

## fig4f_pathologist_annotation_analyses.py

Origin: `analysis/hne/pathologist_annotation_analyses.py` (Owkin).

Changes:
- The annotation loading (DataConfig, AnalysisConfig, S3 reading, standardisation, exclusions,
  spot IDs, centre) comes from `gbmcollaboration.utils.pathology_annotation_loader`, of which the
  script held a copy; `load_and_prepare` merges its output with the ECM table as before.
- ECM table path from `gbmcollaboration.paper.config.ECM_ANNOTATIONS_CSV`.
- Output folder `OUT = output.setup("fig4f_pathologist_annotation_analyses")` instead of
  `results/hne/ecm_pathology_figures/`; Fig 4f also written to the panels folder (`panel="Fig 4f"`).
- Removed the figures not in the paper (agreed with the user): the three-panel centre
  replication (`fig_three_panel_centers`), the dot plot (`fig_dotplot`) and the sparse top-effects
  heatmap (`fig_sparse_top_effects`), with their calls; `savefig` import no longer needed.
- Added at the end: the Fig 4f table is read back and compared with the published CSV
  (`results/hne/ecm_pathology_figures/Fig_residual_heatmap.csv`).
- Per-file "Changes" list removed from the docstring (this file).

Package helpers used:
- `gbmcollaboration.utils.pathology_annotation_loader`: existing package module (committed,
  not paper-specific) holding the annotation loading the script duplicated.
- `paper.config` (paths), `paper.output` (output folder, panels).

Verification (re-run 2026-11-20, ~1 min): Chi2 = 77872.63 (dof 81), Cramér's V = 0.152;
batch-adjusted LR = 130.33 (df 81), p = 4.2e-4; Jaccard UKER vs CHUV = 0.781.
Fig 4f: same rows / columns / order as the published table, max |difference| = 0.

## fig4h_suppfig4cd_neighborhood_survival.ipynb

Origin: Spencer Watson, notebook Updated_figures (local, not in the repository): the RMST sweep under
"Updated ECM Neighborhood Plot" and the tau = 1.5 years check after it, as reproduced in
`analysis/scratch/review_2026-11/source_cells_spencer_2026-11.ipynb` cells 6–9.

Changes (port and this pass):
- The sweep cell is split into its Config block (Config cell), its functions (`clinical_frame`,
  `apply_split`, `run_rmst`, `fmt_p`, `run_sweep`, `forest_plot`, `dagger`, `describe_fit`; one
  functions cell after the Config) and the calls, so the same sweep runs on the ECM neighborhoods
  (Fig 4h) and on the cell-type neighborhoods (Supp Fig 4c–d). Analysis code as in the source cells.
- This pass: the functions were in `gbmcollaboration/paper/neighborhood_rmst.py`; they are used only
  by this notebook, so they are back at its top and the module is deleted.
- Clinical table: `AllDataCenters.load_clinical(exclude_idh1_mut=False)` in place of the local CSV
  export of the processed clinical tables (`clinical_frame`).
- Neighborhood labels from the frozen per-spot annotation CSVs (`config.ECM_ANNOTATIONS_CSV`,
  `config.CT_ANNOTATIONS_CSV`).
- Added: `describe_fit` (sections without a subject, fitted patients outside the 260-patient
  survival cohort, covariates coded as the reference level), the ECM_3 / ECM_1 print, and the
  CT_3 / CT_0 mutual adjustment table.
- Outputs: `OUT = output.setup("fig4h_suppfig4cd_neighborhood_survival")`, panels `Fig 4h`, `Supp Fig 4c`.
- Header: change narration and decision codes removed.

Package helpers used:
- `cohort.survival_cohort()` in `describe_fit`: the paper's 260-patient survival cohort definition
  (loader exclusions + the 13 Charité patients), shared by all survival files.
- `paper.config` (paths, R_HOME), `paper.output`.

Verification (re-run 2026-11-20, 56 s): all printed output identical to the previous run (ECM_3
adjusted −3.50 mo, p = 0.0114; ECM_1 0.2595; ECM_9 0.4743; † ECM_3 −1.7 mo p = 0.0129; CT_3 −3.91,
p = 0.0057; CT_0 +3.30, p = 0.0257; mutual adjustment table). fig4h / suppfig4c CSVs and PNGs and
the mutual-adjustment CSV byte-identical to `figures/paper/15_neighborhood_survival/`.

## suppfig4d_CT0_CT3_composition.ipynb

No original code found. Implements Methods "Neighborhood cell type composition and marker scoring" (CT_0/CT_3 log2FC vs mean of the other eight neighborhoods, union of top-5 enriched/depleted with ≥ 0.3 % composition; compartment pp CT_0 vs all other neighborhoods).
Data: `AllDataCenters.load_cell2location(sections, abundance="q05_cell_abundance_w_sf")` (identical to the deconv_cohort_level_4 obs of load_spatial_adata), frozen k = 10 CT labels.
Verification: lymphoid + dendritic +7.09 pp (text +7.1). Mesenchymal −6.7 pp not reproduced (MES-like −8.3, MES-like + MES2-like −11.6, malignant −16.2).

## fig5_0_TIL_zones.ipynb

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

## fig5a_TIL_zones_map.ipynb

Origin: `_plot_til_zones_ax` / `_til_percent_by_sample` of EFS Analysis_Repo/Notebooks/TIL_Niche_Comparisons_v2_clean.ipynb cell 26 (published style, Spencer Watson); `plot_visium_til` of fig6a_0_NicheNet_runs.ipynb (then NicheNet_v3.ipynb) cells 184–185.

Changes:
- Sections set explicitly: ER_G_631a_vis (TIL-high) and CR_G_310a_vis (TIL-low), identified by matching the published TIL+ islands; the facet cell's automatic top/bottom selection picks ER_G_740a/736a and CR_G_322a/310a with the frozen zones (printed).
- `_plot_til_zones_ax` takes `adata` as an argument; coordinates from obs only; returns coords for the new `add_scale_bar` (1 mm, pixel size from the 100 µm spot pitch).
- `_til_percent_by_sample`: spot count aggregates `("_til", "size")` (the original `(key, "count")` raises KeyError with as_index=False in pandas 2.2).
- `plot_visium_til` verbatim except: no `plt.show()` inside (so the figure can be saved), NearestNeighbors import at top.
- Data: `DataCenter.load_spatial(list_sample_id=..., load_gex=False)` per centre (AllDataCenters.load_spatial with a list fails for centres without listed sections) + frozen `config.TIL_ZONES_CSV`.

Package helpers: `config` (paths), `output`.
Verification: TIL+ counts 660 (ER_G_631a) / 7 (CR_G_310a); maps match the published Fig 5a (published panels are crops).

## fig5b_TIL_survival.ipynb

Origin: Spencer Watson, TIL_Niche_Comparisons (EFS `Analysis_Repo/Notebooks/`), cell "RMST Survival
Analysis", as reproduced in `analysis/scratch/review_2026-11/source_cells_spencer_2026-11.ipynb` cells 16–17.

Changes (port and this pass):
- Spot table from the frozen TIL zone file (`til.til_zones()`) in place of `adata.obs` of the full
  Visium object (the cell only reads `sample` and `TIL_positive`).
- Clinical table: `AllDataCenters.load_clinical(exclude_idh1_mut=False)` (313 rows), the table the
  cell's `clin_out` was; this pass replaces the removed `cohort.clinical_as_ran()`.
- RMST through `rmst()` (survRM2 via rpy2, unadjusted + adjusted, results in months); this pass
  moves `rmst` and `to_r` from `gbmcollaboration/paper/survival.py` into the notebook (used only
  here) and removes `rmst` from the module.
- Figure settings in the Parameters cell; KM plot code as in the source cell; output folder
  `OUT = output.setup("fig5b_TIL_survival")`, panel `Fig 5b`.
- Added: printed checks of the join (sections without a clinical match, missing MGMT, sections
  without TIL+ spots); TIL+ fraction vs age and by MGMT; the sensitivity table (other covariate sets,
  without probe-batch blocks, corrected cohort).
- Bug fixed in an added check (this pass): "Fitted patients outside the paper cohort" compared
  9-character join keys with block IDs and the full cohort, so it always printed 0. It now lists
  the fitted patients outside the 260-patient survival cohort: 10 Charité patients (the same 10 as
  `describe_fit` in fig4h).
- Header: change narration and decision codes removed; comments neutral.

Package helpers used:
- `til.til_zones()`, `til.til_groups()`: the frozen TIL zone table and the TIL-high / TIL-low
  sections, shared with fig5c, fig5de and data_release.
- `cohort.code_covariates()`, `cohort.survival_cohort()`: covariate coding (as reference / complete
  cases) and the 260-patient survival cohort, shared by the survival files.
- `ids.parse_ids`, `paper.config` (R_HOME), `paper.output`.

Verification (re-run 2026-11-20, 58 s): all printed numbers identical to the previous run (n = 207,
54 / 153, +3.93 mo [0.5, 7.4], p = 0.0251; ρ = −0.22; sensitivity table) except the corrected
check line above. fig5b CSVs and PNG byte-identical to `figures/paper/17_TIL_survival/` (PDF differs
only in its creation date).

## fig5c_suppfig5ab_TIL_supplement.ipynb

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

## fig5de_suppfig5cde_ECM_supplement.ipynb

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

## fig5fg_TIL_zone_GSEA.ipynb

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

## supptables_4a_ecm_til_zone_tables.ipynb

Origin: Spencer Watson, EFS Analysis_Repo/Notebooks/TIL_Niche_Comparisons_v2_clean.ipynb cells 70–71
(`ecm_neighborhood_TIL_group_by_zone_pseudobulk`, written 2025-12-12) and the cell-72 matrisome gene list;
`ecm_neighborhood_niche_matrisome_summary` (written 2025-12-10): no code found, reconstructed from its columns
(spot-level fractions per zone / TIL group, float64 mean log-normalised expression per neighborhood).
Inputs the published tables used (AS_PUBLISHED = True): the pre-v4 ECM labels
`Analysis_Repo/Notebooks/ecm_neighborhood_annotations.csv` and 229 sections = the paper's 228 + CR_G_335a_vis
(added to PATIENT_EXCLUSION on 2025-12-12); its TIL zones are recomputed with the fig5_0_TIL_zones rule
(`section_til_zones`, used only here); groups 58 TIL-high / 171 TIL-low. The paper-cohort version
(228 sections, 57/171, v4 labels) is also written (`*_paper_cohort_v4.csv`).
Data: `til.til_zones()`, `CRDataCenter.load_spatial(load_gex=False)` + `load_cell2location` for CR_G_335a_vis,
`AllDataCenters.load_spatial_adata(genes=96 matrisome genes, exclude_patient_ids=False)`.
Verification: both as-published tables byte-identical to the published files. Runtime 12 min (~35 GB).

## supptables_4a_til_zone_celltype_stats.ipynb

No code found (repo history all branches, EFS Final_Notebooks, Analysis_Repo, other EFS notebooks); the tables were written 2025-12-17 to EFS results/visium/til_zones/. Reconstructed: mean of log1p(cell2location means abundance) over the zone's spots vs the other spots of the same TIL group's sections (float32, NaN where a type is absent from a centre's reference), log2 ratio with 1e-6 pseudocount; summary = top 30 curated names by mean_in and by log2 ratio.
Data: `AllDataCenters.load_cell2location(visium.sections())`, `til.til_zones()`, `til.til_groups()` (57/171), `annotation.celltype_key(efs_copy=True)`.
Outputs: SuppTable4a_til_zone_celltype_stats_long.csv, SuppTable4a_til_zone_celltype_summary.csv.
Verification: both byte-identical to the published files. Runtime ~1 min.

## supptables_4b_gsea_TIL_zones.ipynb

Origin: Spencer Watson, EFS Analysis_Repo/Notebooks/TIL_Niche_Comparisons_v2_clean.ipynb cells 60 (gene-set library) and 61 (pseudobulk, Welch / paired t rankings, gseapy prerank), run 2026-02-25; fig5fg_TIL_zone_GSEA.ipynb only reads the frozen table.
Changes: per-section zone means from `neighborhoods.pseudobulk_sums` (regular per-centre `load_spatial`, a few sections at a time, gene space of all 228 sections) instead of the in-memory `load_spatial_adata` matrix (> 100 GB); image coordinates `imagerow`/`imagecol` added as two "genes" (`COORDS_IN_RANKING = True`) because the loader of February 2026 returned them in the matrix (18,084 vars); `processes=` -> `threads=` (gseapy deprecation, results identical); per-comparison CSVs not written (only the combined table); helper prints moved to `print_top`.
Package helpers: `neighborhoods.pseudobulk_sums` (shared with fig4e, avoids the > 100 GB matrix), `neighborhoods.load_spots` (coordinates), `til.til_zones/til_groups`, `visium.sections`.
Verification: 47,957 rows, same columns; ES, lead genes identical; NES/p/FDR identical in 6/7 comparisons, TIL+ between-group max |ΔNES| 7e-4, |ΔFDR| 2e-4 (tied-score gene order); same 21,442 FDR<0.05 calls. Runtime 40 min, < 25 GB.

## fig6a_0_NicheNet_runs.ipynb

Origin: Spencer Watson, `analysis/survival/NicheNet_v3.ipynb` (moved with `git mv`).

Changes (this pass, first markdown cell only): references to the files that run cells 162–174 and
the circos cells updated to `fig2de_celltype_modules.ipynb` and `fig6a_suppfig6_NicheNet_circos.ipynb` (same folder after the rename);
the reference to the removed `gbmcollaboration.paper.nichenet` module and the edit narration removed.

Not re-run (pipeline; outputs of the original run kept).

## fig6a_suppfig6_NicheNet_circos.ipynb

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

## fig6bg_xenium_niches.ipynb

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
