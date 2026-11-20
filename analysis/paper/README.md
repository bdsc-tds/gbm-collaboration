# Analyses of Watson et al., MOSAIC GBM

Every panel, supplementary table and number of the manuscript comes from one file
in this folder. Files run top to bottom in the `gbm_main` environment; numbers follow the
manuscript order. Upstream steps that produced the frozen inputs (clustering, inferCNV,
NicheNet, mutation calling) are in `pipelines/`; they are slow and are not re-run for the
figures. Code not used in the paper is in `../../scratch/`, in its original folders.

Most files are the authors' original notebooks and scripts, moved here with `git mv` and
edited as little as possible: data loading through the package, exploratory and dead cells
removed, bugs fixed. Each file's first cell says what it produces and lists every change.
Read `git log --follow -p <file>` for the diff against the original.

## Shared code (`gbmcollaboration/paper/`)

Code that several analyses repeated now lives in the package, so a file only contains its
own analysis:

| Module | Replaces |
|---|---|
| `config` | paths of frozen inputs, loader flags, constants (TIL threshold, τ) |
| `ids` | patient / block / section parsing (keeps backup sections `_vbu`, `_vb2`, `_bu2`, `_bu3` attached to their block) |
| `cohort` | the clinical loading cell (`clinical_as_ran()`, every row as the notebooks loaded it) and the 273-patient cohort (`clinical()`, `survival_cohort()`) |
| `visium` | `load_spatial_adata` (batched, gene subset), the cell2location loader and merge (`load_c2l`, `add_c2l_to_obs`) |
| `annotation` | the annotation-key cell (`celltype_key`, `add_celltype_key`) |
| `til` | TIL⁺ spots, TIL zones, TIL-high / TIL-low sections (Spencer's TIL_supplement helpers) |
| `survival` | survRM2 RMST, Kaplan–Meier, Cox PH test |
| `neighborhood_rmst` | the RMST sweep of Updated_figures (Fig 4h, Supp Fig 4c) |
| `clinical_survival` | the PH-diagnostics cell of GBM_MS_Summary (Supp Fig 1a–b), run on two cohorts |
| `neighborhoods` | ECM / cell-type neighborhood plotting helpers shared by 11–13 |
| `snrna` | snRNA-seq loading and plotting helpers shared by 03–05 |
| `wes` | oncoplot functions, obs-only h5ad reader, GTF overlap, PERMANOVA |
| `nichenet` | the circos cell of NicheNet_v3 (Fig 6a, Supp Fig 6) |
| `xenium` | tiling and TIL–tile Spearman helpers of niches.ipynb |
| `output` | output folder per file, figure saving, package versions |

Cohort exclusions are the lists in `gbmcollaboration/definitions.py`, applied by the
DataCenter loaders and, for frozen tables that do not go through a loader, by
`gbmcollaboration.apply_exclusions`. To remove a patient from every analysis, add them to
`PATIENT_EXCLUSION` there.

Outputs go to EFS `figures/paper/<file name>/` (set `PAPER_OUT` to change). Every
manuscript panel and supplementary table is also written to **`figures/paper/panels/`**
under its manuscript name (`Fig4h.pdf/.png/.csv`, `SuppFig5a…`, `SuppTable3b.csv`; panels
drawn as several files carry a suffix, e.g. `Fig6b_hBrain`): `output.save(..., panel="Fig 4h")`.
Running every file fills that folder with the whole paper.

## Files

| File | Panels and numbers | Origin |
|---|---|---|
| `table1bc_cohort_description.py` | Table 1b–1c | analysis/cohort_desc/cohort_description.py |
| `table1a_sample_number_per_modality.ipynb` | Table 1a | analysis/cohort_desc/ |
| `fig1b_suppfig1_GBM_MS_Summary.ipynb` | Cohort numbers, Fig 1b, Supp Fig 1a–c | analysis/survival/ (Spencer) |
| `fig2a_suppfig2_figure2.py` | Fig 2a, Supp Fig 2 UMAPs and composition bars | figures/figure2.py |
| `fig2bc_suppfig2_marker_dotplots.ipynb` | Fig 2b–c, Supp Fig 2 dotplots | analysis/annotation/ |
| `fig2de_celltype_modules.ipynb` | Fig 2d–e | pipelines/NicheNet_v3.ipynb cells 162–174 (Spencer) |
| `fig3a_oncoplot.py` | Fig 3a, alteration frequencies | analysis/mutations/mutations_updated_oncoplot.py |
| `fig3b_suppfig3ac_infercnv_wes_analyses.py` | Fig 3b, Supp Fig 3a, 3c, CNV scores | analysis/infercnv/ |
| `suppfig3ef_analyses_wes_sc_infercnv.py` | Supp Fig 3e–f | analysis/infercnv/ |
| `fig3d_correlation_cnv_wes_ecm.py` | Fig 3d | analysis/multimodal_correlation/ |
| `fig3c_suppfig3b_wes_profile_os_rmst.ipynb` | Fig 3c, Supp Fig 3b | analysis/survival/ |
| `fig4cdg_suppfig4f_ECM_niches_v4.ipynb` | Fig 4c, 4d, 4g, Supp Fig 4f, Supp Table 3a | analysis/survival/ (Spencer) |
| `fig4ab_suppfig4abg_Cell_Neighborhoods.ipynb` | Fig 4a–b, Supp Fig 4a, 4b, 4g | analysis/survival/ (Spencer) |
| `fig4e_suppfig4e_colocalization_plots.ipynb` | Fig 4e, Supp Fig 4e, Supp Tables 3b–c | analysis/cell2loc/ |
| `fig4f_pathologist_annotation_analyses.py` | Fig 4f | analysis/hne/ |
| `fig4h_suppfig4cd_neighborhood_survival.ipynb` | Fig 4h and †, Supp Fig 4c–d | Updated_figures (Spencer) |
| `fig5c_suppfig5ab_TIL_supplement.ipynb` | TIL groups, Supp Fig 5a–b, Fig 5c | TIL_supplement (Spencer) |
| `fig5b_TIL_survival.ipynb` | Fig 5b, TIL vs age, sensitivity table | TIL_Niche_Comparisons cell 32 (Spencer) |
| `fig5de_suppfig5cde_ECM_supplement.ipynb` | Fig 5d–e, Supp Fig 5c–e, Supp Table 4a gene table | ECM_supplement (Spencer) |
| `fig5fg_TIL_zone_GSEA.ipynb` | Fig 5f–g | TIL_Niche_Comparisons_v2_clean cell 63 (Spencer) |
| `fig6a_suppfig6_NicheNet_circos.ipynb` | Fig 6a, Supp Fig 6a–b | pipelines/NicheNet_v3.ipynb circos cells (Spencer) |
| `fig6bg_xenium_niches.ipynb` | Fig 6b, 6g | analysis/xenium/niches.ipynb |

Not generated by code: Fig 1a, 1c (schematics) and Fig 6c–f (microscopy).

Not yet in a file:
- **Fig 5a** (TIL zones on ER_G_672a_vis): `plot_visium_til` in `pipelines/NicheNet_v3.ipynb`
  cells 177–185.
- **Supp Fig 3d** (clinical RMST with and without CNV cluster 2) and **Supp Fig 4d** (CT_0 vs
  CT_3 composition, log2 fold change): no code found in the repo or on EFS; needs the author.
- **Fig 4a** runs but does not reproduce the published panel (see below).

| `pipelines/` | Frozen output |
|---|---|
| `infercnv_sc_analyses.py` | results/scRNAseq/infercnv/ (genomic windows, malignant-only object) |
| `generate_cnv_wes_matrix.py` | results/WES/cnv_matrix/cnv_wes_matrix_full_cohorts.csv |
| `mutations.py` | results/WES/mutations/patient_gene_matrix.csv |
| `TIL_zones.ipynb` | results/visium/til_zones/til_zone_annotations.csv |
| `NicheNet_v3.ipynb` | zone DE, NicheNet runs, nichenet_analysis/*/axes_selected.csv |

The ECM and cell-type neighborhood labels are computed in `fig4cdg_suppfig4f_ECM_niches_v4.ipynb` and
`fig4ab_suppfig4abg_Cell_Neighborhoods.ipynb` (clustering cells, off by default: `RUN_CLUSTERING = False`
reads the frozen labels); the WES CNV clusters in `fig3b_suppfig3ac_infercnv_wes_analyses.py`.

## Review answers (D1–D18, Spencer, November 2026)

Where each answer is implemented, with the numbers the scripts print.

| Item | Decision | Implemented in | Result |
|---|---|---|---|
| D1 Fig 5b cohort | Keep the published model; send the corrected cohort for the response letter | `fig5b_TIL_survival.ipynb` | As ran: n = 207 (54 / 153), +3.93 mo, p = 0.0251 (21 backup sections lost at the join, 5 missing MGMT counted unmethylated). Corrected: n = 212 (55 / 157), +4.00 mo [0.6, 7.4], p = 0.022 |
| D2 TIL and age | Fig 5b stays MGMT-adjusted; text adds ρ and the 4-covariate estimate | `fig5b_TIL_survival.ipynb` | ρ = −0.22, p = 0.001 (207 sections); MGMT + age + sex + p53: +2.90 mo [−0.3, 6.1], p = 0.078; sensitivity table incl. probe-batch rows |
| D3 57 / 153 / 18 | TIL-low = 171 incl. the 23 without TIL⁺ spots; zone comparisons 57 vs 148 | `gbmcollaboration/paper/til.py`, `fig5c_suppfig5ab_TIL_supplement`, `fig5de_suppfig5cde_ECM_supplement` | 18 = 21 TIL-negative sections in the Fig 5b model minus 3 without MGMT (printed) |
| D4 Fig 4h source | Updated_figures ECM sweep, ported as it ran | `fig4h_suppfig4cd_neighborhood_survival.ipynb` | Reproduced: ECM_3 adj p = 0.0114, ECM_9 0.4743 (ECM_1 0.2595 vs published 0.2598). 10 fitted patients outside the survival cohort, 5 MGMT and 51 p53 values counted as the reference level (printed) |
| D5 ECM_3 / ECM_1 text | Numbers from the Fig 4h model | `fig4h_suppfig4cd_neighborhood_survival` | ECM_3 −3.5 mo [−6.2, −0.8], p = 0.011; † = ECM_3 only (−1.7 mo at 1.5 y, p = 0.013). **ECM_1 is not significant unadjusted in this model** (+2.7 mo, p = 0.083); the text's +3.0 mo, p = 0.042 came from the MGMT-only cell |
| D6 CT_3 / CT_0 | One model: the sweep | `fig4h_suppfig4cd_neighborhood_survival` | **The sweep as it ran (n = 207) gives both published numbers**: CT_3 −3.9 mo, p = 0.006; CT_0 +3.3 mo [0.4, 6.2], p = 0.026. The review's CT_0 +2.2, p = 0.13 came from a rerun whose join kept the backup sections (n = 225) |
| D7 mutual adjustment | Sentence with the actual numbers | `fig4h_suppfig4cd_neighborhood_survival` | CT_3 + CT_0: −3.0 mo, p = 0.036; CT_0 + CT_3: +2.85 mo, p = 0.052 (τ = 3 y) |
| D8 missing covariates | Each panel keeps its coding; Methods say which | all survival scripts | Each script prints how many fitted patients the coding affects |
| D9 Fig 3c | Regenerated KM; MGMT + age + sex + p53 instead of site; print n | `fig3c_suppfig3b_wes_profile_os_rmst.ipynb` | 253 joined rows (252 patients); D9 model n = 173 (36 / 137; 74 lack p53, 6 MGMT): cluster 2 +5.1 mo [1.5, 8.6], p = 0.005. Published model (age, MGMT, sex, site as one ordinal covariate) reproduced: n = 247, +3.35 mo [0.3, 6.4], p = 0.0326. KM bands now render in the PDF |
| D10 OS range / median | Survival cohort; Schoenfeld as published and on the survival cohort | `fig1b_suppfig1_GBM_MS_Summary.ipynb` | n = 260: 1.2–121.7 mo, KM median 14.9. Schoenfeld as published (N = 287): MGMT 0.2031, age 0.0147, sex 0.8606; survival cohort: 0.384, 0.059, 0.766 |
| D11 Charité 13 | Text: immature follow-up | text only | |
| D12 age ≥ 45 | Give the range | `fig1b_suppfig1_GBM_MS_Summary.ipynb` | 37–89 y; CH_G_152 (37, not flagged) and CH_G_456 (41, `violates_cohort_inclusion_criteria_cohort1` = True, `cohort_violations` = "all_sample_sources_ok"). To exclude: add to `PATIENT_EXCLUSION` in definitions.py |
| D13 clinical variables | Count from the clinical table | `fig1b_suppfig1_GBM_MS_Summary.ipynb` | 251 columns (Supp Table 1c lists 255: the same plus 4 retired os/pfs-days columns; 4 are ID columns) |
| D14 Fig 5e | Keep the 10-gene panel | `fig5de_suppfig5cde_ECM_supplement.ipynb` | Reproduced: n = 192, 177 / 192, p = 1.66e-30 |
| D15 Supp Fig 5a | TIL scatter of Updated_figures | `fig5c_suppfig5ab_TIL_supplement.ipynb` | 57 TIL-high / 171 TIL-low, the 23 without TIL⁺ spots at the axis floor |
| D16 Supp Table 4a genes | One row per gene with the analyses it entered | `fig5de_suppfig5cde_ECM_supplement.ipynb` | CTSD, IGFBP7, LGALS3: commented out everywhere, used by no analysis. IGFBP2: only in the NicheNet pathway filter (NicheNet_v3 `ecm_genes_brain`) |
| D17 versions | Every script prints its package versions | all scripts | |
| D18 code | Supplement notebooks ported into these scripts | this folder | Spencer's source cells: `scratch/review_2026-11/source_cells_spencer_2026-11.ipynb`; review report: `scratch/review_2026-11/report/index.html` |

Also found while porting:
- **TIL⁻ zone, Fig 5c / 5d / Supp Fig 5c–e** (D3): with the 23 sections without TIL⁺ spots
  left out (148 TIL-low), TIL+ / BZ1 / BZ2 results are identical to the published tables;
  TIL⁻ changes (in Fig 5c, T_cell_IFN_response changes significance; Supp Fig 5c calls are
  unchanged; the BH correction across zones moves Supp Fig 5e FDRs by ≤ 0.03).
- **Fig 5d filter**: the Final_Notebooks tables behind Fig 5d kept section-zones with ≥ 20
  spots, Spencer's cell keeps all (`MIN_ZONE_SPOTS = 0`); stars differ (TIL+ ECM_0: * vs
  ***). To confirm which the published panel used.
- **Fig 6a / Supp Fig 6 text vs plotted axes**: BZ1 RGS2 (axis RGS1); TIL+ SRPX2 (SRPX),
  FGF2→SNAI1 (HGF→SNAI2), TGFβ1 / IL-15 from "interferon-responsive tumor cells" (from
  Macrophage_IFN_response), PRELP from TGFβ1/2 (from APOD); BZ2 COL1A2 and LUM come from
  CCN2 and BMP2, not TGFB3. Supp Fig 6a cannot be rebuilt from the notebook state (6 axis
  IDs missing from the current catalog); the frozen axes_selected.csv files are the input.
- **Fig 3a**: percentages are any alteration of the gene over 251 blocks (250 patients;
  ER_G_290a and b both kept, as it ran). "EGFR amplification ~52%" includes EGFR mutations
  (amplification alone 50.6%); "PTEN LOH" is PTEN loss of function (deletion or truncating).
- **Fig 3d**: the published panel is a Kruskal–Wallis test across the four clusters with BH
  FDR (identical to figure3/fig3e.csv), not the "T-test vs total cohort" of the legend.
  Microglia_homeostatic, called depleted in the text, has FDR 0.071.
- **Supp Fig 3c and CNV scores (0.196 etc.)**: computed with ER_G_673a in cluster 2, as a
  re-run of the clustering gives; the frozen cluster table (Fig 3b–d, Supp Fig 3b, 3e–f)
  has it in cluster 0. `fig3b_suppfig3ac_infercnv_wes_analyses.py` prints both.
- **PERMANOVA in the Methods** (pseudo-F 0.58 / 4.56): ran 999 permutations, Methods say
  9,999; needs scikit-bio, not ported.
- **Fig 5d filter (decision for Spencer)**: `fig5de_suppfig5cde_ECM_supplement.ipynb` keeps his
  `MIN_ZONE_SPOTS = 0` (n_low 148) and prints the ≥ 20-spot variant, which reproduces
  Final_Notebooks/focus_niches_significance.csv exactly (n_low 66 / 115 / 129; TIL+ ECM_0 *
  instead of ***). The published stars therefore came from the ≥ 20-spot filter.
- **Cell-type PERMANOVA (Methods "Level 4", pseudo-F 4.56)**: 4.56 is reproduced by level-3
  curated proportions (level 4 gives 4.02); `fig3d_correlation_cnv_wes_ecm.py` prints both.
  ECM PERMANOVA 0.58, p 0.87 reproduces. 999 permutations ran (Methods: 9,999).
- **Table 1a** (`table1a_sample_number_per_modality`): the as-published values are recomputed next to the corrected ones;
  corrected: WES 263 (250 with follow-up), Visium + snRNA + WES 213 (204).
- **Fig 2d–e** (`fig2de_celltype_modules`): reproduce exactly from NicheNet_v3 cells 162/172/174 (298 samples as
  loaded there, not 266; 25 modules; τ = 1.5 y chosen by the grid search). The forest shows
  modules 9, 8, 4, 2, 7, 24, 16; the text's list (7, 9, 11, 16, 24) includes module 11, which
  is not in the forest. † at 3 years: 9, 8, 4, 2 (text: 2, 4, 8). As-ran coding: 6 MGMT and
  75 p53 unknowns counted as the reference level; 20 fitted patients outside the survival
  cohort. Corrected (survival cohort, complete cases): only module 8 stays significant.
- **Fig 2b–c dot plots** (`fig2bc_suppfig2_marker_dotplots`): same genes and cell types as the published CSVs, but
  fraction expressing differs by up to 0.3 (e.g. GFAP in Malignant_AC1 0.959 vs 0.955). Not
  explained by the loader exclusions or data changes; most likely an older marker workbook
  (changed 31 Mar 2026; published files 21 Apr). Needs the author.
- **Fig 2a** (`fig2a_suppfig2_figure2`): same nuclei (1,123,420 from 265 samples), different UMAP layout (no
  published coordinates to compare). As committed, figure2.py could not produce the
  lineage UMAPs (scaled 2,000-HVG object reassigned); fixed. All 12 Supp Fig 2 stacked bars
  reproduce exactly. `.X` holds log1p values that are normalised and logged again (kept as
  ran; `UMAP_INPUT = "unnormalized_counts"` switches).
- **Fig 4** (`fig4*` files): Fig 4d is the "top 40 most variable genes" heatmap (legend: "top 5
  per neighborhood"); same genes and order as published. Fig 4e is `dotplot_composite_score`
  with its defaults (the cell last ran with other settings). Fig 4f is the residual heatmap
  of `fig4f_pathologist_annotation_analyses.py` (identical data); Supp Fig 4e the per-patient
  bars of `fig4e_suppfig4e_colocalization_plots` (identical). Supp Tables 3b–c reproduce. Global NMI 0.140, per-sample median
  0.237. Supp Table 3a lists 123 genes; the clustering panel is 117 (COL16A1, CTSD, IGFBP2,
  IGFBP7, LGALS3, PLOD3 commented out). ECM cell 46 also tested `imagerow` / `imagecol` as
  genes (fixed; no signature gene changes).
- **Fig 4a does not reproduce**: the published panel shows CT_0 and 36 cell types with scANVI
  names; no saved code produces it. The closest settings give 38 cell types (+ Mural_LEPR,
  Plasma_nos) and another row order. Needs the author.
- **Supp Table numbering**: the Fig 5f GSEA table is EFS supplementary_tables/5b_... but
  is cited as Supp Table 4b.
- **Fig 1b text mixes two cohorts**: MGMT +7.1, age −5.4, female +3.1 mo come from the
  GBM_MS_Summary model on every subject of the clinical tables (N = 287); "p53 p = 0.085"
  is the same model on the 260-patient survival cohort (p53 on N = 287: +3.7 mo, p = 0.026).
  `fig1b_suppfig1_GBM_MS_Summary.ipynb` prints both.
