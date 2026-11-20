# Analyses of Watson et al., MOSAIC GBM

Every panel, supplementary table and number of the manuscript comes from one file in this
folder. Files are named after the panels they make (`fig4h_…`, `suppfig3d_…`,
`table1a_…`, `supptables_4a_…`) and run top to bottom in the `gbm_main` environment.

Most files are the authors' notebooks and scripts, moved here with `git mv` and edited as
little as possible (regular loaders, exploratory cells removed, bugs fixed). Each file opens
with what it produces, where it comes from and how it compares with the published output.
What was changed in each file, and why, is in [`CHANGES.md`](CHANGES.md); `git log --follow -p
<file>` shows the diff against the original.

## Layout

```
analysis/
  paper/
    figNN_…, suppfigNN_…, tableNN_…   one file per group of panels
    figNN_0_…                         upstream steps that wrote frozen inputs (slow, not re-run)
    supptables_…                      supplementary tables not drawn by a figure file
    data_release/                     per-patient data release (see its README)
    CHANGES.md                        changes made to each file
  scratch/                            code not used in the paper, in its original folders
```

Upstream steps sit next to the figure they feed, with `_0_` in the name:
`fig3a_0_mutation_calling.py` (WES mutation matrix), `fig3b_0_infercnv_sc.py` (snRNA
inferCNV), `fig3b_0_cnv_wes_matrix.py` (WES CNV matrix), `fig5_0_TIL_zones.ipynb` (TIL zones),
`fig6a_0_NicheNet_runs.ipynb` (zone DE and NicheNet). The ECM, cell-type and WES CNV clusterings
are cells of their figure files, off by default (`RUN_CLUSTERING = False` reads the frozen labels);
`ecm_neighborhood_annotations_v4.csv` (Supp Table 4a) is written by the clustering cells of
`fig4cdg_suppfig4f_ECM_niches_v4.ipynb`.

## Data loading and shared code

Data come from the regular loaders, `gbmcollaboration.AllDataCenters` / `DataCenter`, with
the exclusion lists of `gbmcollaboration/definitions.py`. Four options were added to them
for the paper analyses:

| Loader | Use |
|---|---|
| `AllDataCenters.load_clinical(...)` | the processed clinical tables of the three centres, concatenated, with a `center` column |
| `load_singlecell(..., obs_only=True)` | cell metadata only (seconds instead of the full ~270 GB object) |
| `load_spatial_adata(genes=..., list_sample_id=...)` | only the listed genes, read in batches (the dense cohort is > 100 GB) |
| `load_cell2location(list_sample_id)` | cell2location abundances per spot (the per-section h5ad files on S3) |

`gbmcollaboration/paper/` holds only what several analyses share: `config` (paths of frozen
inputs, constants), `ids` (patient / block / section IDs, keeping backup sections `_vbu`,
`_vb2`, `_bu2`, `_bu3` with their block), `cohort` (the 273-patient cohort, 260 in survival
analyses), `annotation` (cell-type key), `visium`, `til` (TIL zones, TIL-high / TIL-low),
`neighborhoods`, `snrna` and `output`. A function used by one file is defined at the top of
that file. To remove a patient from every analysis, add them to `PATIENT_EXCLUSION` in
`definitions.py`.

## Outputs and source data

Each file writes to EFS `figures/paper/<file name>/` (`PAPER_OUT` to change). Every
manuscript panel is also written to **`figures/paper/panels/`** under its manuscript name
(`Fig4h.png`, `.pdf`, and `.csv` with the plotted data; `SuppTable3b.csv`), so running every
file fills that folder with the whole paper. UMAPs have no data CSV, and survival panels give
aggregate results (RMST, n, CI, p), never per-patient times.

## Files

| File | Panels and numbers | Origin |
|---|---|---|
| `table1a_sample_number_per_modality.ipynb` | Table 1a | analysis/cohort_desc/ |
| `table1bc_cohort_description.py` | Table 1b–1c | analysis/cohort_desc/cohort_description.py |
| `supptables_1bc_pathology_clinical_variables.ipynb` | Supp Tables 1b–1c (pathology labels, clinical variables) | annotation loading of colocalization_plots; a commented line of scratch/survival/survival_patient_features.ipynb |
| `fig1b_suppfig1_GBM_MS_Summary.ipynb` | Cohort numbers, Fig 1b, Supp Fig 1a–c | analysis/survival/ (Spencer) |
| `fig2a_suppfig2_figure2.py` | Fig 2a, Supp Fig 2 UMAPs and composition bars | figures/figure2.py |
| `fig2bc_suppfig2_marker_dotplots.ipynb` | Fig 2b–c, Supp Fig 2 dot plots | analysis/annotation/ |
| `fig2de_celltype_modules.ipynb` | Fig 2d–e | NicheNet_v3 module cells (Spencer) |
| `fig3a_0_mutation_calling.py` | WES mutation matrix (upstream) | analysis/mutations/mutations.py |
| `fig3a_oncoplot.py` | Fig 3a, alteration frequencies | analysis/mutations/ |
| `fig3b_0_infercnv_sc.py`, `fig3b_0_cnv_wes_matrix.py` | snRNA inferCNV, WES CNV matrix (upstream) | analysis/infercnv/ |
| `fig3b_suppfig3ac_infercnv_wes_analyses.py` | Fig 3b, Supp Fig 3a, 3c, CNV scores | analysis/infercnv/ |
| `fig3c_suppfig3b_wes_profile_os_rmst.ipynb` | Fig 3c, Supp Fig 3b | analysis/survival/ |
| `fig3d_correlation_cnv_wes_ecm.py` | Fig 3d, PERMANOVA | analysis/multimodal_correlation/ |
| `suppfig3d_wes_cluster_covariate_rmst.ipynb` | Supp Fig 3d | reconstruction (no code found) |
| `suppfig3ef_analyses_wes_sc_infercnv.py` | Supp Fig 3e–f | analysis/infercnv/ |
| `fig4ab_suppfig4abg_Cell_Neighborhoods.ipynb` | Fig 4a–b, Supp Fig 4a, 4b, 4g | analysis/survival/ (Spencer) |
| `fig4cdg_suppfig4f_ECM_niches_v4.ipynb` | Fig 4c, 4d, 4g, Supp Fig 4f, Supp Table 3a | analysis/survival/ (Spencer) |
| `fig4e_suppfig4e_colocalization_plots.ipynb` | Fig 4e, Supp Fig 4e, Supp Tables 3b–c | analysis/cell2loc/ |
| `fig4f_pathologist_annotation_analyses.py` | Fig 4f | analysis/hne/ |
| `fig4h_suppfig4cd_neighborhood_survival.ipynb` | Fig 4h and †, Supp Fig 4c | Updated_figures (Spencer) |
| `suppfig4d_CT0_CT3_composition.ipynb` | Supp Fig 4d | reconstruction (no code found) |
| `fig5_0_TIL_zones.ipynb` | TIL zones (upstream) = Supp Table 4a til_zone_annotations | TIL_Niche_Comparisons_v2_clean cells 20–24 (Spencer) |
| `fig5a_TIL_zones_map.ipynb` | Fig 5a | TIL_Niche_Comparisons_v2_clean facet cell (Spencer) |
| `fig5b_TIL_survival.ipynb` | Fig 5b, TIL vs age, sensitivity table | TIL_Niche_Comparisons RMST cell (Spencer) |
| `fig5b_TIL_spline.ipynb` | TIL burden as a continuous variable: spline RMST curve, spline tests, survRM2 cutoff scan (response to reviewers) | continuous_hr.py / continuous_rmst.ipynb (Jonathan Bac) |
| `fig5c_suppfig5ab_TIL_supplement.ipynb` | TIL groups, Supp Fig 5a–b, Fig 5c | TIL_supplement (Spencer) |
| `fig5de_suppfig5cde_ECM_supplement.ipynb` | Fig 5d–e, Supp Fig 5c–e, gene use of the Supp Table 3a genes | ECM_supplement (Spencer) |
| `fig5fg_TIL_zone_GSEA.ipynb` | Fig 5f–g | TIL_Niche_Comparisons_v2_clean (Spencer) |
| `supptables_4a_ecm_til_zone_tables.ipynb` | Supp Table 4a ecm_neighborhood_TIL_group_by_zone_pseudobulk, ecm_neighborhood_niche_matrisome_summary | TIL_Niche_Comparisons_v2_clean cells 70–72 (Spencer); matrisome summary rebuilt from its columns |
| `supptables_4a_til_zone_celltype_stats.ipynb` | Supp Table 4a til_zone_celltype_stats_long, til_zone_celltype_summary | reconstruction (no code found) |
| `supptables_4b_gsea_TIL_zones.ipynb` | Supp Table 4b (zone GSEA run behind Fig 5f–g) | TIL_Niche_Comparisons_v2_clean cells 60–61 (Spencer) |
| `fig6a_0_NicheNet_runs.ipynb` | zone DE, NicheNet runs, selected axes (upstream) | NicheNet_v3 (Spencer) |
| `fig6a_suppfig6_NicheNet_circos.ipynb` | Fig 6a, Supp Fig 6a–b | NicheNet_v3 circos cells (Spencer) |
| `fig6bg_xenium_niches.ipynb` | Fig 6b, 6g | analysis/xenium/niches.ipynb |

Not produced by code: Fig 1a, 1c (schematics), Fig 6c–f (microscopy), and the supplementary
tables that are gene lists or vendor panels (2a marker workbook, 3a matrisome genes, 5 Xenium
panels). Every other supplementary table regenerates and matches the published file (Supp
Table 4b: one of seven comparisons differs by < 1e-3 in NES / FDR, from ties in the ranking).

## Review answers (D1–D18, Spencer, November 2026)

| Item | Decision | File | Result |
|---|---|---|---|
| D1 Fig 5b cohort | Keep the published model; corrected cohort for the response letter | `fig5b_TIL_survival` | As ran: n = 207 (54 / 153), +3.93 mo, p = 0.0251 (21 backup sections lost at the join, 5 missing MGMT counted unmethylated). Corrected: n = 212 (55 / 157), +4.00 mo [0.6, 7.4], p = 0.022 |
| D2 TIL and age | Fig 5b stays MGMT-adjusted; text adds ρ and the 4-covariate estimate | `fig5b_TIL_survival` | ρ = −0.22, p = 0.001; MGMT + age + sex + p53: +2.90 mo [−0.3, 6.1], p = 0.078; sensitivity table incl. probe-batch rows |
| D3 57 / 153 / 18 | TIL-low = 171 incl. the 23 without TIL⁺ spots; zone comparisons 57 vs 148 | `til`, `fig5c…`, `fig5de…` | 18 = the 21 TIL-negative sections of the Fig 5b model minus 3 without MGMT |
| D4 Fig 4h source | Updated_figures sweep, as it ran | `fig4h…` | ECM_3 adj p = 0.0114, ECM_9 0.4743, ECM_1 0.2595 (published 0.2598). 10 fitted patients outside the survival cohort; 5 MGMT and 51 p53 values counted as the reference level |
| D5 ECM_3 / ECM_1 text | Numbers from the Fig 4h model | `fig4h…` | ECM_3 −3.5 mo [−6.2, −0.8], p = 0.011; † = ECM_3 only (−1.7 mo at 1.5 y, p = 0.013). ECM_1 unadjusted +2.7 mo, p = 0.083 in this model (text's +3.0, p = 0.042 is from an MGMT-only cell) |
| D6 CT_3 / CT_0 | One model | `fig4h…` | The sweep as it ran gives both published numbers: CT_3 −3.9 mo, p = 0.006; CT_0 +3.3 mo [0.4, 6.2], p = 0.026 |
| D7 mutual adjustment | Sentence with the numbers | `fig4h…` | CT_3 + CT_0: −3.0 mo, p = 0.036; CT_0 + CT_3: +2.85 mo, p = 0.052 (τ = 3 y) |
| D8 missing covariates | Each panel keeps its coding | all survival files | Each prints how many fitted patients the coding affects |
| D9 Fig 3c | MGMT + age + sex + p53 instead of site; print n | `fig3c…` | n = 173 (36 / 137): +5.1 mo [1.5, 8.6], p = 0.005. Published model reproduced: n = 247, +3.35 mo [0.3, 6.4], p = 0.0326 |
| D10 OS range / median | Survival cohort; Schoenfeld both ways | `fig1b…` | n = 260: 1.2–121.7 mo, KM median 14.9. Schoenfeld N = 287: MGMT 0.2031, age 0.0147, sex 0.8606; N = 260: 0.384, 0.059, 0.766 |
| D11 Charité 13 | Text: immature follow-up | text | |
| D12 age ≥ 45 | Give the range | `fig1b…` | 37–89 y; CH_G_152 (37) and CH_G_456 (41, flagged `violates_cohort_inclusion_criteria_cohort1`) |
| D13 clinical variables | Count from the clinical table | `fig1b…`, `supptables_1bc…` | 251 columns (Supp Table 1c: 255 incl. 4 retired os/pfs-days columns) |
| D14 Fig 5e | Keep the 10-gene panel | `fig5de…` | n = 192, 177 / 192, p = 1.66e-30 |
| D15 Supp Fig 5a | TIL scatter | `fig5c…` | 57 / 171; the 23 without TIL⁺ spots at the axis floor |
| D16 Supp Table 4a genes | Gene use table | `fig5de…` | CTSD, IGFBP7, LGALS3 used by no analysis; IGFBP2 only in the NicheNet pathway filter |
| D17 versions | Every file prints its versions | all | |
| D18 code | Supplements ported | this folder | Spencer's source cells and the review report: `analysis/scratch/review_2026-11/` |

## Findings for the authors

Reproduced, but text or legend should change:
- **Fig 1b**: +7.1 / −5.4 / +3.1 mo come from all 287 subjects of the clinical tables; "p53
  p = 0.085" is the same model on the 260-patient survival cohort (on 287: +3.7 mo, p = 0.026).
- **Fig 2e**: the forest shows modules 9, 8, 4, 2, 7, 24, 16; the text cites 7, 9, 11, 16, 24
  (module 11 is not in the forest). † at 3 y: 9, 8, 4, 2 (text: 2, 4, 8). The modules were
  built on 298 samples (no loader exclusions). Survival cohort, complete cases: only module 8
  stays significant.
- **Fig 3a**: "EGFR amplification ~52%" counts amplification or mutation (amplification alone
  50.6%); "PTEN LOH" is PTEN loss of function; denominator 251 blocks (250 patients).
- **Fig 3d**: Kruskal–Wallis across the four clusters with BH FDR, not the legend's t-test.
- **Supp Fig 3c and CNV scores**: computed with ER_G_673a in cluster 2 (as a re-clustering
  gives); the frozen cluster table has it in cluster 0. Both are printed.
- **PERMANOVA (Methods)**: the cell-type pseudo-F 4.56 uses level-3 proportions (level 4:
  4.02); 999 permutations ran, not 9,999.
- **Fig 4d**: top 40 most variable genes, not "top 5 per neighborhood".
- **Fig 4a**: closest reconstruction has 38 cell types; dropping Mural_LEPR and Plasma_nos gives
  the published 36 (likely a cell2location run without those factors, not on S3 / EFS).
- **Fig 5c–d, Supp Fig 5c–e**: TIL⁻ zone with the 23 no-TIL⁺ sections left out (D3);
  TIL+ / BZ1 / BZ2 identical to the published tables.
- **Fig 5d**: the published stars use section-zones with ≥ 20 spots; Spencer's code keeps all.
  Both are printed; to decide.
- **Fig 6 text vs plotted axes**: BZ1 RGS2 (axis RGS1); TIL+ SRPX2 (SRPX), FGF2→SNAI1
  (HGF→SNAI2), TGFβ1 / IL-15 from Macrophage_IFN_response, PRELP from APOD; BZ2 COL1A2 and LUM
  from CCN2 and BMP2. The TIL+ axis selection cannot be rebuilt from the current catalog.
- **TIL burden as a continuous variable** (`fig5b_TIL_spline`, 218 patients): the spline RMST
  curve rises steeply up to a TIL fraction of ~0.04 and then flattens; its inflection (0.037,
  75th percentile) coincides with the Fig 5b cutoff (0.0304). Non-linearity is not significant
  (likelihood-ratio p = 0.21 raw, 0.67 log scale; the module's Wald p = 0.45). survRM2 cutoff scan:
  ΔRMST is positive for every cutoff from the 50th to the 90th percentile, largest around the
  75th, which is the only cutoff significant after MGMT adjustment (+4.0 mo, p = 0.022; nominal,
  cutoffs not corrected for multiplicity).
- **Supp Fig 3d** (reconstruction): MGMT, age and sex unchanged with cluster 2 as covariate,
  but p53 moves from p = 0.082 to 0.033.
- **Supp Fig 4d** (reconstruction): CT_0 +7.1 pp lymphoid/DC reproduces; CT_3 −6.7 pp
  mesenchymal does not with any definition tried.
- **Table 1a**: published WES / all-modality counts came from two bugs (clinical table read as
  WES; backup Visium sections unmatched). Corrected: WES 263 (250 with follow-up),
  Visium + snRNA + WES 213 (204).
- **Supp Table 4a, ECM tables (decision needed)**: the published
  ecm_neighborhood_TIL_group_by_zone_pseudobulk and ecm_neighborhood_niche_matrisome_summary were
  built on an earlier ECM clustering (`ecm_neighborhood_annotations.csv`, not the v4
  neighborhoods of the paper; 13% of spots share a label) and include section CR_G_335a_vis
  (excluded since 12 Dec 2025; 58 / 171 groups). `supptables_4a_ecm_til_zone_tables.ipynb`
  reproduces them as published and also writes the paper-cohort version (v4 labels, 228
  sections, `*_paper_cohort_v4.csv`), which differs markedly (max |Δ log2FC| 3.1). Replace or
  footnote.
- **Supp Table 3a**: the 117-gene panel written by `fig4cdg…` differs from the published list in
  ordering (COL16A1) and one gene (TNR).
- **Supp Table numbering**: the text's numbering differs from the files (e.g. the Fig 5f GSEA
  table is cited as 4b).
