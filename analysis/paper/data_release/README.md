# MOSAIC GBM data release: Watson et al. (2027)

Derived, per-patient data for the manuscript *"Multimodal Spatial Atlas of Glioblastoma
Reveals Conserved Perivascular ECM Border Zones Associated with Immune Infiltration"*.
The scripts in this folder rebuild every released file from the MOSAIC GBM source data
and check the manuscript's numbers against them.

**Release location (shared SageMaker EFS):**
`/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/exports/mosaic_gbm_watson2027/`
(`_work/` holds caches, logs and the full-cohort build, and is not part of the release).

**Released patients: deceased only.** Only patients with a recorded death are shared: 241
of the paper's 273 patients, plus 3 of the 5 Xenium validation patients. The scripts first
build every file for the full cohort (`_work/full/`, used to check the manuscript's
numbers), then `07_release_filter.py` writes the deceased subset here. Numbers in the
released files therefore differ slightly from the manuscript, which analysed all patients.

## Source data

MOSAIC GBM datasets referenced in `gbmcollaboration/definitions.py`: one dataset per centre
(CHUV, UKER, Charité) and modality, behind S3 access point
`pa-bi4jtrdsa151p4d3m-…-s3alias` (AWS account 767397995824). Processing versions:
snRNA v5, Visium v7, clinical v3, WES v2. Xenium: shared EFS `data/xenium/processed`.

## Cohort

The paper analysed 273 patients with IDH-wildtype GBM (CHUV 65, UKER 135, Charité 73),
one primary tumour block each. Loader exclusions: IDH-mutant, predefined patient
exclusions, recurrence samples. Recurrence samples are not part of any analysis in the
paper and are not released.
Patient ER_G_290 has two baseline blocks. Block `a` (low tumour content) is dropped, so
`ER_G_290b` is used. Both of its snRNA blocks are kept in `snrna/` (266 samples), so
`ER_G_290a` appears only there.

| Modality | Paper cohort | Released (deceased) | Notes |
|---|---|---|---|
| Clinical | 273 | 241 | Paper: 260 with mature follow-up. Released: 236 of those, plus 5 Charité patients with immature follow-up (OS NA) |
| snRNA-seq | 265 (266 samples) | 236 (237 samples) | One patient profiled on two blocks |
| Visium | 228 | 202 | QC-passing sections; 21 of the 228 are backup sections (`_vbu`, `_vb2`, `_bu2`, `_bu3`) |
| WES mutations | 263 | 231 | |
| WES CNV | 271 (264 clustered) | 239 (235 clustered) | |
| Xenium | 5 (6 sections per panel) | 3 (4 + 3 sections) | Validation patients outside the 273 |
| Visium + snRNA + WES mutations | 213 | 190 | |

`cohort/patient_manifest.csv` lists, for each patient, which modalities exist and which
block/section was used.

## Identifiers and pseudonymisation

- `sample_id` is the tumour block (e.g. `CH_G_151a`). `patient_id` is the same without
  the block letter (`CH_G_151`). These coded IDs are the only identifiers in the release.
- The prefix encodes the centre. Released centre codes are `CH-1` (CHUV), `CR-2`
  (Charité) and `ER-3` (UKER).
- IDs also appear **inside composite strings**: h5ad `obs_names` (`<sample_id>__<region>`),
  Visium `visium_section` (`<sample_id>_vis` / `_vbu`…), and Xenium `tile_id` / `section`
  (`<sample_id>_x1_tile_<x>_<y>`). `inventory.csv` lists every field that carries an ID.
  Apply one mapping across all files and fields.
- Not released: dates, free text, sex, age (exact or binned), sequencing reads, VCFs, germline
  or HLA data. Sex and age are withheld for privacy reasons: combined with centre and
  survival they would raise re-identification risk. Analyses that adjust for them (e.g. RMST
  models with age/sex covariates) cannot be fully re-run from the release alone.

## Files

### clinical/clinical_binned.csv
One row per released patient (241).

| Column | Values |
|---|---|
| `os_months` | Overall survival in months (1 decimal) |
| `os_event` | 1 for every released patient (deceased only) |
| `survival_eligible` | Whether the patient is in the paper's survival analyses. False for 5 Charité patients with immature follow-up, whose `os_months` is NA |
| `mgmt` | Methylated / Unmethylated |
| `p53` | Positive / Negative (IHC); NA when not recorded |
| `centre` | CH-1 / CR-2 / ER-3 |
| `til_group` | See below |

`til_group` is TIL-high or TIL-low, NA without Visium. TIL+ spot = summed cell2location
T/NK abundance ≥ 2. TIL-high = TIL+ spot fraction ≥ 75th percentile across the paper's 228
sections (cut-off 0.0304); TIL-low = the other sections, including the 23 with no TIL+
spot (paper: 57 / 171; released: 51 / 151). The paper's comparisons between TIL zones
leave the sections without TIL+ spots out, since border zones are defined around TIL+
spots: in the released files these are the TIL-low sections with no `TIL+` region.

### snrna/
- `pseudobulk_counts_sample_x_celltype.h5ad`
  - X: SoupX-corrected counts summed per sample × level-4 cell type (72 types), rounded
    to integers. Doublet-labelled nuclei are excluded.
  - obs: `n_cells`, `n_genes_detected`, cell-type levels 1–4.
  - var: `is_hvg`, `is_named_gene`, `named_gene_sources`.
- `n_cells_sample_x_celltype.csv`
- `n_genes_detected_sample_x_celltype.csv`: genes with ≥ 1 count over the full transcriptome.

### visium/
- `pseudobulk_counts_sample_x_ecm_neighborhood.h5ad` (ECM_0–ECM_9)
- `pseudobulk_counts_sample_x_til_zone.h5ad` (TIL+, BZ1, BZ2, TIL-)
  - X: raw Space Ranger counts summed over spots.
  - obs: `n_spots`, `til_group`, `visium_section`.
- `n_spots_sample_x_{ecm_neighborhood,til_zone}.csv`
- `c2l_proportions_sample_x_{ecm_neighborhood,til_zone}.csv`: cell2location level-4
  abundances summed over the region's spots and normalised to sum to 1.
  Long format: `sample_id`, `region`, `n_spots`, then one column per cell type.
  `Fibro_Cycling`, `T_other` and `cDC_CD207` are absent from some centres' cell2location
  references and are 0 there.

Spot-level region labels are the paper's frozen annotations (Supp Table 4a).

### wes/
- `cnv_log2ratio_sample_x_bin.csv`: mean CNVkit log2 copy ratio per autosomal genomic bin
  (1,295 bins).
- `cnv_bins.csv`: chromosome, start, end of each bin (GRCh38). Bins are windows of 200
  genes with a 10-gene step.
- `cnv_clusters.csv`: WES CNV Leiden cluster (0–3).
- `oncogenic_alterations_sample_x_gene.csv`: 0/1 per sample. Columns ending `_mut` are
  oncogenic SNV/indels (including TERT promoter C228T/C250T), `_amp` are amplifications,
  and `_lof` are loss of function (deletion or truncating variant).

### xenium/ (per panel: `hBrain_v1.2`, `hAtlas_v1.1`)
- `tile_counts_<panel>.h5ad`: SPLIT-corrected transcript counts summed per 55 µm tile.
  - obs: tile grid coordinates, `n_cells`, `til_count` (T cells + NK), and
    `in_paper_analysis` (≥ 5 cells).
- `tile_celltype_counts_<panel>.csv`: RCTD level-4 cell-type counts per tile.

Segmentation is `10x_mm_5um`, which is what produced Fig 6.

### genes/named_genes.csv
Genes named in the manuscript text or figures, with their source:
- Fig 2 markers (Supp 2a)
- matrisome panel (Supp 3a)
- NicheNet circos ligands, receptors and targets (Fig 6a, Supp 6)
- manuscript text

The pseudobulks release the union of these genes with 3,000 highly variable genes
(`seurat_v3`), computed separately per modality.

### Other
- `cohort/patient_manifest.csv`
- `inventory.csv`: every file, its dimensions and the fields that carry IDs.
- `REPORT.md`: handover summary of the release and of manuscript-vs-data findings.
- `sanity_report.md` / `sanity_checks.csv`: every check, with the values compared.

## Rebuilding

Run from this folder (`analysis/paper/data_release/`) with the `gbm_main` conda environment and S3 access, in this order
(01-06 and 00 write the full-cohort build to `_work/full/`):

```
python 01_named_genes.py
python 02_clinical.py
python 03_wes.py
python 04_singlecell_pseudobulk.py   # ~1.4M nuclei; needs ~250 GB RAM
python 05_visium_pseudobulk.py       # ~870k spots
python 06_xenium_tiles.py
python 00_cohort_manifest.py         # reads the outputs above
python 07_release_filter.py          # deceased subset -> release folder
python 98_inventory.py
python 99_sanity_checks.py
```
