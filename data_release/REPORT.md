# Watson et al. (2027) data release: handover report

This report covers what was exported for the MOSAIC GBM manuscript and how its numbers
were checked against the data. File formats are described in `README.md`. Every check is
listed, with the values it compared, in `sanity_report.md` / `sanity_checks.csv`.

## 1. What was released

Only deceased patients are released: any patient with a recorded death, 241 of the
paper's 273 (236 of the 260 in the survival analyses, plus 5 Charité patients whose
follow-up was too short for the survival analyses; their OS is NA). Recurrence samples are
not part of any analysis in the paper and are not released. Patient counts below are
paper cohort → released.

| Modality | Content | Patients |
|---|---|---|
| Clinical | OS (months), event, MGMT, p53, centre, TIL group | 273 → 241 |
| snRNA-seq | Pseudobulk counts per sample × 72 level-4 cell types, plus n_cells and n_genes | 265 → 236 (237 samples) |
| Visium | Pseudobulk counts, n_spots and cell2location proportions per sample × ECM neighbourhood (10) and × TIL zone (4) | 228 → 202 |
| WES | CNV log2 ratio per sample × 1,295 autosomal bins; CNV cluster; oncogenic alteration table (32 genes) | CNV 271 → 239; clustered 264 → 235; mutations 263 → 231 |
| Xenium | Counts and RCTD cell-type counts per 55 µm tile, two panels | 5 → 3 (validation patients outside the 273) |

Analyses re-run on the released subset will give numbers that differ slightly from the
manuscript, which analysed every patient; the conclusions should hold.

Gene sets: 3,000 seurat_v3 highly variable genes per modality, plus the 523 genes named in
the manuscript, its figures and the supplementary tables (`genes/named_genes.csv`).

**Pseudonymisation.** The only identifiers are coded sample and patient IDs. They also
occur inside composite strings. `inventory.csv` lists every field that carries an ID.
Not released: dates, free text, sex, age (exact or binned), reads, VCFs, germline or HLA data.
**Sex and age are deliberately withheld for privacy reasons**: together with centre and
survival they would raise re-identification risk. Consequently, covariate-adjusted models
that use age or sex (RMST, Cox) cannot be fully reproduced from the release alone.

## 2. How the data were checked

- Every export was built from the same loaders and frozen annotation files as the
  paper's notebooks.
- Visium regions use the frozen spot annotations of Supp Table 4a.
- The Xenium tiles reproduce the per-section tile counts of the Fig 6 source data exactly
  (12/12 sections).
- In the TIL-zone export, the median T/NK proportion decreases from TIL+ (13%) through
  BZ1 (7%) and BZ2 (6%) to TIL− (5%), as expected.

## 3. Manuscript numbers vs data

Status meanings:
- **OK:** the paper's number matches the data.
- **Precise:** the number is right, but the text should say what it counts.
- **Inconsistent:** the number is wrong; the cause is given where it was found.
- **Untraced:** the number could not be reproduced.

### OK
- 273 patients; 260 with complete follow-up.
- 228 Visium sections.
- 72 level-4 cell types.
- Collagen higher in TIL-associated than TIL− spots in 177 of 192 patients. Here 192 =
  patients with ≥ 20 spots in each compartment.
- Xenium tile counts in Fig 6b/6g.

### Precise: number traceable, definition should be stated

| Claim | Paper | Data | What it actually is |
|---|---|---|---|
| Patients per centre (CHUV/UKER/Charité) | 72 / 152 / 73 (sum 297) | 65 / 135 / 73 (sum 273) | Clinical rows **before** relapse samples were excluded |
| snRNA-seq passing QC | 266 | 266 samples, 265 patients | Patient ER_G_290 contributed two blocks |
| Xenium validation | "six samples" | 6 sections per panel, 5 patients | Two patients contribute two sections to a panel |
| Fig 3a alteration frequencies | 52/26/21/19/18/16 % | 50.4/25.6/20.8/19.2/18.0/16.4 % | Denominator = 250 patients (mutation calls ∩ follow-up). "PTEN LOH" = PTEN loss of function (deletion or truncating variant) |
| CNV inference on "~1.2M malignant cells" | ~1.2M | 1,158,994 non-doublet nuclei, 747,131 malignant | ~1.2M is all nuclei (malignant + reference) |
| Fig 5b TIL survival group sizes | not stated | 55 vs 163 in the notebook | 218 survival-eligible sections; the 75th-percentile cut-off is **recomputed on these 218** (0.0355 instead of 0.0304), so CH_G_200a and CR_G_318a move from high to low. No-TIL+ sections counted as low. State n in the legend |

### Inconsistent: cause identified

**WES patient numbers.**
- The paper gives 261 (Results), 264 (Methods), and 274 / 261 (Supp 1a).
- In the committed `sample_number_per_modality.ipynb`, `load_wes()` loads the
  clinical table, not WES. So "274 WES" and "261 WES_surv" are clinical row counts, and
  the Results' 261 inherits this.
- Real numbers:
  - 263 patients with mutation calls.
  - 271 with CNV profiles.
  - 264 CNV-clustered; this matches the Methods' 264.

**Patients with all modalities (192).**
- Running the notebook's logic reproduces Supp 1a exactly, including 192. Two bugs
  combine:
  1. The WES loader bug above.
  2. 21 Visium backup sections (suffixes `_vbu`, `_vb2`, `_bu2`, `_bu3`) are not mapped
     back to their block, so those patients never match across modalities.
- Correct counts:
  - 213 patients with Visium, snRNA and WES mutations.
  - 204 when follow-up is also required.
- Supp Table 1a should be regenerated from `cohort/patient_manifest.csv`.

**Xenium segmentation.** The Methods state Proseg v3.0.10. The Fig 6 source data were
generated with the 10x multimodal segmentation (`10x_mm_5um`). The Methods should be
changed, or the figure regenerated.

**CNV genomic bins (Fig 3b).**
- `generate_cnv_wes_matrix.py` drops the wrong bins. The saved matrix has no chr21 bin
  and keeps chrX (which clustering removes later).
- So Fig 3b effectively excludes chr21. The release contains the 1,295 autosomal bins
  that exist.

**Cohort description.**
- "All patients aged ≥ 45": 2 patients are younger (youngest 37.1).
- OS range / median: the paper says 0.1–121 months / 15.12; the data give 1.2–121.7 /
  KM median 14.90. Specify whether the median is Kaplan–Meier or raw.
- "161 clinical features": Supp Table 1c has 255 rows.

**Text.**
- The matrisome gene list is called Supp Table 3a in the Results but Supp Table 5 in
  the Methods. Supp 5 is the Xenium panels, so 3a is correct.
- "Supp Table 4b–4c" for ECM neighbourhood genes and cell types should be 3b–3c; 4b is
  GSEA.
- ECM RMST covariates differ: MGMT, age, sex, p53 in the Results; age + MGMT in the
  Methods.
- Fig 3c covariates differ: age + MGMT in the legend; all covariates in the text.
- Typo: "most differentially associated with these zones in 47" should read
  "in TIL-low samples".

### Untraced: to ask the analysis author

**TIL-high / TIL-low / no TIL+ = 57 / 153 / 18: traced and confirmed by Spencer (D3).**
- 57 = TIL-high over all 228 sections (cut-off 0.0304); TIL-low = the other 171,
  including the 23 sections without TIL+ spots.
- 153 = the TIL-low arm of the Fig 5b model (207 patients: 54 vs 153, MGMT-adjusted
  ΔRMST +3.93 months, p = 0.0251).
- 18 = the 21 TIL-negative sections in that model minus the 3 without an MGMT result.
- Comparisons between TIL zones (Fig 5c/5d, Supp Fig 5c–e) leave the 23 sections without
  TIL+ spots out, because border zones are defined around TIL+ spots (57 vs 148).
- The released `til_group` is TIL-high / TIL-low for every Visium section (57 / 171 in the
  paper cohort, 51 / 151 released).

**n = 253 for the CNV-cluster RMST: reproduced.** 261 clinical rows (Charité exclusion
applied) ∩ 265 CNV-clustered samples = 253 rows. These are 252 patients (ER_G_290
contributes both blocks); 247 are modelled. Covariates: age, MGMT, sex, site. Cluster 2:
+3.35 months, p = 0.0326, as in the text. Earlier "252" counted ER_G_290 once.

**Survival models: all re-run and reproduced.** See RESOLUTION_PLAN.md §4:
- TIL (with an age-confounding caveat)
- ECM_3 (robust)
- CT_3 / CT_0 (the text mixes two models)
- the OS range and median (taken from a 287-patient cohort that includes IDH-mutant
  patients)

## 4. Release decisions to be aware of

- ER_G_290 has two baseline blocks. Block `a` (low tumour content) is excluded from the
  clinical table and the manifest, but both blocks are kept in the snRNA export, which
  has 266 samples.
- snRNA counts are SoupX-corrected, so not integer; they are rounded to integers after
  summing. Nuclei labelled as doublets are excluded.
- In the cell2location proportions, `Fibro_Cycling`, `T_other` and `cDC_CD207` are
  absent from some centres' reference files and are set to 0 there.
- The 5 deceased Charité patients excluded from the paper's survival analyses (immature
  follow-up at data lock) keep `survival_eligible` = False and OS = NA.
- Centre codes (`CH-1`, `CR-2`, `ER-3`) and the ID prefixes encode the centre. Codes will
  be replaced by the documented two-letter pseudonyms once agreed; the manuscript itself
  names the three centres.
