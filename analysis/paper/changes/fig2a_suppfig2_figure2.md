# fig2a_suppfig2_figure2.py

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
