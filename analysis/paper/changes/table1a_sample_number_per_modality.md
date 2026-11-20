# table1a_sample_number_per_modality.ipynb

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
