"""Check every count/claim in the manuscript against the released data.

For each claim: recompute from the exports (+ loaders where needed); if it differs,
try to reproduce the paper's number from a documented alternative (pre-exclusion
counts, modality intersections, notebook bugs) so the paper can be corrected with
the precise definition.

Writes sanity_report.md and sanity_checks.csv to the release folder.
"""

import json

import numpy as np
import pandas as pd

import config as C
from common import block_of, definitions_list, load_clinical, out, patient_of, release_out, til_labels

ROWS = []
EFS = C.EFS


def check(section, claim, paper, data, status, note=""):
    """status: OK | PRECISE (paper traceable, wording/definition to be made precise)
    | INCONSISTENT (paper contradicts data/itself) | UNTRACED"""
    ROWS.append(dict(section=section, claim=claim, paper=str(paper), data=str(data), status=status, note=note))
    print(f"[{status}] {section}: {claim} | paper={paper} | data={data}")


def main():
    man = pd.read_csv(out("cohort", "patient_manifest.csv"), index_col="patient_id")
    clin = pd.read_csv(out("clinical", "clinical_binned.csv"), index_col="sample_id")
    surv_excl = set(definitions_list("SURVIVAL_EXCLUSION"))
    til = til_labels()

    # ------------------------------------------------------------------ cohort
    n_pat = len(man)
    check("Cohort", "Patients in cohort", 273, n_pat, "OK" if n_pat == 273 else "INCONSISTENT")

    cc = man.centre.value_counts().to_dict()
    from gbmcollaboration import CHUVDataCenter, CRDataCenter, UKERDataCenter

    rows_pre_relapse = {}
    for name, dc in (("CH-1", CHUVDataCenter), ("ER-3", UKERDataCenter), ("CR-2", CRDataCenter)):
        d = dc.load_tabular(exclude_idh1_mut=True, exclude_patient_ids=True, exclude_relapse=False)
        rows_pre_relapse[name] = len(d["clinical"]["processed gbm clinical"])
    check(
        "Methods/Cohort", "Patients per centre CHUV/UKER/Charité", "72 / 152 / 73 (sum 297)",
        f"{cc['CH-1']} / {cc['ER-3']} / {cc['CR-2']} (sum {n_pat})", "PRECISE",
        "Paper numbers = clinical sample rows after IDH/patient exclusions but BEFORE relapse-sample "
        f"exclusion: {rows_pre_relapse['CH-1']} / {rows_pre_relapse['ER-3']} / {rows_pre_relapse['CR-2']}. "
        "Report patients after all exclusions.",
    )

    n_surv = int(man.survival_eligible.sum())
    check("Cohort", "Patients with complete follow-up (survival analyses)", 260, n_surv,
          "OK" if n_surv == 260 else "INCONSISTENT", "273 minus the 13 Charité patients with pending follow-up (SURVIVAL_EXCLUSION).")

    clin_full = load_clinical()
    age = clin_full.drop(index=C.DUPLICATE_BLOCK_DROP, errors="ignore")["age_at_diagnosis"]
    n_u45 = int((age < 45).sum())
    check("Methods/Cohort", "All patients aged >=45 at diagnosis", ">=45", f"min {age.min():.1f}; {n_u45} patients <45",
          "INCONSISTENT" if n_u45 else "OK", "Inclusion criterion statement contradicts data; correct or state exceptions.")

    os_m = clin.loc[clin.survival_eligible, "os_months"]
    from lifelines import KaplanMeierFitter

    km = KaplanMeierFitter().fit(os_m, clin.loc[clin.survival_eligible, "os_event"])
    check("Results/Fig 1b", "OS range / median (months)", "0.1-121 / 15.12",
          f"{os_m.min():.1f}-{os_m.max():.1f} / KM median {km.median_survival_time_:.2f}, raw median {os_m.median():.2f}",
          "OK" if abs(km.median_survival_time_ - 15.12) < 0.2 or abs(os_m.median() - 15.12) < 0.2 else "INCONSISTENT",
          "Specify KM median vs. raw median of follow-up.")

    n_feat = len(pd.read_csv(EFS / "supplementary_tables/1c_clinical_variables.csv", header=None)) if (EFS / "supplementary_tables/1c_clinical_variables.csv").exists() else None
    check("Results", "Clinical features compiled (Table 1c)", 161, f"{n_feat} rows in Supp Table 1c",
          "INCONSISTENT", "Either filter Table 1c to the 161 harmonised features or update the text.")

    # ------------------------------------------------------------------ modalities
    n_sn_s = int(man.snrna_block.dropna().str.count(";").add(1).sum())
    n_sn_p = int(man.has_snrna.sum())
    two = man.index[man.snrna_block.fillna("").str.contains(";")].tolist()
    check("Results/Methods", "snRNA-seq samples passing QC", 266, f"{n_sn_s} samples / {n_sn_p} patients",
          "PRECISE", f"Patient(s) with 2 snRNA blocks: {two}. Say '266 samples from {n_sn_p} patients'.")

    n_vis = int(man.has_visium.sum())
    check("Results/Methods", "Visium samples passing QC", 228, n_vis, "OK" if n_vis == 228 else "INCONSISTENT")

    n_mut, n_cnv, n_cl = int(man.has_wes_mutation.sum()), int(man.has_wes_cnv.sum()), int(man.has_wes_cnv_cluster.sum())
    n_mut_surv = int((man.has_wes_mutation & man.survival_eligible).sum())
    n_clin_rows = len(clin_full)
    check("Results 'WES data from 261 patients' / Methods 'WES passed QC for 264' / Supp 1a WES=274, WES_surv=261",
          "WES patients", "261 / 264 / 274 / 261",
          f"mutation calls {n_mut}; CNV profiles {n_cnv}; CNV-clustered {n_cl}; mutation∩survival {n_mut_surv}",
          "INCONSISTENT",
          f"Supp 1a 'WES' ({n_clin_rows}) and 'WES_surv' ({n_clin_rows - len(surv_excl & set(clin_full.index))}) are clinical-table "
          "row counts: load_wes() in the committed cohort notebook loads the clinical table, not WES. The Results' 261 "
          f"inherits this. 264 matches the number of CNV-clustered patients ({n_cl}). Use {n_mut} (oncoplot/mutations) "
          f"and {n_cl} (CNV clusters).")

    # Reproduce Supp 1a exactly with the notebook's logic (sample-level; '_vis' only stripped; WES:=clinical)
    vis_buggy = pd.Index(til.visium_section).str.replace("_vis", "")
    sn_blocks = pd.Index(sum(man.snrna_block.dropna().str.split(";").tolist(), []))
    clin_rows = pd.Index(clin_full.index)
    sets = {"Visium": set(vis_buggy), "snRNAseq": set(sn_blocks), "WES": set(clin_rows),
            "WES_surv": set(clin_rows) - surv_excl}
    supp1a = pd.read_csv(EFS / "supplementary_tables/1a_patient_count_per_modality.csv", index_col=0)["Count"] \
        if (EFS / "supplementary_tables/1a_patient_count_per_modality.csv").exists() else None
    repro = {k: len(v) for k, v in sets.items()}
    combos = {
        "Visium_snRNAseq": ["Visium", "snRNAseq"],
        "Visium_WES": ["Visium", "WES"],
        "Visium_WES_surv": ["Visium", "WES_surv"],
        "snRNAseq_WES": ["snRNAseq", "WES"],
        "snRNAseq_WES_surv": ["snRNAseq", "WES_surv"],
        "Visium_snRNAseq_WES": ["Visium", "snRNAseq", "WES"],
        "Visium_snRNAseq_WES_surv": ["Visium", "snRNAseq", "WES_surv"],
    }
    for combo, parts in combos.items():
        repro[combo] = len(set.intersection(*[sets[p] for p in parts]))
    match = supp1a is not None and all(supp1a.get(k) == v for k, v in repro.items())
    all_mod = int((man.has_visium & man.has_snrna & man.has_wes_mutation & man.survival_eligible).sum())
    all_mod_nosurv = int((man.has_visium & man.has_snrna & man.has_wes_mutation).sum())
    check("Results/Table 1a", "Patients with all modalities", 192,
          f"{all_mod} (Visium∩snRNA∩WES∩survival); {all_mod_nosurv} without survival requirement",
          "INCONSISTENT",
          f"Supp 1a reproduced exactly by notebook logic: {match} ({repro}). Two bugs: (1) 'WES' = clinical table; "
          "(2) Visium backup-section suffixes (_vbu/_vb2/_bu2/_bu3, 21 samples) not stripped, so those patients "
          "never intersect. Regenerate Table 1a from patient_manifest.csv.")
    pd.Series(repro).to_csv(C.WORK_DIR / "supp1a_reproduced_with_notebook_logic.csv")

    # ------------------------------------------------------------------ WES frequencies
    mut = pd.read_csv(out("wes", "oncogenic_alterations_sample_x_gene.csv"), index_col=0)
    onco = mut.loc[mut.index.isin(clin.index[clin.survival_eligible])]
    paper_freq = {"EGFR_amp": 52, "CDK6_amp": 26, "PTEN_lof": 21, "TERTp_mut": 19, "CDK4_amp": 18, "NTRK2_amp": 16}
    f_all, f_on = (mut.mean() * 100).round(1), (onco.mean() * 100).round(1)
    check("Results/Fig 3a", "Alteration frequencies EGFR amp/CDK6 amp/PTEN LOH/TERTp/CDK4 amp/NTRK2 amp (%)",
          "/".join(map(str, paper_freq.values())),
          f"oncoplot set n={len(onco)}: " + "/".join(str(f_on[k]) for k in paper_freq) +
          f"; all {len(mut)}: " + "/".join(str(f_all[k]) for k in paper_freq),
          "PRECISE", "Oncoplot = mutation calls ∩ clinical with follow-up (mutations_updated_oncoplot.py). "
          "State the denominator; 'PTEN LOH' is PTEN loss-of-function (deletion or truncating variant).")

    cl = man.dropna(subset=["cnv_cluster"])
    cands = {
        "cluster∩survival": cl.survival_eligible,
    }
    cfull = clin.reindex(cl.clinical_block)
    cfull.index = cl.index
    cands["+MGMT,age"] = cl.survival_eligible & cfull.mgmt.notna()
    cands["+MGMT,age,sex,p53"] = cands["+MGMT,age"] & cfull.p53.notna()
    cnts = {k: int(v.sum()) for k, v in cands.items()}
    check("Results/Fig 3c", "n for CNV-cluster RMST", 253, cnts, "PRECISE",
          "253 = CNV-clustered samples joined to clinical rows by block (252 patients; ER_G_290 has two blocks); "
          "247 enter the published model (age, MGMT, sex, site). Revised panel adjusts for MGMT, age, sex, p53 and "
          "prints the model n (Spencer, D9).")

    # ------------------------------------------------------------------ TIL
    g = til
    hi, lo_pos, zero = int((g.til_group == "TIL-high").sum()), int(((g.til_group == "TIL-low") & (g.n_til_pos > 0)).sum()), int((g.n_til_pos == 0).sum())
    check("Results/Fig 5", "TIL-high / TIL-low / no TIL+", "57 / 153 / 18", f"{hi} / {lo_pos + zero} ({zero} without TIL+ spots)",
          "PRECISE",
          f"57 = TIL-high over all 228 sections (cut-off {g.attrs['cutoff']:.4f}); TIL-low = the other {lo_pos + zero}, "
          f"incl. {zero} without TIL+ spots, which TIL-zone comparisons (Fig 5c/5d, Supp Fig 5c-e) leave out "
          f"({hi} vs {lo_pos}). 153 = TIL-low arm of the Fig 5b model (n = 207); 18 = the 21 TIL-negative sections "
          "of that model minus 3 without MGMT (Spencer, D3). Released til_group = 57 / 171.")
    s = man[man.has_visium & man.survival_eligible]
    sg = til.loc[block_of(s.visium_section)]
    check("Fig 5b", "TIL survival analysis group sizes", "54 / 153 (n = 207)",
          f"all sections mapped, with follow-up: high {int((sg.til_group=='TIL-high').sum())}, low {int((sg.til_group=='TIL-low').sum())} "
          f"(of which {int((sg.n_til_pos==0).sum())} without TIL+ spots)", "PRECISE",
          "Published model (TIL_Niche_Comparisons cell 32) drops the 21 backup sections at the clinical join and codes "
          "missing MGMT as unmethylated; kept as run (Spencer, D1), corrected n = 212 for the response letter.")

    col = EFS / "Final_Notebooks/collagen_engaged_vs_cold_stats.csv"
    if col.exists():
        cs = pd.read_csv(col)
        check("Results/Fig 5e", "Collagen higher in TIL-associated vs TIL- spots", "177 of 192 patients",
              cs.to_dict("records")[-1] if len(cs) else "?", "OK",
              "192 here = patients with >=20 spots in both compartments (not the 192 all-modality patients); say so.")

    # ------------------------------------------------------------------ snRNA
    pb_obs = pd.read_csv(out("snrna", "n_cells_sample_x_celltype.csv"), index_col=0)
    n_types = pb_obs.shape[1]
    c2l_types = pd.read_csv(out("visium", "c2l_proportions_sample_x_ecm_neighborhood.csv"), nrows=1).shape[1] - 3
    check("Results/Methods", "Level-4 cell types", 72, f"{n_types} non-doublet level-4 labels in snRNA; {c2l_types} in cell2location reference",
          "OK" if 72 in (n_types, c2l_types) else "PRECISE")
    import anndata as ad

    sobs = ad.read_h5ad(out("snrna", "pseudobulk_counts_sample_x_celltype.h5ad"), backed="r").obs
    n_mal = int(sobs.loc[sobs.celltype_level1_curated == "Malignant", "n_cells"].sum())
    n_tot = int(sobs.n_cells.sum())
    check("Results", "scRNA CNV inference across ~1.2M malignant cells from 266 tumours", "~1.2M malignant",
          f"{n_mal:,} malignant of {n_tot:,} non-doublet nuclei", "PRECISE" if abs(n_tot - 1.2e6) < 0.1e6 else ("OK" if n_mal > 1.1e6 else "INCONSISTENT"),
          "~1.2M matches all non-doublet nuclei (inferCNV input = malignant + reference), not malignant cells; say '~1.2M nuclei'.")

    # ------------------------------------------------------------------ Xenium
    for p in C.XENIUM_PANELS:
        f = out("xenium", f"tile_celltype_counts_{p}.csv")
        if f.exists():
            t = pd.read_csv(f, usecols=["section", "sample_id", "in_paper_analysis", "til_count"])
            mine = t[t.in_paper_analysis & (t.til_count > 0)].groupby("section").size()
            fig = pd.read_csv(EFS / f"figures/xenium/{p}/{C.XENIUM_SEGMENTATION}/ecm_spearman_tiles.csv")
            fig = fig.drop_duplicates("donor").set_index("donor")["n_tiles"]
            same = mine.reindex(fig.index).eq(fig).all()
            check("Fig 6b/6g", f"Tiles analysed per section ({p}; >=5 cells, TIL>0)", fig.to_dict(), mine.to_dict(),
                  "OK" if same else "INCONSISTENT", "Export reproduces the figure's tile set exactly." if same else "")
            check("Results/Methods Xenium", f"Validation samples ({p})", "six samples",
                  f"{t.section.nunique()} sections from {t.sample_id.nunique()} patients", "PRECISE",
                  "Six sections per panel from five patients (CH_G_534a/535a contribute two sections in one panel).")
    check("Methods Xenium", "Segmentation", "Proseg v3.0.10", f"Fig 6b/6g generated with {C.XENIUM_SEGMENTATION} (figures/xenium_fig_6.zip)",
          "INCONSISTENT", "Methods should state 10x multimodal segmentation (5 um expansion) for Fig 6, or regenerate with Proseg.")

    # ------------------------------------------------------------------ data-processing issues found
    check("WES CNV matrix", "Autosomal bins chr1-22 (Fig 3b)", "chr1-22", "chr21 bin absent (bin 1274); chrX kept in saved matrix, removed before clustering",
          "INCONSISTENT", "generate_cnv_wes_matrix.py drops the wrong bins; Fig 3b effectively excludes chr21.")
    for claim, paper, note in [
        ("Matrisome gene list location", "Supp Table 3a (Results) vs Supp Table 5 (Methods)", "Supp 5 is the Xenium panels; use 3a."),
        ("ECM neighbourhood genes / cell types", "'Supp Table 4b-4c'", "Files are 3b/3c; 4b is GSEA."),
        ("ECM RMST covariates", "Results: MGMT, age, sex, p53; Methods: age + MGMT", "Align Methods with the model actually run."),
        ("Fig 3c covariates", "legend: age + MGMT; text: all covariates", "Align."),
        ("Typo", "'most differentially associated with these zones in 47'", "Should read 'in TIL-low samples'."),
    ]:
        check("Text", claim, paper, "-", "INCONSISTENT", note)

    df = pd.DataFrame(ROWS)
    df.to_csv(release_out("sanity_checks.csv"), index=False)
    with open(release_out("sanity_report.md"), "w") as fh:
        fh.write("# Sanity checks: manuscript vs released data\n\n")
        fh.write("Status: OK = matches; PRECISE = paper number traceable, definition must be stated precisely; "
                 "INCONSISTENT = contradicts data or itself; UNTRACED = could not reproduce.\n\n")
        fh.write(df.to_markdown(index=False))
    print(df.status.value_counts())


if __name__ == "__main__":
    main()
