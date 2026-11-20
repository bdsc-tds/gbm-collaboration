"""Cohort manifest: one row per patient (273), modality availability and the block used.

Run after 02-05 (reads their outputs). Matching is done at patient level
(block letter dropped) after normalising Visium section suffixes
(_vis/_vbu/_vb2/_bu2/_bu3), so that a modality profiled on a different block of the
same patient is still counted (the block used is recorded per modality).
"""

import pandas as pd

import config as C
from common import block_of, out, patient_of, til_labels


def by_patient(blocks, name):
    blocks = pd.Index(blocks)
    s = pd.Series(blocks, index=patient_of(blocks), name=f"{name}_block")
    dup = s.index[s.index.duplicated()].unique()
    if len(dup):
        print(f"[{name}] patients with >1 block: {list(dup)} -> {s[s.index.isin(dup)].to_dict()}")
    return s.groupby(level=0).agg(lambda x: ";".join(sorted(x)))


def main():
    clin = pd.read_csv(out("clinical", "clinical_binned.csv"), index_col="sample_id")
    m = clin[["patient_id", "centre", "survival_eligible", "til_group"]].reset_index().set_index("patient_id")
    m = m.rename(columns={"sample_id": "clinical_block"})

    sn = pd.read_csv(out("snrna", "n_cells_sample_x_celltype.csv"), index_col=0)
    til = til_labels()
    cnv = pd.read_csv(out("wes", "cnv_log2ratio_sample_x_bin.csv"), index_col=0, usecols=[0])
    cl = pd.read_csv(out("wes", "cnv_clusters.csv"), index_col=0)
    mut = pd.read_csv(out("wes", "oncogenic_alterations_sample_x_gene.csv"), index_col=0, usecols=[0])

    for name, blocks in [
        ("snrna", sn.index), ("visium", til.index), ("wes_mutation", mut.index),
        ("wes_cnv", cnv.index), ("wes_cnv_cluster", cl.index),
    ]:
        m = m.join(by_patient(blocks, name))
        m[f"has_{name}"] = m[f"{name}_block"].notna()
    m["visium_section"] = til.set_index(patient_of(til.index))["visium_section"].reindex(m.index)
    m["cnv_cluster"] = cl.set_index(patient_of(cl.index))["cnv_cluster"].reindex(m.index)
    m["n_snrna_nuclei"] = sn.sum(axis=1).groupby(patient_of(sn.index)).sum().reindex(m.index)

    # samples present in a modality but whose patient is not in the clinical cohort
    for name, blocks in [("snrna", sn.index), ("visium", til.index), ("wes_mutation", mut.index), ("wes_cnv", cnv.index)]:
        extra = sorted(set(patient_of(blocks)) - set(m.index))
        if extra:
            print(f"[{name}] patients not in clinical cohort: {extra}")

    m.to_csv(out("cohort", "patient_manifest.csv"))
    print(m.filter(like="has_").sum().to_dict())


if __name__ == "__main__":
    main()
