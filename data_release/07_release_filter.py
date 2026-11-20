"""Write the released subset of the full-cohort build.

Only deceased patients may be shared: any patient with a recorded death in the clinical
tables. That is 241 of the 273 cohort patients (including 5 Charité patients that the
paper's survival analyses exclude for immature follow-up) and 3 of the 5 Xenium
validation patients (CH_G_532 is alive, CH_G_534 has no recorded status). For those 5, os_months stays NA (as in the paper) and
os_event is 1.

Every file of the build (_work/full) is subset to rows of released patients and written
to the release folder with the same relative path. Files without patient rows (gene list,
genomic bins) are copied. til_group in every file comes from the clinical table, so it
is defined once.
"""

import shutil

import anndata as ad
import numpy as np
import pandas as pd

import config as C
from common import out, patient_of, release_out

ID_RE = r"^[A-Z]{2}_G_\d+"


def released_patients() -> pd.Index:
    """Patients with a recorded death in the clinical tables.

    Read without the patient exclusions, because the Xenium validation patients
    (CH_G_532-536) are in PATIENT_EXCLUSION (pending final submission) yet have clinical
    rows. Only patients with rows in the build can end up in the release.
    """
    from gbmcollaboration import AllDataCenters

    clin = AllDataCenters.load_clinical(exclude_idh1_mut=False)
    if not C.RELEASE_DECEASED_ONLY:
        return pd.Index(patient_of(clin.index).unique())
    return pd.Index(clin.loc[clin["os_event"] == 1, "Subject Id"].astype(str).unique())


def filter_csv(rel, keep, til):
    """Every build CSV has its row ID in the first column: a sample, tile or patient ID,
    or a gene/bin for files without patient rows."""
    df = pd.read_csv(out(rel), index_col=0)
    ids = df.index.astype(str)
    if not ids.str.match(ID_RE).any():
        shutil.copy(out(rel), release_out(rel))
        return len(df), len(df), "copied (no patient rows)"
    pat = ids.to_series(index=df.index) if df.index.name == "patient_id" else pd.Series(
        patient_of(ids.str.extract(r"^([A-Z]{2}_G_\d+[a-z])", expand=False)), index=df.index)
    sub = df[pat.isin(keep).to_numpy()].copy()
    if "til_group" in sub.columns:
        sub["til_group"] = pat[pat.isin(keep)].map(til).to_numpy()
    if rel.endswith("clinical_binned.csv"):
        sub["os_event"] = 1.0
    sub.to_csv(release_out(rel))
    return len(df), len(sub), ""


def filter_h5ad(rel, keep, til):
    a = ad.read_h5ad(out(rel))
    pat = pd.Series(patient_of(a.obs["sample_id"].astype(str)), index=a.obs_names)
    b = a[pat.isin(keep).to_numpy()].copy()
    if "til_group" in b.obs.columns:
        b.obs["til_group"] = pat[b.obs_names].map(til).to_numpy()
    b.write_h5ad(release_out(rel), compression="gzip")
    return a.n_obs, b.n_obs, ""


def main():
    keep = released_patients()
    clin = pd.read_csv(out("clinical", "clinical_binned.csv"), index_col="sample_id")
    til = clin.set_index("patient_id")["til_group"]
    print(f"released cohort patients: {clin.patient_id.isin(keep).sum()} of {clin.patient_id.nunique()} "
          f"(deceased only: {C.RELEASE_DECEASED_ONLY})")
    rows = []
    for f in sorted(C.BUILD_DIR.rglob("*")):
        if f.is_dir():
            continue
        rel = str(f.relative_to(C.BUILD_DIR))
        release_out(rel)  # creates the folder
        if f.suffix == ".csv":
            n_in, n_out, note = filter_csv(rel, keep, til)
        elif f.suffix == ".h5ad":
            n_in, n_out, note = filter_h5ad(rel, keep, til)
        else:
            continue
        rows.append(dict(file=rel, rows_full=n_in, rows_released=n_out, note=note))
        print(f"  {rel}: {n_in} -> {n_out} {note}")
    pd.DataFrame(rows).to_csv(C.WORK_DIR / "release_filter_summary.csv", index=False)


if __name__ == "__main__":
    main()
