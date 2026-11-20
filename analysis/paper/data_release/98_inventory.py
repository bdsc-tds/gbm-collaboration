"""Inventory of released files: shape, size, and which fields carry patient/sample IDs
(what engineering must remap when pseudonymising)."""

import anndata as ad
import pandas as pd

import config as C
from common import release_out as out

ID_PATTERN = r"[A-Z]{2}_G_\d+"


def id_fields_df(df: pd.DataFrame, index_label: str):
    fields = []
    if df.index.astype(str).str.contains(ID_PATTERN).any():
        fields.append(f"{index_label} (index)")
    for c in df.columns:
        if not pd.api.types.is_numeric_dtype(df[c]) and df[c].astype(str).str.contains(ID_PATTERN).any():
            fields.append(c)
    return fields


def main():
    rows = []
    for f in sorted(C.OUT_DIR.rglob("*")):
        if f.is_dir() or "_work" in f.parts or f.name.startswith("sanity") or f.name in ("inventory.csv", "README.md", "REPORT.md"):
            continue
        rel = f.relative_to(C.OUT_DIR)
        size_mb = round(f.stat().st_size / 1e6, 1)
        if f.suffix == ".h5ad":
            a = ad.read_h5ad(f, backed="r")
            ids = [f"obs:{x}" for x in id_fields_df(a.obs, "obs_names")]
            rows.append(dict(file=str(rel), format="h5ad", rows=a.n_obs, cols=a.n_vars, size_mb=size_mb, id_fields=";".join(ids)))
        elif f.suffix == ".csv":
            df = pd.read_csv(f, index_col=0, low_memory=False)
            rows.append(dict(file=str(rel), format="csv", rows=len(df), cols=df.shape[1] + 1, size_mb=size_mb,
                             id_fields=";".join(id_fields_df(df, df.index.name or "first column"))))
    inv = pd.DataFrame(rows)
    inv.to_csv(out("inventory.csv"), index=False)
    print(inv.to_markdown(index=False))


if __name__ == "__main__":
    main()
