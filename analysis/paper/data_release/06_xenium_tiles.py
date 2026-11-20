"""Xenium 55 um tiles (Fig 6b/6g validation cohort).

Tiles of analysis/paper/fig6bg_xenium_niches.ipynb (Fig 6b, 6g) for both panels:
segmentation 10x_mm_5um (what produced figures/xenium_fig_6.zip), SPLIT
(split_fully_purified) corrected counts, RCTD (class-aware) level-4 labels.

Outputs (xenium/), per panel:
  tile_counts_<panel>.h5ad
      X   : SPLIT-corrected transcript counts summed over the cells of each tile
            (all panel genes; panels are targeted so no HVG filtering)
      obs : sample_id (patient block), section, panel, tile_x, tile_y (grid index),
            x_um, y_um (tile origin, um), n_cells, til_count (T_cell + NK, level 3),
            in_paper_analysis (n_cells >= 5, as in Fig 6b/6g)
  tile_celltype_counts_<panel>.csv
      tile_id, sample_id, section, n_cells, til_count, one column per level-4 cell type
"""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import config as C
from common import out

from gbmcollaboration.utils import xenium_io_utils as xio

XENIUM_LEVELS = ["segmentation", "condition", "panel", "donor", "sample"]
LEVEL_KEYS = [f"celltype_level{i}_curated_first_type" for i in range(1, 5)]
READ_KWARGS = dict(
    anndata=True, cells_as_circles=False, cells_boundaries=False, cells_boundaries_layers=False,
    nucleus_boundaries=False, cells_labels=False, nucleus_labels=False, transcripts=False,
    morphology_mip=False, morphology_focus=False, aligned_images=False, gex_only=True,
    pool_mode="thread", max_workers=6,
)


def load_panel(panel):
    paths, annot_paths = xio.discover_xenium_paths(
        base_dir=C.XENIUM_BASE, normalisation="lognorm", reference="allcenters_10k_ct",
        method="rctd_class_aware", level=None, load_correction=True,
        split_level=C.XENIUM_SPLIT_LEVEL, segmentations_filter=[C.XENIUM_SEGMENTATION],
        conditions_filter=["gbm"], panels_filter=[panel],
    )
    xds = xio.read_count_correction_samples(paths, [C.XENIUM_CORRECTION])
    xds["raw"] = xio.read_xenium_samples(paths["raw"], **READ_KWARGS)
    xio.read_annotations(xds, ["raw", C.XENIUM_CORRECTION], annot_paths, None, max_workers=6)
    for k, xd in xds[C.XENIUM_CORRECTION].items():
        xd.obsm["spatial"] = xds["raw"][k][xd.obs_names].obsm["spatial"]
    a = sc.concat(dict(xds[C.XENIUM_CORRECTION].items()), label="dataset_id", join="outer")
    a.obs[XENIUM_LEVELS] = pd.DataFrame(a.obs["dataset_id"].tolist(), index=a.obs.index, columns=XENIUM_LEVELS)
    a = a[a.obs[LEVEL_KEYS].notna().all(axis=1)].copy()
    return a


def tiles(a: ad.AnnData, panel: str):
    xy = a.obsm["spatial"]
    tx = np.floor(xy[:, 0] / C.XENIUM_TILE_UM).astype(int)
    ty = np.floor(xy[:, 1] / C.XENIUM_TILE_UM).astype(int)
    section = a.obs["sample"].astype(str)
    sec_short = section.str.extract(r"(CH_G_\d+[a-z]_x[0-9r])", expand=False)
    tile_id = sec_short + "_tile_" + pd.Series(tx, index=a.obs_names).astype(str) + "_" + pd.Series(ty, index=a.obs_names).astype(str)
    a.obs["tile_id"] = tile_id.to_numpy()
    a.obs["is_til"] = a.obs[C.XENIUM_TIL_COL].isin(C.XENIUM_TIL_TYPES)

    meta = a.obs.groupby("tile_id", observed=True).agg(
        section=("sample", "first"), n_cells=("tile_id", "size"), til_count=("is_til", "sum")
    )
    meta["section"] = meta.section.astype(str).str.extract(r"(CH_G_\d+[a-z]_x[0-9r])", expand=False)
    meta.insert(0, "sample_id", meta.section.str.extract(r"(CH_G_\d+[a-z])", expand=False))
    meta["panel"] = panel
    txy = meta.index.to_series().str.extract(r"_tile_(-?\d+)_(-?\d+)$").astype(int)
    meta["tile_x"], meta["tile_y"] = txy[0], txy[1]
    meta["x_um"], meta["y_um"] = meta.tile_x * C.XENIUM_TILE_UM, meta.tile_y * C.XENIUM_TILE_UM
    meta["in_paper_analysis"] = meta.n_cells >= C.XENIUM_MIN_CELLS_PER_TILE

    ct = a.obs.groupby("tile_id", observed=True)[C.XENIUM_CELLTYPE_COL].value_counts().unstack(fill_value=0)
    ct = ct.loc[meta.index]

    pb = sc.get.aggregate(a, by="tile_id", func="sum", layer="counts")
    pb = ad.AnnData(X=pb.layers["sum"], obs=meta.loc[pb.obs_names].copy(), var=pd.DataFrame(index=a.var_names))
    pb.uns["description"] = f"Xenium {panel} {C.XENIUM_SEGMENTATION} {C.XENIUM_CORRECTION} counts summed per {C.XENIUM_TILE_UM:g} um tile"
    return pb, meta.join(ct)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="+", default=C.XENIUM_PANELS)
    args = ap.parse_args()
    for panel in args.panels:
        a = load_panel(panel)
        a.layers["counts"] = a.X.copy()
        print(panel, a.shape, "sections", a.obs["sample"].nunique())
        pb, table = tiles(a, panel)
        pb.write_h5ad(out("xenium", f"tile_counts_{panel}.h5ad"), compression="gzip")
        table.index.name = "tile_id"
        table.to_csv(out("xenium", f"tile_celltype_counts_{panel}.csv"))
        print(panel, "tiles", pb.n_obs, "in analysis", int(table.in_paper_analysis.sum()),
              "sections", table.section.nunique(), "patients", table.sample_id.nunique())


if __name__ == "__main__":
    main()
