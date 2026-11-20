"""Visium pseudobulks per sample x ECM neighbourhood and per sample x TIL zone.

Inputs: raw Space Ranger counts (AllDataCenters.load_spatial_adata, normalization=
"raw_count", QC-passing 228 samples), frozen spot annotations used in the paper
(ECM neighbourhood v4 and TIL zones = Supp Table 4a), cell2location abundances
(obsm['means_cell_abundance_w_sf'] per sample on S3, as read in ECM_niches_v4).

Outputs (visium/):
  pseudobulk_counts_sample_x_ecm_neighborhood.h5ad
  pseudobulk_counts_sample_x_til_zone.h5ad
      X   : summed raw counts (int, sparse); genes = seurat_v3 HVGs (on spots) U named genes
      obs : sample_id, visium_section, region, n_spots, til_group
  n_spots_sample_x_ecm_neighborhood.csv / n_spots_sample_x_til_zone.csv
  c2l_proportions_sample_x_ecm_neighborhood.csv / c2l_proportions_sample_x_til_zone.csv
      cell2location abundances summed over the region's spots and row-normalised to 1
      (long format: sample_id, region, n_spots, then one column per level-4 cell type)
"""

import gc
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from gbmcollaboration import AllDataCenters

import config as C
from common import block_of, out, sum_by_group, til_labels

HVG_SUBSAMPLE = 300_000
C2L_PREFIX = "meanscell_abundance_w_sf_"


def load_c2l(samples):
    """Cell2location abundances (means) per spot, indexed by <section>_<barcode>, one column per cell type."""
    c2l = AllDataCenters.load_cell2location(sorted(samples)).drop(columns="total_cells")
    c2l.index = [f"{s}_{b}" for s, b in c2l.index]
    c2l.columns = c2l.columns.str.replace(C2L_PREFIX, "", regex=False)
    return c2l


def pseudobulk(counts, spots: pd.DataFrame, region_col, var, til, name):
    grp = spots["visium_section"] + "|" + spots[region_col].astype(str)
    S, gidx = sum_by_group(counts, grp)
    obs = pd.DataFrame(index=gidx)
    obs[["visium_section", "region"]] = obs.index.to_series().str.split("|", expand=True).to_numpy()
    obs.insert(0, "sample_id", block_of(obs.visium_section))
    obs["n_spots"] = grp.value_counts().reindex(gidx).to_numpy()
    obs["til_group"] = til["til_group"].reindex(obs.sample_id).to_numpy()
    obs.index = obs.sample_id + "__" + obs.region
    keep = (var.is_hvg | var.is_named_gene).to_numpy()
    X = S[:, keep]
    X.data = np.rint(X.data)
    pb = ad.AnnData(X=X.astype(np.int64), obs=obs, var=var[keep].copy())
    pb.uns["region_definition"] = region_col
    pb.write_h5ad(out("visium", f"pseudobulk_counts_sample_x_{name}.h5ad"), compression="gzip")
    obs.pivot(index="sample_id", columns="region", values="n_spots").fillna(0).astype(int).to_csv(
        out("visium", f"n_spots_sample_x_{name}.csv")
    )
    print(name, pb.shape, "samples", obs.sample_id.nunique(), "regions", sorted(obs.region.unique()))
    return obs


def c2l_props(c2l: pd.DataFrame, spots: pd.DataFrame, region_col, name):
    grp = spots["visium_section"] + "|" + spots[region_col].astype(str)
    summed = c2l.groupby(grp.to_numpy()).sum()
    prop = summed.div(summed.sum(axis=1), axis=0)
    idx = prop.index.to_series().str.split("|", expand=True)
    prop.insert(0, "n_spots", grp.value_counts().reindex(prop.index).to_numpy())
    prop.insert(0, "region", idx[1].to_numpy())
    prop.insert(0, "sample_id", block_of(idx[0]))
    prop.sort_values(["sample_id", "region"]).to_csv(
        out("visium", f"c2l_proportions_sample_x_{name}.csv"), index=False, float_format="%.6g"
    )


def load_raw_counts(sections, batch=25):
    """Raw counts for the given Visium sections, loaded per centre in batches and kept
    sparse (load_spatial_adata builds dense frames for the whole cohort and runs out of
    memory). Same loader/filters as load_spatial_adata; genes = intersection across
    centres, as there."""
    from gbmcollaboration import CHUVDataCenter, CRDataCenter, UKERDataCenter

    parts = []
    for prefix, dc in (("CH", CHUVDataCenter), ("ER", UKERDataCenter), ("CR", CRDataCenter)):
        secs = sorted(s for s in sections if s.startswith(prefix))
        for i in range(0, len(secs), batch):
            obj = dc.load_spatial(
                exclude_low_quality_visium=True, **C.LOADER_FLAGS, list_sample_id=secs[i : i + batch],
                normalization="raw_count", load_deconvolution=False,
            )
            df = obj.genecounts_df
            genes = df.columns.difference(["imagerow", "imagecol"])
            ids = pd.Index([f"{s}_{b}" for s, b in df.index])
            parts.append((sp.csr_matrix(df[genes].to_numpy(dtype=np.float32)), pd.Index(genes), ids))
            del obj, df
            gc.collect()
            print(f"  {prefix}: loaded {min(i + batch, len(secs))}/{len(secs)} sections")
    common = parts[0][1]
    for _, g, _ in parts[1:]:
        common = common.intersection(g)
    common = common.sort_values()
    mats = [m[:, g.get_indexer(common)] for m, g, _ in parts]
    return sp.vstack(mats).tocsr(), common, pd.Index(np.concatenate([ids for _, _, ids in parts]))


def load_spots():
    ecm = pd.read_csv(C.ECM_NEIGHBORHOOD_CSV, index_col="spot_id")
    tz = pd.read_csv(C.TIL_ZONE_CSV, index_col="spot_id")
    assert ecm.index.equals(tz.index) or set(ecm.index) == set(tz.index)
    return pd.DataFrame(
        {
            "visium_section": ecm["sample"].to_numpy(),
            "ecm_neighborhood": "ECM_" + ecm["ECM_neighborhood"].astype(str).to_numpy(),
            "til_zone": tz.loc[ecm.index, "TIL_zone"].to_numpy(),
        },
        index=ecm.index,
    )


def run_c2l(spots):
    c2l = load_c2l(set(spots.visium_section))
    n_missing = len(spots.index.difference(c2l.index))
    print(f"spots absent from cell2location files: {n_missing}")
    if n_missing:
        raise SystemExit("cell2location spot IDs do not match annotations")
    c2l = c2l.reindex(spots.index)
    na_cols = c2l.columns[c2l.isna().any()].tolist()
    print(f"cell types missing in some samples (set to 0 there): {na_cols}")
    c2l = c2l.fillna(0)
    c2l_props(c2l, spots, "ecm_neighborhood", "ecm_neighborhood")
    c2l_props(c2l, spots, "til_zone", "til_zone")


def main():
    spots = load_spots()
    if "--c2l-only" in sys.argv:
        return run_c2l(spots)
    counts, var_names, spot_ids = load_raw_counts(set(spots.visium_section))
    print("spots", len(spot_ids), "genes", len(var_names))

    missing = spots.index.difference(spot_ids)
    if len(missing):
        raise SystemExit(f"{len(missing)} annotated spots missing from counts")
    extra = spot_ids.difference(spots.index)
    print(f"spots in counts without annotation (dropped): {len(extra)}")
    pos = spot_ids.get_indexer(spots.index)
    counts = counts[pos]

    # HVGs on spots
    rng = np.random.default_rng(0)
    idx = np.sort(rng.choice(counts.shape[0], size=min(HVG_SUBSAMPLE, counts.shape[0]), replace=False))
    sub = ad.AnnData(counts[idx], var=pd.DataFrame(index=var_names))
    sc.pp.highly_variable_genes(sub, n_top_genes=C.N_HVG, flavor=C.HVG_FLAVOR)
    named = pd.read_csv(out("genes", "named_genes.csv")).set_index("gene")["sources"]
    var = pd.DataFrame(index=var_names)
    var["is_hvg"] = sub.var["highly_variable"].to_numpy()
    var["is_named_gene"] = var.index.isin(named.index)
    var["named_gene_sources"] = named.reindex(var.index).to_numpy()
    print(f"HVG {var.is_hvg.sum()}, named present {var.is_named_gene.sum()}/{len(named)}")

    til = til_labels()
    pseudobulk(counts, spots, "ecm_neighborhood", var, til, "ecm_neighborhood")
    pseudobulk(counts, spots, "til_zone", var, til, "til_zone")

    run_c2l(spots)


if __name__ == "__main__":
    main()
