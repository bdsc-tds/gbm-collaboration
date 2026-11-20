"""snRNA-seq pseudobulk: raw counts summed per sample x level-4 cell type.

Outputs (snrna/):
  pseudobulk_counts_sample_x_celltype.h5ad
      X   : summed raw (unnormalized) counts, int, sparse
      obs : sample_id, celltype_level4_curated, celltype_level{3,2,1}_curated,
            n_cells, n_genes_detected
      var : is_hvg, is_named_gene, named_gene_sources
      genes = top-N seurat_v3 HVGs (computed on cells) U genes named in the manuscript
  n_cells_sample_x_celltype.csv          wide table (0 = cell type absent)
  n_genes_detected_sample_x_celltype.csv wide table, genes with >=1 count over the whole
                                         transcriptome (not only released genes)

Cell-type labels are re-derived from celltype_level4_scanvi through the annotation key,
as in every figure notebook. Doublet labels are excluded (as in the paper's
proportion analyses).
"""

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import config as C
from common import annotation_key, block_of, nnz_by_group, out, sum_by_group

CACHE = C.WORK_DIR / "sc_cohort_cache.h5ad"
HVG_SUBSAMPLE = 400_000


def load():
    if CACHE.exists():
        print("reading cache", CACHE)
        return ad.read_h5ad(CACHE)
    from gbmcollaboration import AllDataCenters

    return AllDataCenters.load_singlecell(**C.LOADER_FLAGS)


def main():
    adata = load()
    print(adata)
    counts = adata.layers["unnormalized_counts"]
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)

    key = annotation_key()
    obs = pd.DataFrame(index=adata.obs_names)
    obs["sample_id"] = block_of(adata.obs["orig.ident"].astype(str)).to_numpy()
    l4s = adata.obs[C.CELLTYPE_SOURCE_COL].astype(str)
    unmapped = sorted(set(l4s) - set(key.index))
    if unmapped:
        sys.exit(f"Unmapped level-4 labels: {unmapped}")
    for lvl in (4, 3, 2, 1):
        obs[f"celltype_level{lvl}_curated"] = key.loc[l4s, f"celltype_level{lvl}_curated"].to_numpy()

    is_doublet = obs[C.CELLTYPE_COL].str.contains("oublet", case=False)
    print(f"cells {len(obs)}, samples {obs.sample_id.nunique()}, doublet cells dropped {int(is_doublet.sum())}")
    keep = ~is_doublet.to_numpy()
    obs, counts = obs[keep], counts[keep]

    # --- HVGs on a random cell subsample (seurat_v3 on raw counts)
    rng = np.random.default_rng(0)
    idx = np.sort(rng.choice(counts.shape[0], size=min(HVG_SUBSAMPLE, counts.shape[0]), replace=False))
    sub = ad.AnnData(counts[idx], var=pd.DataFrame(index=adata.var_names))
    sc.pp.highly_variable_genes(sub, n_top_genes=C.N_HVG, flavor=C.HVG_FLAVOR)
    hvg = sub.var["highly_variable"]

    named = pd.read_csv(out("genes", "named_genes.csv")).set_index("gene")["sources"]
    var = pd.DataFrame(index=adata.var_names)
    var["is_hvg"] = hvg.reindex(var.index).fillna(False).to_numpy()
    var["is_named_gene"] = var.index.isin(named.index)
    var["named_gene_sources"] = named.reindex(var.index).to_numpy()
    keep_genes = (var.is_hvg | var.is_named_gene).to_numpy()
    missing_named = sorted(set(named.index) - set(var.index))
    print(f"genes released: {keep_genes.sum()} (HVG {var.is_hvg.sum()}, named present "
          f"{var.is_named_gene.sum()}/{len(named)}); named absent from snRNA: {len(missing_named)}")

    # --- pseudobulk
    group = obs["sample_id"] + "|" + obs[C.CELLTYPE_COL]
    S, gidx = sum_by_group(counts, group)
    n_cells = group.value_counts().reindex(gidx)
    n_genes = nnz_by_group(counts, group).reindex(gidx)

    pobs = pd.DataFrame(index=gidx)
    pobs[["sample_id", C.CELLTYPE_COL]] = pobs.index.to_series().str.split("|", expand=True).to_numpy()
    lv = obs.drop_duplicates(C.CELLTYPE_COL).set_index(C.CELLTYPE_COL)
    for lvl in (3, 2, 1):
        pobs[f"celltype_level{lvl}_curated"] = lv.loc[pobs[C.CELLTYPE_COL], f"celltype_level{lvl}_curated"].to_numpy()
    pobs["n_cells"] = n_cells.to_numpy()
    pobs["n_genes_detected"] = n_genes.to_numpy()
    pobs.index = pobs.sample_id + "__" + pobs[C.CELLTYPE_COL]

    X = S[:, keep_genes]
    X.data = np.rint(X.data)
    pb = ad.AnnData(X=X.astype(np.int64), obs=pobs, var=var[keep_genes].copy())
    pb.uns["description"] = "snRNA-seq raw counts summed per sample x celltype_level4_curated"
    pb.uns["hvg"] = f"{C.HVG_FLAVOR} n_top_genes={C.N_HVG} on {len(idx)} random cells (seed 0)"
    pb.uns["named_genes_absent"] = missing_named
    pb.write_h5ad(out("snrna", "pseudobulk_counts_sample_x_celltype.h5ad"), compression="gzip")

    for col, name in (("n_cells", "n_cells"), ("n_genes_detected", "n_genes_detected")):
        wide = pobs.pivot(index="sample_id", columns=C.CELLTYPE_COL, values=col).fillna(0).astype(int)
        wide.to_csv(out("snrna", f"{name}_sample_x_celltype.csv"))
    print(pb)
    print(f"samples {pobs.sample_id.nunique()}, cell types {pobs[C.CELLTYPE_COL].nunique()}")


if __name__ == "__main__":
    main()
