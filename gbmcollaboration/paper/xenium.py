"""Xenium tile analysis helpers (Fig 6b, 6g).

Moved verbatim from analysis/xenium/niches.ipynb cell 0 (now analysis/21_xenium_niches.ipynb):
cells are binned in a regular grid of tiles per section (`assign_grid_tiles`), tiles are
pseudobulked and log-normalised (`pseudobulk_lognorm`), and the TIL count of a tile is
correlated with tile-level categories per section (`spearman_tils`, BH within section).
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def pseudobulk_lognorm(adata, by):
    pb = sc.get.aggregate(adata, by=by, func='sum',layer='counts' )
    pb.X = pb.layers['sum']
    sc.pp.normalize_total(pb)
    sc.pp.log1p(pb)
    return pb


def assign_grid_tiles(adata, sample_key='sample', spatial_key='spatial', tile_size=55.0):
    """
    Assigns each cell to a non-overlapping spatial grid tile.
    Adds 'tile_id' to adata.obs.
    """
    x = adata.obsm[spatial_key][:, 0]
    y = adata.obsm[spatial_key][:, 1]

    # Calculate discrete grid indices
    tile_x = np.floor(x / tile_size).astype(int).astype(str)
    tile_y = np.floor(y / tile_size).astype(int).astype(str)

    # Create unique tile IDs. We prepend the sample ID
    samples = adata.obs[sample_key].astype(str)
    tile_ids = samples + "_tile_" + tile_x + "_" + tile_y

    adata.obs['tile_id'] = pd.Categorical(tile_ids)


def spearman_tils(pb, categories, sample_key, til_count_key='TIL_count',
                  min_tiles=10, filter_zero_til=True):
    """
    Computes Spearman correlation per sample and applies FDR correction per sample.
    """
    # 1. Prepare data
    pb_df = pb.obs.copy() # Work directly on obs DataFrame
    if filter_zero_til:
        pb_df = pb_df[pb_df[til_count_key] > 0]

    # 2. Correlate per sample
    results = []
    for sample_id, group in pb_df.groupby(sample_key, observed=True):
        if len(group) < min_tiles:
            continue

        for cat in categories:
            rho, p_val = spearmanr(group[til_count_key], group[cat], nan_policy='omit')
            results.append({
                sample_key: sample_id,
                'category': cat,
                'n_tiles': len(group),
                'rho': rho,
                'p': p_val
            })

    res_df = pd.DataFrame(results)


    # 3. Correct p-values per sample
    def fdr_correct(group):
        # Drop NaNs before correction
        valid = group['p'].notna()
        p_adj = np.full(len(group), np.nan)
        if valid.any():
            _, corrected, _, _ = multipletests(group.loc[valid, 'p'], method='fdr_bh')
            p_adj[valid] = corrected
        group['p_adj'] = p_adj
        return group

    res_df = res_df.groupby(sample_key, group_keys=False).apply(fdr_correct)

    return res_df
