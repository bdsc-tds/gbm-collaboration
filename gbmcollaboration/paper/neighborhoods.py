"""ECM and cell-type neighborhoods: code shared by fig4cdg_suppfig4f_ECM_niches_v4,
fig4ab_suppfig4abg_Cell_Neighborhoods and fig4e_suppfig4e_colocalization_plots.

- ``load_spots``: spots, coordinates and cell2location q05 abundances (``deconv_cohort_level_4``)
  of the Visium cohort, as the ``adata.obs`` of ``AllDataCenters.load_spatial_adata`` without
  the expression matrix (``AllDataCenters.load_spatial(load_gex=False)``). The gene-subset
  path of ``load_spatial_adata`` (``genes=``) returns neither coordinates nor deconvolution
  (its spots and spot order are the same: checked in fig4cdg).
- ``pseudobulk_sums``: per section and label (e.g. ECM neighborhood) sums of every gene of the
  log-normalised matrix, read a few sections at a time with each centre's ``load_spatial``.
  The full dense matrix (870,077 spots x 18,082 genes) does not fit in memory, and
  ``load_spatial_adata(genes=...)`` would hold all genes of all spots at once.
- ``c2l_curated_counts``: maps the cell2location mean abundances
  (``meanscell_abundance_w_sf_<level 4>`` columns) to a curated level of the annotation key,
  drops doublet / other types and sums per curated type.
"""

import gc

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from gbmcollaboration import AllDataCenters

from . import config as C
from .visium import C2L_PREFIX

_PREFIXES = {"mosaic_chuv": "CH", "mosaic_uker": "ER", "mosaic_cr": "CR"}


def load_spots(deconvolution_slot_name: str = "deconv_cohort_level_4") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Spots of the Visium cohort with the loader exclusions of the notebooks.

    Returns ``obs`` (index ``<section>_<barcode>``; columns ``sample``, ``imagerow``,
    ``imagecol``) and ``q05`` (cell2location q05 abundances per level-4 scANVI type, the
    columns that ``load_spatial_adata`` put in ``adata.obs``), both in the spot order of
    ``load_spatial_adata`` (the order of the frozen neighborhood files).
    """
    loaded = AllDataCenters.load_spatial(exclude_low_quality_visium=True, **C.LOADER_FLAGS,
                                         normalization="log_norm_count", radius=1,
                                         deconvolution_slot_name=deconvolution_slot_name,
                                         load_gex=False, load_deconvolution=True)
    df = pd.concat([o.celltype_df for o in loaded.values()])
    index = pd.Index([f"{s}_{b}" for s, b in df.index], name="spot_id")
    obs = pd.DataFrame({"sample": df.index.get_level_values(0).astype(str),
                        "imagerow": df["imagerow"].to_numpy(), "imagecol": df["imagecol"].to_numpy()}, index=index)
    q05 = df.drop(columns=["imagerow", "imagecol"])
    q05.index = index
    return obs, q05


def _sums_one_batch(dc, sections, labels, normalization):
    obj = dc.load_spatial(exclude_low_quality_visium=True, **C.LOADER_FLAGS, list_sample_id=list(sections),
                          normalization=normalization, radius=1, load_deconvolution=False)
    df = obj.genecounts_df
    genes = [g for g in df.columns if g not in ("imagerow", "imagecol")]
    X = df[genes].to_numpy()
    spot = np.array([f"{s}_{b}" for s, b in df.index])
    sample = df.index.get_level_values(0).astype(str).to_numpy()
    lab = labels.reindex(spot).to_numpy()
    del obj, df
    sums, n = {}, {}
    for key in pd.unique(pd.Series(list(zip(sample, lab)))):
        if pd.isna(key[1]):
            continue
        m = (sample == key[0]) & (lab == key[1])
        sums[key] = np.nan_to_num(X[m]).sum(axis=0, dtype=np.float64)
        n[key] = int(m.sum())
    del X
    gc.collect()
    return pd.DataFrame(sums, index=genes).T, pd.Series(n)


def pseudobulk_sums(labels: pd.Series, normalization: str = "log_norm_count", section_ids=None,
                    batch: int = 4, n_jobs: int = 8) -> tuple[pd.DataFrame, pd.Series]:
    """Sums of every gene per (section, label), and the number of spots of each.

    ``labels``: one label per spot (index ``<section>_<barcode>``). Genes are those of
    ``load_spatial_adata``: within a centre, the union over its sections (0 where a section
    lacks a gene); across centres, the intersection. Sums are accumulated in float64.
    Returns (sums: rows MultiIndex (sample, label), columns genes; n: spots per row).
    Batches of ``batch`` sections are read in ``n_jobs`` processes (a few GB each).
    """
    if section_ids is None:
        section_ids = sorted(labels.index.str.rsplit("_", n=1).str[0].unique())
    jobs = []
    for name, dc in AllDataCenters.data_centers.items():
        prefix = _PREFIXES[name]
        secs = sorted(s for s in section_ids if s.startswith(prefix))
        jobs += [(prefix, dc, secs[i:i + batch]) for i in range(0, len(secs), batch)]
    res = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_sums_one_batch)(dc, secs, labels, normalization) for _, dc, secs in jobs)
    centre_genes = {}
    for (prefix, _, _), (s, _) in zip(jobs, res):
        centre_genes.setdefault(prefix, set()).update(s.columns)
    genes = sorted(set.intersection(*centre_genes.values()))
    sums = pd.concat([s.reindex(columns=genes, fill_value=0.0) for s, _ in res])
    n = pd.concat([k for _, k in res])
    sums.index = n.index = pd.MultiIndex.from_tuples(sums.index, names=["sample", "label"])
    return sums.sort_index(), n.sort_index()


def c2l_curated_counts(obs: pd.DataFrame, celltype_key: pd.DataFrame, curation_col: str = "celltype_level4_curated",
                       exclude_patterns=("doublet", "other"), c2l_prefix: str = C2L_PREFIX, verbose: bool = True):
    """cell2location mean abundances summed per curated cell type (spots x curated types).

    Level-4 columns ``<c2l_prefix><type>`` of ``obs`` that are in the annotation key, mapped
    to ``curation_col``; curated names containing one of ``exclude_patterns``
    (case-insensitive) are dropped before summing.
    """
    c2l_cols = [c for c in obs.columns if c.startswith(c2l_prefix)]
    lvl4_names = pd.Index([c.replace(c2l_prefix, "") for c in c2l_cols])
    counts = obs[c2l_cols].copy()
    counts.columns = lvl4_names

    # Map to curated names
    lvl4_in_key = lvl4_names.intersection(celltype_key.index)
    counts = counts[lvl4_in_key]
    curated_map = celltype_key.loc[lvl4_in_key, curation_col].astype(str).fillna("Unknown")

    # Filter out excluded cell types before aggregation
    keep_mask = ~curated_map.str.lower().apply(lambda x: any(pat in x for pat in exclude_patterns))
    counts = counts.loc[:, keep_mask.values]
    curated_map = curated_map[keep_mask]
    if verbose:
        print(f"Excluded cell types: {(~keep_mask).sum()}")
        print(f"Retained cell types: {keep_mask.sum()}")
    return counts.T.groupby(curated_map).sum().T
