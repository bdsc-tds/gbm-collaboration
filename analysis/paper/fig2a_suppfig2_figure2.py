"""Fig 2a and the Supp Fig 2 UMAPs and per-patient composition bars (snRNA-seq atlas).

Produces (OUT = figures/paper/fig2a_suppfig2_figure2/, and figures/paper/panels/):
- Fig 2a: UMAP of all annotated nuclei coloured by level-4 cell type (fig2a_umap_level4) and
  its coordinates (fig2a_umap_coordinates.csv).
- Supp Fig 2a / 2d / 2g: UMAPs of malignant (BBKNN), immune and stromal + neuronal nuclei,
  levels 1-4 (supp2a_malignant_umap_<lvl>_bbknn, supp2d_immune_umap_<lvl>,
  supp2g_stromal_neuronal_umap_<lvl>).
- Supp Fig 2b / 2e / 2h: per-patient composition stacked bars of the same lineages, levels 1-4
  (supp2b_malignant_stacked_bar_<lvl>, supp2e_immune_stacked_bar_<lvl>,
  supp2h_stromal_neuronal_stacked_bar_<lvl>).
The Fig 2b-c and Supp Fig 2c/f/i dot plots are in fig2bc_suppfig2_marker_dotplots.ipynb.

Origin: figures/figure2.py (Owkin). Published panels: EFS figures/figure2/ (old names fig2a,
fig2c-h; fig2c-h are now Supp Fig 2a, 2b, 2d, 2e, 2g, 2h).

Comparison: level-1 and level-4 labels are those of the published panels (same nuclei, same
Fig 2a); with the repo annotation key some level-2/3 names differ (Mono_Mac, cDC1 / cDC2).
The stacked-bar tables and UMAP coordinates are identical to the previous run
(figures/paper/03_figure2/).
"""

# %% -- Import packages
import gc
import time

import bbknn
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from gbmcollaboration import AllDataCenters
from gbmcollaboration.paper import annotation, output, snrna

# %% -- Parameters
OUT = output.setup("fig2a_suppfig2_figure2")

SEED = 123
FILTERING_SAMPLES = ["ER_G_743a"]
BATCH_KEY = "orig.ident"
LEVELS = ["level1", "level2", "level3", "level4"]

T0 = time.time()


# -- helper functions
def compute_pca_umap(
    adata,
    batch_key=None,
    n_neighbors=15,
    min_dist=0.3,
    seed=SEED,
):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, random_state=seed)

    if batch_key is None:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    else:
        bbknn.bbknn(adata, batch_key=batch_key)

    sc.tl.umap(adata, min_dist=min_dist, random_state=seed)
    return adata


def downsample_stratified(
    adata,
    group_keys,
    frac=0.5,
    seed=123,
    min_cells_per_group=10,
):
    """
    Stratified downsampling for AnnData.
    Keeps a fraction of cells per stratum defined by group_keys.
    Uses integer positions (robust to string obs_names).
    """
    if isinstance(group_keys, str):
        group_keys = [group_keys]

    rng = np.random.default_rng(seed)

    # Build group IDs (one per cell)
    group_id = adata.obs[group_keys].astype(str).agg("||".join, axis=1)

    # Work with integer positions to avoid string-index issues
    pos = np.arange(adata.n_obs, dtype=int)
    tmp = pd.DataFrame({"pos": pos, "_group_id": group_id.values})

    keep_pos = []

    for gid, sub in tmp.groupby("_group_id"):
        idx = sub["pos"].to_numpy(dtype=int)
        n = idx.size

        if n < min_cells_per_group:
            chosen = idx  # keep all small strata
        else:
            k = int(np.ceil(frac * n))
            k = max(1, k)
            chosen = rng.choice(idx, size=k, replace=False)

        keep_pos.append(chosen)

    keep_pos = np.concatenate(keep_pos).astype(int)
    keep_mask = np.zeros(adata.n_obs, dtype=bool)
    keep_mask[keep_pos] = True

    return adata[keep_mask].copy()


def save_umap(
    adata,
    color,
    stem,
    title=None,
    legend="right margin",
    size=4,
    palette=None,
    cmap=None,
    panel=None,
):
    sc.pl.umap(
        adata,
        color=color,
        title=title,
        legend_loc=legend,
        show=False,
        size=size,
        palette=palette,   # catégories
        cmap=cmap,         # continu
    )
    output.save(OUT, stem, fig=plt.gcf(), panel=panel)
    plt.close()


def stacked_bar(obs, reference_column, title, legend_title, stem, panel=None):
    """STRUCTURE-like stacked bar: per-patient fractions, patients ordered by Ward clustering
    of their composition."""
    # 1) counts -> fractions per patient
    prop = (
        obs
        .groupby([BATCH_KEY, reference_column])
        .size()
        .unstack(fill_value=0)
    )
    prop = prop.div(prop.sum(axis=1), axis=0)

    # 2) order patients by hierarchical clustering on composition
    Z = linkage(
        pdist(prop.values, metric="euclidean"),
        method="ward",
    )
    ordered_idx = leaves_list(Z)
    prop_ord = prop.iloc[ordered_idx]

    # 3) plot (no x labels, packed bars)
    fig, ax = plt.subplots(figsize=(12, 3))

    prop_ord.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        width=1.0,     # packed like STRUCTURE
        linewidth=0,
    )

    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("Fraction")
    ax.set_title(title)

    ax.legend(
        title=legend_title,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
    )

    plt.tight_layout()
    output.save(OUT, stem, fig=fig, data=prop_ord, panel=panel)
    plt.close(fig)
    return prop_ord


def lineage_panels(cells, name, letters, bar_title, legend_title, umap_title, batch_key=None):
    """Supp Fig 2 stacked bars (levels 1-4), then the lineage UMAP (levels 1-4)."""
    bar_letter, umap_letter = letters
    for lvl in LEVELS:
        print(f"{name} stack bar {lvl}")
        stacked_bar(cells.obs, f"celltype_{lvl}_curated", bar_title, legend_title,
                    f"supp2{bar_letter}_{name}_stacked_bar_{lvl}", panel=f"Supp Fig 2{bar_letter} {lvl}")
    print(f"Starting {name} UMAP")
    cells = compute_pca_umap(cells, batch_key=batch_key, seed=SEED)
    suffix = "_bbknn" if batch_key else ""
    for lvl in LEVELS:
        print(f"{name}-umap {lvl}")
        save_umap(cells, color=f"celltype_{lvl}_curated", stem=f"supp2{umap_letter}_{name}_umap_{lvl}{suffix}",
                  title=umap_title, panel=f"Supp Fig 2{umap_letter} {lvl}")


# %% -- Load data
adata = AllDataCenters.load_singlecell(
    exclude_idh1_mut=True,
    exclude_patient_ids = True,
    exclude_relapse = True,
    )
print(f"Loaded {adata.n_obs:,} nuclei x {adata.n_vars:,} genes, {adata.obs[BATCH_KEY].nunique()} samples")

# create simplified annotations from lvl1 to lvl4
snrna.add_curated_levels(adata, annotation.celltype_key(), LEVELS)

# %% -- RESET embeddings (force recompute)
for key in ["X_pca", "X_umap"]:
    if key in adata.obsm:
        del adata.obsm[key]

for key in ["pca", "neighbors", "umap"]:
    if key in adata.uns:
        del adata.uns[key]

# -- Drop samples with too few cells
adata = adata[
    ~adata.obs[BATCH_KEY].isin(FILTERING_SAMPLES)
]
# -- Drop "Other" and "Doublet"
adata_filt = adata[
    ~adata.obs["celltype_level1_curated"].isin(["Other", "Doublet", "Doublet_TNK_nos"])
].copy()
del adata
gc.collect()

# -- Drop the Malignant_other_1
n_before = adata_filt.n_obs
adata_filt = adata_filt[
    ~adata_filt.obs["celltype_level4_curated"].isin(["Malignant_other_1"])
]
print(f"Malignant_other_1 filter removed {n_before - adata_filt.n_obs} nuclei")
print(f"Annotated nuclei: {adata_filt.n_obs:,} from {adata_filt.obs[BATCH_KEY].nunique()} samples, "
      f"{adata_filt.obs['celltype_level4_curated'].nunique()} level-4 cell types")
print(adata_filt.obs["celltype_level1_curated"].value_counts().to_string())

############################################################
############### ALL CELLS ##################################
############################################################
# %% Fig 2a — Single cell UMAP, level 4
adata_umap = compute_pca_umap(adata_filt.copy(), seed=SEED)

save_umap(
    adata_umap,
    color="celltype_level4_curated",
    stem="fig2a_umap_level4",
    title="Single-cell MOSAIC GBM UMAP",
    panel="Fig 2a",
)
umap_coords = pd.DataFrame(adata_umap.obsm["X_umap"], index=adata_umap.obs_names, columns=["UMAP1", "UMAP2"])
umap_coords[BATCH_KEY] = adata_umap.obs[BATCH_KEY].astype(str).values
umap_coords["celltype_level4_curated"] = adata_umap.obs["celltype_level4_curated"].astype(str).values
output.save_table(OUT, "fig2a_umap_coordinates", umap_coords)
del adata_umap, umap_coords
gc.collect()
print(f"Fig 2a done ({(time.time() - T0) / 60:.0f} min)")

############################################################
############### MALIGNANT ##################################
############################################################
# %% Supp Fig 2a-b — Subset tumor cells (keep samples only if more than 30 cells)
tumor_cells = adata_filt[adata_filt.obs["celltype_level1_curated"] == "Malignant"].copy()
min_cells = 30
counts = tumor_cells.obs[BATCH_KEY].value_counts()
keep = counts[counts >= min_cells].index
tumor_cells_f = tumor_cells[tumor_cells.obs[BATCH_KEY].isin(keep)].copy()

# Downsample: keep 50% per patient x level4 because the adata_filt is too big to use bbknn
tumor_cells_ds = downsample_stratified(
    tumor_cells_f,
    group_keys=[BATCH_KEY, "celltype_level4_curated"],  # patient x level4
    frac=0.5,
    seed=SEED,
    min_cells_per_group=10,
)
print("Tumor cells before:", tumor_cells.n_obs)
print("Tumor cells after :", tumor_cells_ds.n_obs)
del tumor_cells, tumor_cells_f
gc.collect()

# BBKNN pipeline on the downsampled data (patient-normalized)
lineage_panels(tumor_cells_ds, "malignant", ("b", "a"), "Tumor cell composition by patient",
               "Tumor cell state", "Tumor cells UMAP (BBKNN)", batch_key=BATCH_KEY)
del tumor_cells_ds
gc.collect()

############################################################
############### IMMUNE #####################################
############################################################
# %% Supp Fig 2d-e — immune
immune_cells = adata_filt[
    adata_filt.obs["celltype_level1_curated"].isin(["Myeloid", "Lymphoid"])
].copy()
lineage_panels(immune_cells, "immune", ("e", "d"), "Immune cell composition by patient",
               "Immune cell state", "Immune cell UMAP")
del immune_cells
gc.collect()

############################################################
############### STROMAL / NEURONAL #########################
############################################################
# %% Supp Fig 2g-h — stromal neuronal
sn_cells = adata_filt[
    adata_filt.obs["celltype_level1_curated"].isin(["Glia", "Neuron", "Vascular"])
].copy()
lineage_panels(sn_cells, "stromal_neuronal", ("h", "g"), "Stromal + Neuronal composition by patient",
               "Stromal/Neuronal state", "Stromal + Neuronal UMAP")
del sn_cells
gc.collect()
print(f"Done ({(time.time() - T0) / 60:.0f} min)")

# %%
print(output.versions(("scanpy", "anndata", "bbknn", "umap-learn", "pynndescent", "scikit-learn",
                       "pandas", "numpy")).to_string(index=False))
