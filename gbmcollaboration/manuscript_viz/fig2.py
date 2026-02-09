"""Scripts to generate subfigures of figure 2 in pdf for the manuscript."""

# %% -- Import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from pathlib import Path
import anndata as ad
import bbknn
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform
import plotly.graph_objects as go

from gbmcollaboration import AllDataCenters
from gbmcollaboration import (
    PATIENTS_IDH1_MUT,
)

SEED = 123
FILTERING_SAMPLES = ["ER_G_743a"]

# -- Global figure parameters (publication-ready)
BLUE_ORANGE_SOFT = LinearSegmentedColormap.from_list(
    "blue_orange_soft",
    ["#4575b4", "#e0f3f8", "#fee090", "#d73027"],
)

sc.set_figure_params(
    dpi=300,              # notebook
    dpi_save=700,        # export high resolution
    fontsize=12,
    frameon=False,
    vector_friendly=True # avoid artefacts
)

plt.rcParams.update({
    "pdf.fonttype": 42,   # TrueType fonts (Illustrator-friendly)
    "ps.fonttype": 42,
})

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


def save_umap(
    adata,
    color,
    path,
    title=None,
    legend="right margin",
    size=4,
    palette=None,
    cmap=None,
    dpi=1000,
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
    plt.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close()

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

############################################################
############### UMAPs ##############################
############################################################
# %% Params
RESULTS_DIR = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/figures/figure2")
SIMPLIFIED_ANNOTATION_MAPPING = Path("/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/resources/SingleCell_Annotation_Key_v3.csv")
REFERENCE_COLUMN = "celltype_level4_scanvi"
BATCH_KEY = "orig.ident"

# -- Create results dir
RESULTS_DIR.mkdir(parents=True,exist_ok=True)

# %% -- Load data
adata = AllDataCenters.load_singlecell(
    exclude_idh1_mut=True,
    exclude_patient_ids = True,
    exclude_relapse = True,
    )

# create simplified annotations from lvl1 to lvl4
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"Mapping {lvl}")
    ct_mapping = (
        pd.read_csv(SIMPLIFIED_ANNOTATION_MAPPING)
        .set_index(REFERENCE_COLUMN)[f"celltype_{lvl}_curated"])
    adata.obs[f"celltype_{lvl}_curated"] = adata.obs[REFERENCE_COLUMN].astype(str).replace(ct_mapping)

# %% -- RESET embeddings (force recompute)
for key in ["X_pca", "X_umap"]:
    if key in adata.obsm:
        del adata.obsm[key]

for key in ["pca", "neighbors", "umap"]:
    if key in adata.uns:
        del adata.uns[key]

# %% -- Drop samples with too few cells
adata = adata[
    ~adata.obs[BATCH_KEY].isin(FILTERING_SAMPLES)
]

# %% Drop "Other" and "Doublet"
adata_filt = adata[
    ~adata.obs["celltype_level1_curated"].isin(["Other", "Doublet", "Doublet_TNK_nos"])
].copy()

# drop the malignant_other_1
adata_filt = adata_filt[
    ~adata_filt.obs["celltype_level4_curated"].isin(["malignant_other_1"])
]

# %% Fig 2a — Single cell UMAP detailed
adata_filt = compute_pca_umap(adata_filt)

for lvl in ["level1", "level2", "level3", "level4"]:
    save_umap(
        adata_filt,
        color=f"celltype_{lvl}_curated",
        path=RESULTS_DIR / f"fig2a_umap_detailed_{lvl}.pdf",
        title="Single-cell MOSAIC GBM UMAP",
    )

# %% Subset tumor cells
# Subset tumor cells (keep samples only if more than 30 cells)
tumor_cells = adata[adata.obs["celltype_level1_curated"] == "Malignant"].copy()
min_cells = 30
counts = tumor_cells.obs[BATCH_KEY].value_counts()
keep = counts[counts >= min_cells].index
tumor_cells_f = tumor_cells[tumor_cells.obs[BATCH_KEY].isin(keep)].copy()

# Downsample: keep 50% per patient x level4 because the adata is too big to use bbknn
tumor_cells_ds = downsample_stratified(
    tumor_cells_f,
    group_keys=[BATCH_KEY, "celltype_level4_curated"],  # patient x level4
    frac=0.5,
    seed=SEED,
    min_cells_per_group=10,
)
print("Tumor cells before:", tumor_cells.n_obs)
print("Tumor cells after :", tumor_cells_ds.n_obs)

# %%
print("Starting tumor-UMAP")
# Now run BBKNN pipeline on downsampled data
tumor_cells_ds = compute_pca_umap(
    tumor_cells_ds,
    batch_key=BATCH_KEY  # BBKNN patient-normalized
)

# %% -- plot
# be careful to add or remove the bbknn suffix in the file name
# when bbknn was used or not
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"tumor-umap {lvl}")
    save_umap(
        tumor_cells_f,
        color=f"celltype_{lvl}_curated",
        path=RESULTS_DIR / f"fig2c_tumor_umap_{lvl}_bbknn.pdf",
        title="Tumor cells UMAP (BBKNN)",
    )

# %% Fig 2d — Tumor stacked barplot (STRUCTURE-like)
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"Tumor stack bar {lvl}")
    REFERENCE_COLUMN = f"celltype_{lvl}_curated"

    tumor_prop = (
        tumor_cells_f.obs
        .groupby([BATCH_KEY, REFERENCE_COLUMN])
        .size()
        .unstack(fill_value=0)
    )

    # Convert to fractions
    tumor_prop = tumor_prop.div(tumor_prop.sum(axis=1), axis=0)

    # Hierarchical clustering on composition profiles
    Z = linkage(
        pdist(tumor_prop.values, metric="euclidean"),
        method="ward",
    )

    # Order patients to maximize similarity between neighbors
    ordered_idx = leaves_list(Z)
    tumor_prop_ord = tumor_prop.iloc[ordered_idx]

    fig, ax = plt.subplots(figsize=(12, 3))

    tumor_prop_ord.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        width=1.0,
        linewidth=0,
    )

    # Remove x labels and ticks
    ax.set_xticks([])
    ax.set_xlabel("")

    ax.set_ylabel("Fraction")
    ax.set_title("Tumor cell composition by patient")

    # Optional: legend outside
    ax.legend(
        title="Tumor cell state",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
    )

    plt.tight_layout()
    plt.savefig(
        RESULTS_DIR / f"fig2d_tumor_stacked_bar_{lvl}.pdf",
        dpi=1000,
        bbox_inches="tight",
    )
    plt.close()


# %% -- Fig2e immune umap
immune_cells = adata[
    adata.obs["celltype_level1_curated"].isin(["Myeloid", "Lymphoid"])
].copy()

immune_cells = compute_pca_umap(immune_cells)

for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"Immune-umap {lvl}")
    save_umap(
        immune_cells,
        color=f"celltype_{lvl}_curated",
        path=RESULTS_DIR / f"fig2e_immune_umap_{lvl}.pdf",
        title="Immune cell UMAP",
    )

# %% Fig 2f — Immune stacked barplot (STRUCTURE-like)
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"Immune stack bar {lvl}")
    REFERENCE_COLUMN = f"celltype_{lvl}_curated"
    # 1) counts -> fractions per patient
    immune_prop = (
        immune_cells.obs
        .groupby([BATCH_KEY, REFERENCE_COLUMN])
        .size()
        .unstack(fill_value=0)
    )
    immune_prop = immune_prop.div(immune_prop.sum(axis=1), axis=0)

    # 2) order patients by hierarchical clustering on composition
    Z = linkage(
        pdist(immune_prop.values, metric="euclidean"),
        method="ward",
    )
    ordered_idx = leaves_list(Z)
    immune_prop_ord = immune_prop.iloc[ordered_idx]

    # 3) plot (no x labels, packed bars)
    fig, ax = plt.subplots(figsize=(12, 3))

    immune_prop_ord.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        width=1.0,     # packed like STRUCTURE
        linewidth=0,
    )

    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("Fraction")
    ax.set_title("Immune cell composition by patient")

    ax.legend(
        title="Immune cell state",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
    )

    plt.tight_layout()
    plt.savefig(
        RESULTS_DIR / f"fig2f_immune_stacked_bar_{lvl}.pdf",
        dpi=1000,
        bbox_inches="tight",
    )
    plt.close()

# %% -- Fig2g stromal neuronal umap
sn_cells = adata[
    adata.obs["celltype_level1_curated"].isin(["Glia", "Neuron", "Vascular"])
].copy()

sn_cells = compute_pca_umap(sn_cells)

# %%
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"stromal-umap {lvl}")
    save_umap(
        sn_cells,
        color=f"celltype_{lvl}_curated",
        path=RESULTS_DIR / f"fig2g_stromal_neuronal_umap_{lvl}.pdf",
        title="Stromal + Neuronal UMAP",
    )

# %% Fig 2h — Stromal + Neuronal stacked barplot (STRUCTURE-like)
for lvl in ["level1", "level2", "level3", "level4"]:
    print(f"Stromal stack bar {lvl}")
    REFERENCE_COLUMN = f"celltype_{lvl}_curated"

    # 1) counts -> fractions per patient
    sn_prop = (
        sn_cells.obs
        .groupby([BATCH_KEY, REFERENCE_COLUMN])
        .size()
        .unstack(fill_value=0)
    )
    sn_prop = sn_prop.div(sn_prop.sum(axis=1), axis=0)

    # 2) order patients by hierarchical clustering on composition
    Z = linkage(
        pdist(sn_prop.values, metric="euclidean"),
        method="ward",
    )
    ordered_idx = leaves_list(Z)
    sn_prop_ord = sn_prop.iloc[ordered_idx]

    # 3) plot (no x labels, packed bars)
    fig, ax = plt.subplots(figsize=(12, 3))

    sn_prop_ord.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        width=1.0,
        linewidth=0,
    )

    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("Fraction")
    ax.set_title("Stromal + Neuronal composition by patient")

    ax.legend(
        title="Stromal/Neuronal state",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
    )

    plt.tight_layout()
    plt.savefig(
        RESULTS_DIR / f"fig2h_stromal_neuronal_stacked_bar_{lvl}.pdf",
        dpi=1000,
        bbox_inches="tight",
    )
    plt.close()

# %% Fig 2i — Cell co-occurrence
# --- Co-occurrence cell types
def dice_score_abundance(u, v):
    eps = 1e-5
    u_presence = u * (u < 1) + (u >= 1)  # max(u, 1)
    v_presence = v * (v < 1) + (v >= 1)
    num = 2 * np.sum(u_presence * v_presence)
    denom = np.sum(u_presence + v_presence)
    dice = num / (denom + eps)
    return dice

def plot_correlation_matrix(
    df,
    column_list,
    method="pearson",
    min_corr_threshold=None,
    cmap="coolwarm",
    save_path=None,
    dpi=1000,
    figsize=(15, 13),
    title=None,
    annot=False,
    fmt=".2f",
    mask_diagonal=True,
    verbose=True,
):
    # Ensure list-like (Index ok, string not ok)
    if isinstance(column_list, str):
        raise TypeError(
            "column_list doit être une liste/Index de colonnes (ex: df.columns), "
            "pas une string."
        )

    # --- compute correlation
    if method == "dice":
        corr_matrix = df.loc[:, column_list].corr(method=dice_score_abundance)
        vmin, vmax = 0, 1
        cbar_label = "Dice score"
    else:
        corr_matrix = df.loc[:, column_list].corr(method=method)
        vmin, vmax = -1, 1
        cbar_label = method.capitalize()

    # Coerce numeric and handle NaNs
    corr_matrix = corr_matrix.apply(pd.to_numeric, errors="coerce")

    # Drop all-NaN rows/cols (can happen with constant columns)
    corr_matrix = corr_matrix.dropna(axis=0, how="all").dropna(axis=1, how="all")

    if verbose:
        print(f"[plot_correlation_matrix] initial shape: {corr_matrix.shape}")
        print(f"[plot_correlation_matrix] NaNs: {corr_matrix.isna().sum().sum()}")

    # --- optional filtering (but never allow empty result)
    if min_corr_threshold is not None and corr_matrix.shape[0] > 0:
        offdiag = corr_matrix.where(~np.eye(len(corr_matrix), dtype=bool))
        mask = offdiag.abs().ge(min_corr_threshold).any()

        filtered_cols = corr_matrix.columns[mask].tolist()

        if len(filtered_cols) == 0:
            if verbose:
                print(
                    f"[plot_correlation_matrix] WARNING: filtre min_corr_threshold={min_corr_threshold} "
                    "retire toutes les colonnes. On désactive le filtre."
                )
        else:
            corr_matrix = corr_matrix.loc[filtered_cols, filtered_cols]
            if verbose:
                print(f"[plot_correlation_matrix] after filter shape: {corr_matrix.shape}")

    # If still empty or 1x1, plot anyway but warn
    if corr_matrix.shape[0] == 0:
        raise ValueError(
            "La matrice de corrélation est vide (0x0). "
            "Baisse min_corr_threshold, ou vérifie que df contient des colonnes variables."
        )

    # --- hierarchical clustering (need at least 2 variables)
    if corr_matrix.shape[0] > 1:
        # Convert correlation -> distance
        dist = 1 - corr_matrix.fillna(0)
        linkage_matrix = linkage(squareform(dist, checks=False), method="ward")
        order = leaves_list(linkage_matrix)
        corr_matrix = corr_matrix.iloc[order, order]

    # --- mask diagonal if desired (often nicer)
    mask_plot = None
    if mask_diagonal and corr_matrix.shape[0] > 1:
        mask_plot = np.eye(len(corr_matrix), dtype=bool)

    # --- plot
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        corr_matrix,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        square=True,
        linewidths=0.5,
        mask=mask_plot,
        annot=annot,
        fmt=fmt,
        cbar_kws={"label": cbar_label},
    )

    if title is not None:
        ax.set_title(title, fontsize=16)

    plt.tight_layout()

    ax.set_xlabel("")
    ax.set_ylabel("")

    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()

    return corr_matrix

# %% -- Fig2i co-occurence
for lvl in ['level4']:#["level1", "level2", "level3", "level4"]:
    print(f"Stromal stack bar {lvl}")
    REFERENCE_COLUMN = f"celltype_{lvl}_curated"

    # Abundance matrix: patients × cell types
    abundance_df = (
        adata_filt.obs
        .groupby([BATCH_KEY, REFERENCE_COLUMN])
        .size()
        .unstack(fill_value=0)
    )

    # Convert to fractions (recommended)
    abundance_df = abundance_df.div(abundance_df.sum(axis=1), axis=0)

    print(f"Number of columns: {len(abundance_df.columns.tolist())}")
    topN = 50
    top_cols = abundance_df.mean(axis=0).sort_values(ascending=False).head(topN).index
    abundance_top = abundance_df[top_cols]

    plot_correlation_matrix(
        abundance_top,
        column_list=abundance_top.columns,
        method="spearman",
        min_corr_threshold=0.3,
        cmap=BLUE_ORANGE_SOFT,
        title="Spearman score of cell type co-occurrence across patients",
        save_path=RESULTS_DIR / f"fig2i_cooccurrence_matrix_{lvl}.pdf",
        dpi=1000,
        annot=False,   # True if few types
        verbose=True,  # print shape + warnings
    )

# %%
