"""Data processing utility functions."""

from __future__ import annotations

import warnings
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import coo_array, coo_matrix
from sklearn.linear_model import LinearRegression
from tqdm import tqdm


def get_common_list(
    lists_of_elem: list[list[str]],
    how: str = "union",
) -> list[str]:
    """Get the common elements in multiple lists."""
    if how == "intersection":
        return sorted(set().intersection(*map(set, lists_of_elem)))
    elif how == "union":
        return sorted(set().union(*map(set, lists_of_elem)))
    else:
        raise ValueError("how must be either 'intersection' or 'union'")


def normalize_coords(coords: pd.DataFrame) -> np.ndarray:
    """Normalize the coordinates of the spots.

    In practice, the spots are arranged on a hexagonal grid. In the output of space
    ranger we get two sets of coordinates: the image coordinates (hexagonal grid, the
    'true' coordinates) and the array coordinates (rectangular grid, the 'numerotation'
    of the spots). This function normalizes the image coordinates to the array coords.
    To do so, we fit a linear regression between the two sets of coordinates and divide
    the image coordinates by the slope of the regression.

    Parameters
    ----------
    coords : pd.DataFrame
        The coordinates of the spots. Must contain the columns "imagerow", "imagecol",
        "array_col", and "array_row".

    Returns
    -------
    np.ndarray
        The normalized coordinates.
    """
    img_row = coords["imagerow"].to_numpy()
    img_col = coords["imagecol"].to_numpy()
    array_col = coords["array_col"].to_numpy()
    array_row = coords["array_row"].to_numpy()
    corr_coef = np.abs(np.corrcoef(array_row.flatten(), img_row.flatten())[0, 1])

    # If corr_coef < 0.9, we consider arrays and images columns inverted
    if corr_coef < 0.9:
        reg_col = LinearRegression().fit(array_col.reshape(-1, 1), img_row)
        reg_row = LinearRegression().fit(array_row.reshape(-1, 1), img_col)
        img_coords = np.vstack([img_col, img_row]).T
    else:
        reg_row = LinearRegression().fit(array_row.reshape(-1, 1), img_row)
        reg_col = LinearRegression().fit(array_col.reshape(-1, 1), img_col)
        img_coords = np.vstack([img_row, img_col]).T

    # Columns is the primary dimension (i.e. the flat axis of the hexagon)
    if (array_col.max() - array_col.min()) > (array_row.max() - array_row.min()):
        unit = 2 * reg_col.coef_
    # Row is the primary dimension
    else:
        unit = 2 * reg_row.coef_
    # We also need to multiply by 2 because array_row/array_col steps by 2 in each
    # direction

    return img_coords / np.abs(unit)


def compute_neighbors(
    coords: np.ndarray,
    radius: int,
    epsilon: float = 0.1,
) -> pd.DataFrame:
    """Compute neihborhoods.

    Parameters
    ----------
    coords : np.ndarray
        The normalized coordinates of the spots. Distance between spots ~1.
    radius : int
        The radius of the neighborhood. The neighborhood is a hexagon of side length
        `radius`.
    epsilon : float
        The tolerance for the distance.

    Returns
    -------
    pd.DataFrame
        The neighborhood matrix. Each row corresponds to a spot and contains the
        indices of the neighbors of the spot. If a spot has less than x neighbors, the
        remaining entries are filled with NaN.
        When radius >=2, the neighbors are filled in a way that each layer of neighbors
        is filled in the adequate range of the distance matrix.
    """
    radius_ = radius + epsilon
    dist_mtx = scipy.spatial.distance.cdist(coords, coords)
    dist_mtx += np.eye(len(coords)) * (radius_ + 1)

    n_spot = len(coords)
    n_neighbors = 6 * radius * (radius + 1) // 2
    out = np.zeros((n_spot, n_neighbors))
    out.fill(np.nan)
    for j in range(0, radius):
        if j == 1:
            # specifically for the second ring of neighbors, we distinguish between
            # the neighbors at a distance of sqrt(3) and the neighbors at a distance
            # of 2. This is necessary for an application in cellcellcom.
            neighbors_j = np.where(
                (dist_mtx <= (np.sqrt(3) + epsilon)) & (dist_mtx > (1 + epsilon))
            )
            for i in range(n_spot):
                to_fill = neighbors_j[1][neighbors_j[0] == i]
                out[i, 6 : 6 + len(to_fill)] = to_fill
            neighbors_j = np.where(
                (dist_mtx <= (2 + epsilon)) & (dist_mtx > (np.sqrt(3) + epsilon))
            )
            for i in range(n_spot):
                to_fill = neighbors_j[1][neighbors_j[0] == i]
                out[i, 12 : 12 + len(to_fill)] = to_fill

        else:
            neighbors_j = np.where(
                (dist_mtx <= (j + 1 + epsilon)) & (dist_mtx > (j + epsilon))
            )
            n_neighbors_j_minus_1 = 6 * (j + 1) * (j) // 2
            for i in range(n_spot):
                to_fill = neighbors_j[1][neighbors_j[0] == i]
                out[i, n_neighbors_j_minus_1 : n_neighbors_j_minus_1 + len(to_fill)] = (
                    to_fill
                )
    return pd.DataFrame(out, index=np.arange(n_spot))


def compute_adjacency_matrix(
    neighbors: tuple[np.ndarray[Any, np.dtype[np.signedinteger[Any]]], ...],
) -> coo_matrix:
    """Compute the adjacency matrix."""
    data = np.ones(len(neighbors[0]))
    return coo_array((data, (neighbors[0], neighbors[1])))


def convert_obsm_to_adata(adata: ad.AnnData, obsm_key: str) -> ad.AnnData:
    """Convert obsm matrices to adata.

    This function converts an obsm matrix to an AnnData object in order to use it with liana functions.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the obsm matrix.
    obsm_key : str
        Key of the obsm matrix to convert to adata.

    Returns
    -------
    ad.AnnData
        AnnData object with obsm matrix as X.

    Raises
    ------
    KeyError
        If obsm_key is not in adata.obsm.
    """
    if obsm_key not in adata.obsm.keys():
        raise KeyError(f"No key with name '{obsm_key}' in adata.obsm")
    obsm = adata.obsm[obsm_key]
    obsm_df = pd.DataFrame(obsm, index=adata.obs.index)
    obsm_adata = ad.AnnData(
        obsm_df,
        # obs=adata.obs,
        var=pd.DataFrame(index=obsm_df.columns),
        uns=adata.uns,
        obsm=adata.obsm,
        obsp=adata.obsp,
    )
    return obsm_adata


def get_subset_gene_list_index(
    gene_list: np.ndarray, gene_sublist: np.ndarray
) -> np.ndarray:
    """Return position of gene_sublist's genes in gene_list."""
    out = np.array([np.where(gene_list == x) for x in gene_sublist]).squeeze()

    return out


def compute_neighbors_as_multiindex(
    coords: pd.DataFrame, radius: int, sample_id: str
) -> pd.DataFrame:
    """Compute the neighbors of the spots."""
    out = compute_neighbors(coords.values, radius=radius)
    out["sample_id"] = [sample_id] * len(out)
    out.index = coords.index
    out = out.set_index("sample_id", append=True)
    out = out.reorder_levels(["sample_id", "x", "y"])
    return out


def load_coords_as_multiindex(
    store: zarr.hierarchy.Group,
    sample_id: str,
) -> pd.DataFrame:
    """Load the coordinates from a zarr store as a multiindex dataframe."""
    coords_to_load = ["imagerow", "imagecol", "array_col", "array_row"]
    out = pd.DataFrame(
        {
            key: read_elem(store["tables"]["table"]["obs"][key])
            for key in coords_to_load
        },
        index=read_elem(store["tables"]["table"]["obs"]["_index"]),
    ).reset_index(names="barcode")
    out["x"] = out["imagerow"]
    out["y"] = out["imagecol"]
    out["sample_id"] = [sample_id] * out.shape[0]
    out = out.set_index(["sample_id", "x", "y"])

    norm_coords = normalize_coords(out)
    out["norm_x"], out["norm_y"] = norm_coords[:, 0], norm_coords[:, 1]
    return out


def compute_neighbors_anndata(
    adata: ad.AnnData,
    radius: int,
    sample_id: str,
    include_spot_itself: bool = False,
) -> pd.DataFrame:
    """Compute the neighbors of the spots in an AnnData object.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the coordinates of the spots.
    radius : int
        The radius of the neighborhood.
    sample_id : str
        The id of the sample.
    include_spot_itself : bool
        Whether to include the spot itself in the neighborhood as the first neighbor.

    Returns
    -------
    pd.DataFrame
        The neighbors of the spots.
    """
    if (
        "array_row_norm" not in adata.obs.columns
        or "array_col_norm" not in adata.obs.columns
    ):
        warnings.warn(
            "array_row_norm or array_col_norm not in adata.obs.columns, "
            "will fallback to imagerow and imagecol",
            stacklevel=2,
        )
        df_coords = adata.obs.loc[:, ["imagerow", "imagecol", "array_row", "array_col"]]
        df_coords["x"] = df_coords["imagerow"].values
        df_coords["y"] = df_coords["imagecol"].values
        df_coords["sample_id"] = [sample_id] * df_coords.shape[0]
        df_coords = df_coords.set_index(["sample_id", "x", "y"])

        norm_coords = normalize_coords(df_coords)
        df_coords["norm_x"], df_coords["norm_y"] = norm_coords[:, 0], norm_coords[:, 1]

        neighbors_idx = compute_neighbors_as_multiindex(
            df_coords[["norm_x", "norm_y"]].loc[sample_id],
            radius=radius,
            sample_id=sample_id,
        )
        return neighbors_idx

    df_coords = adata.obs.loc[:, ["array_row_norm", "array_col_norm"]]
    df_coords["norm_x"] = df_coords["array_row_norm"].values
    df_coords["norm_y"] = df_coords["array_col_norm"].values
    df_coords["x"] = df_coords["array_row_norm"].values
    df_coords["y"] = df_coords["array_col_norm"].values
    df_coords["sample_id"] = [sample_id] * df_coords.shape[0]
    df_coords = df_coords.set_index(["sample_id", "x", "y"])

    neighbors_idx = compute_neighbors_as_multiindex(
        df_coords[["array_row_norm", "array_col_norm"]].loc[sample_id],
        radius=radius,
        sample_id=sample_id,
    )
    if include_spot_itself:
        neighbors_df = pd.DataFrame(index=neighbors_idx.index)
        for col in range(len(neighbors_idx.columns)):
            neighbors_df[col + 1] = neighbors_idx[col]
        neighbors_df[0] = np.arange(len(neighbors_idx))
        neighbors_df = neighbors_df.sort_index(axis=1)
        neighbors_idx = neighbors_df

    return neighbors_idx


def compute_neighbors_anndata_cohort(
    adata: ad.AnnData,
    radius: int,
    column_sample_id: str,
    include_spot_itself: bool = False,
) -> pd.DataFrame:
    """Compute the neighbors of the spots in an AnnData object for a cohort.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object containing the coordinates of the spots.
    radius : int
        The radius of the neighborhood.
    column_sample_id : str
        The column of the adata.obs dataframe containing the sample ids.
    include_spot_itself : bool
        Whether to include the spot itself in the neighborhood as the first neighbor.

    Returns
    -------
    pd.DataFrame
        The neighbors of the spots.
    """
    list_sample_ids = adata.obs[column_sample_id].unique()
    neighbors_idx = []
    for sample_id in tqdm(list_sample_ids):
        adata_sample = adata[adata.obs[column_sample_id] == sample_id]
        neighbors_idx.append(
            compute_neighbors_anndata(
                adata_sample, radius, sample_id, include_spot_itself
            )
        )
    return pd.concat(neighbors_idx)
