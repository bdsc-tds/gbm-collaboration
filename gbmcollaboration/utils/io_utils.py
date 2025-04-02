"""I/O utility functions."""

from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import issparse
from env_utils import Path


def glob(path_pattern: str) -> list[Path]:
    """Perform a glob operation on a string, can be an s3 address or local path.

    Input string should contain a / and a wildcard (*) character after the last /.
    """
    if "/" not in path_pattern:
        raise ValueError(
            "Pattern should have the format: <s3 or local path to dataset>/*.zarr"
        )
    root, pattern = path_pattern.rsplit("/", 1)
    if "*" in root:
        raise ValueError(
            "Pattern should have the format: <s3 or local path to dataset>/*.zarr"
        )

    return list(Path(root).glob(pattern))


def get_object_paths(objects: str | list[str] | list[Path]):
    """Get the paths of the objects to load."""
    if isinstance(objects, str):
        return sorted(glob(objects))
    if isinstance(objects, list):
        list_path: list[str] = []
        for obj in objects:
            list_path = list_path + glob(str(obj))
        return sorted(list_path)


def read_elem_with_path(store: zarr.hierarchy.Group, path: list[str]) -> Any:
    """Recursive function to read an element from a zarr store with a given path.

    Parameters
    ----------
    store : zarr.hierarchy.Group
        The store to read from.
    path : list[str]
        The path to the element to read. E.g. ["tables","table","var","_index"]

    Returns
    -------
    Any
        The element read from the store.
    """
    if len(path) == 1:
        return read_elem(store[path[0]])
    else:
        return read_elem_with_path(store[path[0]], path[1:])


def extract_counts_matrix(store: zarr.Group, key: str) -> np.ndarray:
    """Extract the counts matrix from a zarr store.

    Parameters
    ----------
    store : zarr.Group
        The zarr store.
    key : str
        The key to the counts matrix.

    Returns
    -------
    np.ndarray
        The filtered counts matrix.
    """
    matrix = read_elem(store["tables"]["table"]["layers"][key])
    if issparse(matrix):
        matrix = matrix.toarray()
    return matrix


def extract_counts_dataframe(
    store: zarr.Group, key: str, list_valid_spot: list[str] | None = None
) -> pd.DataFrame:
    """Extract the counts dataframe from a zarr store, filtering by UMI threshold.

    Parameters
    ----------
    store : zarr.Group
        The zarr store.
    key : str
        The key to the counts matrix.
    list_valid_spot : list[str] | None
        List of valid spots to filter the deconvolution dataframe.

    Returns
    -------
    pd.DataFrame
        The filtered counts dataframe.
    """
    matrix = extract_counts_matrix(store=store, key=key)
    obs_names = read_elem(store["tables"]["table"]["obs"]["_index"])
    var_names = read_elem(store["tables"]["table"]["var"]["_index"])

    if list_valid_spot:  # Filter by UMI threshold
        valid_spot_mask = [k in list_valid_spot for k in obs_names]
        filtered_obs_names = obs_names[valid_spot_mask]
        filtered_matrix = matrix[valid_spot_mask]
        return pd.DataFrame(
            filtered_matrix, index=filtered_obs_names, columns=var_names
        )

    return pd.DataFrame(matrix, index=obs_names, columns=var_names)


def extract_deconvolution_dataframe(
    store: zarr.Group,
    refactor_dict: dict | None = None,
    list_valid_spot: list[str] | None = None,
    deconvolution_slot_name: str = "deconvolution",
    keep_original_celltypes: bool = False,
) -> list[pd.DataFrame]:
    """Extract the deconvolution dataframe from a zarr store, filtering by UMI threshold.

    Parameters
    ----------
    store : zarr.Group
        The zarr store.
    refactor_dict : dict | None
        Optional dictionary to refactor the cell types.
    list_valid_spot : list[str] | None
        List of valid spots to filter the deconvolution dataframe.
    deconvolution_slot_name : str
        The name of the deconvolution wanted when several are available.
        By default, "deconvolution" is used, as it was the original name in
        Mosaic visium pipeline.
    keep_original_celltypes : bool
        Whether to keep the original cell types when refactoring. Default is False.

    Returns
    -------
    list[pd.DataFrame]
        The filtered deconvolution dataframe and the fractions calculated from it.
    """
    deconv_df = read_elem(store["tables"]["table"]["obsm"][deconvolution_slot_name])
    if list_valid_spot:  # Filter by UMI threshold
        obs_names = read_elem(store["tables"]["table"]["obs"]["_index"])
        valid_spot_mask = [k in list_valid_spot for k in obs_names]
        deconv_df = deconv_df[valid_spot_mask]

    if refactor_dict is not None:  # Refactor the cell types
        for new_type, old_types in refactor_dict.items():
            deconv_df[new_type] = deconv_df[old_types].sum(axis=1)
        if not keep_original_celltypes:
            deconv_df.drop(
                columns=[
                    old
                    for old_list in refactor_dict.values()
                    for old in old_list
                    if old not in refactor_dict.keys()
                ],
                inplace=True,
            )
    deconv_frac_df = deconv_df.div(deconv_df.sum(axis=1), axis=0)
    return [deconv_df, deconv_frac_df]


def extract_coords_dataframe(store: zarr.Group) -> pd.DataFrame:
    """Extract the coordinates dataframe from a zarr store.

    Parameters
    ----------
    store : zarr.Group
        The zarr store.

    Returns
    -------
    pd.DataFrame
        The coordinates dataframe.
    """
    return pd.DataFrame(
        {
            "imagerow": read_elem(store["tables"]["table"]["obs"]["imagerow"]),
            "imagecol": read_elem(store["tables"]["table"]["obs"]["imagecol"]),
        },
        index=read_elem(store["tables"]["table"]["obs"]["_index"]),
    )


def check_path(
    path: str | list[str | Path],
):
    """Check if any file matches the Path pattern.

    Parameters
    ----------
    path : str | list[str | Path]
        Path to data

    Raises
    ------
    ValueError
        If no files are found matching the path pattern.
    """
    if isinstance(path, list):
        for p in path:
            # Convert Path to string if necessary
            p_str = str(p)
            matching_files = glob(p_str)
            if not matching_files:
                raise ValueError(f"No files found matching the path pattern: {p_str}")
    else:
        # Convert Path to string if necessary
        path_str = str(path) if isinstance(path, Path) else path
        matching_files = glob(path_str)
        if not matching_files:
            raise ValueError(f"No files found matching the path pattern: {path_str}")


def check_duplicated_block(
    samples: list[str],
):
    """Check that there is no duplicated blocks in the dataset, corresponding to the same sample.

    Parameters
    ----------
    samples : list[str]
        List of samples
    """
    short_sample_id = [
        i.replace("_vis", "").replace("_vbu", "").replace("_bu3", "") for i in samples
    ]
    counter = Counter(short_sample_id)
    duplicates = [item for item, count in counter.items() if count > 1]
    if duplicates:
        print(
            f"Warning: duplicate sample found: {duplicates}.\nPlease check which one has better QC (between *_vis *_vbu *_bu3) and\nremove the other using ignore_sample_ids argument of the VisiumCohort class."
        )
