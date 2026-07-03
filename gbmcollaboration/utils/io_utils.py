"""I/O utility functions."""

from __future__ import annotations

from collections import Counter, namedtuple
from operator import attrgetter
from typing import Any

import anndata
import h5py
import boto3
import numpy as np
import pandas as pd
import fsspec
import s3fs
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import issparse
from .env_utils import Path


def s3_open_h5(s3_path, backend="fsspec"):
    if backend == "fsspec":
        s3f = fsspec.open(s3_path, mode="rb").open()
        h5f = h5py.File(s3f, mode="r")
    elif backend == "s3fs":
        fs = s3fs.S3FileSystem(anon=False)
        s3f = fs.open(s3_path, "rb")
        h5f = h5py.File(s3f, mode="r")
    return h5f


def _read_anndata_dict(f, field, field_name):
    if isinstance(field, bool) and field is True:
        field_dict = read_elem(f[field_name])
    else:
        field_dict = {}
        for k in field:
            field_dict[k] = read_elem(f[f"{field_name}/{k}"])
    return field_dict

def s3_read_h5ad_partial(
    s3_path,
    X=None,
    obs=True,
    var=True,
    uns=[],
    obsm=[],
    varm=[],
    layers=[],
    raw=False,
):


    h5f = s3_open_h5(s3_path)
    if raw:
        h5f = h5f["raw"]

    if X:
        X = read_elem(h5f["X"])
    if obs and "obs" in h5f:
        obs = read_elem(h5f["obs"])
    if var and "var" in h5f:
        var = read_elem(h5f["var"])
    if uns and "uns" in h5f:
        uns = _read_anndata_dict(h5f, uns, "uns")
    if obsm and "obsm" in h5f:
        obsm = _read_anndata_dict(h5f, obsm, "obsm")
    if varm and "varm" in h5f:
        varm = _read_anndata_dict(h5f, varm, "varm")
    if layers and "layers" in h5f:
        layers = _read_anndata_dict(h5f, layers, "layers")

    adata = anndata.AnnData(
        X=X,
        obs=obs,
        var=var,
        uns=uns,
        obsm=obsm,
        varm=varm,
        layers=layers,
    )
    return adata


def s3_ls(
    bucket,
    path="",
    start=None,
    end=None,
    recursive=True,
    list_dirs=True,
    list_objs=True,
    limit=None,
):
    """
    https://stackoverflow.com/a/51372405
    Iterator that lists a bucket's objects under path, (optionally) starting with
    start and ending before end.
    If recursive is False, then list only the "depth=0" items (dirs and objects).
    If recursive is True, then list recursively all objects (no dirs).
    Args:
        bucket:
            a str or boto3.resource('s3').Bucket().
        path:
            a directory in the bucket.
        start:
            optional: start key, inclusive (may be a relative path under path, or
            absolute in the bucket)
        end:
            optional: stop key, exclusive (may be a relative path under path, or
            absolute in the bucket)
        recursive:
            optional, default True. If True, lists only objects. If False, lists
            only depth 0 "directories" and objects.
        list_dirs:
            optional, default True. Has no effect in recursive listing. On
            non-recursive listing, if False, then directories are omitted.
        list_objs:
            optional, default True. If False, then directories are omitted.
        limit:
            optional. If specified, then lists at most this many items.
    Returns:
        an iterator of S3Obj.
    Examples:
        # set up
        >>> s3 = boto3.resource('s3')
        ... bucket = s3.Bucket('bucket-name')
        # iterate through all S3 objects under some dir
        >>> for p in s3list(bucket, 'some/dir'):
        ...     print(p)
        # iterate through up to 20 S3 objects under some dir, starting with foo_0010
        >>> for p in s3list(bucket, 'some/dir', limit=20, start='foo_0010'):
        ...     print(p)
        # non-recursive listing under some dir:
        >>> for p in s3list(bucket, 'some/dir', recursive=False):
        ...     print(p)
        # non-recursive listing under some dir, listing only dirs:
        >>> for p in s3list(bucket, 'some/dir', recursive=False, list_objs=False):
        ...     print(p)
    """
    s3 = boto3.resource("s3")

    if isinstance(bucket, str):
        bucket = bucket = s3.Bucket(bucket)

    S3Obj = namedtuple("S3Obj", ["key", "mtime", "size", "ETag"])

    kwargs = dict()
    if start is not None:
        if not start.startswith(path):
            start = os.path.join(path, start)
        # note: need to use a string just smaller than start, because
        # the list_object API specifies that start is excluded (the first
        # result is *after* start).
        kwargs.update(Marker=_prev_str(start))
    if end is not None:
        if not end.startswith(path):
            end = os.path.join(path, end)
    if not recursive:
        kwargs.update(Delimiter="/")
        if not path.endswith("/"):
            path += "/"
    kwargs.update(Prefix=path)
    if limit is not None:
        kwargs.update(PaginationConfig={"MaxItems": limit})

    paginator = bucket.meta.client.get_paginator("list_objects")
    for resp in paginator.paginate(Bucket=bucket.name, **kwargs):
        q = []
        if "CommonPrefixes" in resp and list_dirs:
            q = [S3Obj(f["Prefix"], None, None, None) for f in resp["CommonPrefixes"]]
        if "Contents" in resp and list_objs:
            q += [
                S3Obj(f["Key"], f["LastModified"], f["Size"], f["ETag"])
                for f in resp["Contents"]
            ]
        # note: even with sorted lists, it is faster to sort(a+b)
        # than heapq.merge(a, b) at least up to 10K elements in each list
        q = sorted(q, key=attrgetter("key"))
        if limit is not None:
            q = q[:limit]
            limit -= len(q)
        for p in q:
            if end is not None and p.key >= end:
                return
            yield p


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
    # keep only the first 9 characters
    short_sample_id = [
        i[:9] for i in samples
    ]
    counter = Counter(short_sample_id)
    duplicates = [item for item, count in counter.items() if count > 1]
    if duplicates:
        print(
            f"Warning: duplicate sample found: {duplicates}.\nPlease check which one has better QC (between *_vis *_vbu *_bu3) and\nremove the other using ignore_sample_ids argument of the VisiumCohort class."
        )
