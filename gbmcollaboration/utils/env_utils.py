"""Utilities to check the environment."""

import os
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Type, Union, List, Sequence, Dict
from warnings import warn

from cloudpathlib import S3Path


def filter_by_depth(
    bucket_name: str,
    prefix: str | None,
    s3_client: Any,  # noqa
    max_depth: int = 3,
) -> List[str]:
    """Filter object keys and prefixes by depth in a S3 bucket.

    Parameters
    ----------
    bucket_name : str
        Name of the S3 bucket.
    prefix : str
        Prefix of the objects' keys in the bucket. If ``None``, the ``inspect_bucket`` function
        will try to list objects at the root of the S3 bucket. This may raise an error if you are
        not authorized to do so. If ``prefix`` is a string, it must end with ``/``.
    s3_client : botocore.client.S3
        Client used to list objects.
    max_depth : int, default=3
        Maximum depth. Note that ``max_depth=1`` corresponds to the root of the S3 bucket.

    Returns
    -------
    List[str]
        Prefixes and object keys whose depth is less than ``max_depth``.
    """
    if prefix is not None:
        current_depth = 1 + prefix.count("/")
        if current_depth == max_depth:
            return [prefix]

    response = s3_client.list_objects(
        Bucket=bucket_name,
        Prefix=prefix if prefix is not None else "",
        Delimiter="/",
        MaxKeys=1_000,
    )
    is_truncated = response.get("IsTruncated")
    if is_truncated:
        warn(f"Prefix {prefix} contains more than 1 000 keys/objects. The list has been truncated!")
    contents = response.get("Contents")
    if contents is None:
        keys = []
    else:
        _all_keys = [obj.get("Key") for obj in contents]
        keys = list(filter(lambda k: (not k.endswith("/")), _all_keys))
    common_prefixes = response.get("CommonPrefixes")
    if common_prefixes is not None:
        prefixes = [obj.get("Prefix") for obj in common_prefixes]
        other_prefixes = [
            filter_by_depth(
                bucket_name=bucket_name,
                prefix=(prefix + f"{this_prefix.split('/')[-2]}/")
                if prefix is not None
                else f"{this_prefix.split('/')[-2]}/",
                max_depth=max_depth,
                s3_client=s3_client,
            )
            for this_prefix in prefixes
        ]
        return keys + [elt for sub_lst in other_prefixes for elt in sub_lst]
    else:
        return keys


def get_manifest() -> Dict[str, Any]:
    """
    Get the manifest file.

    Note: this only works from inside SageMaker.

    Returns
    -------
    Dict[str, Any]
        The manifest file as a dictionary.

    Warns
    -----
    UserWarning
        If manifest is not found, returns an empty dictionary.
    """
    manifest_path = Path("/home/sagemaker-user/.config/abstra/project.json")
    if manifest_path.exists():
        with open(manifest_path, "r") as f:
            return json.load(f)
    else:
        warn("Manifest not found, returning empty manifest")
        return {}


def get_s3_path(dataset_id: str, path: str | None = None):
    """Build the full s3 path.

    Parameters
    ----------
    dataset_id : str
        Dataset id
    path : str | None
        Path in the dataset
    """
    manifest = get_manifest()
    if dataset_id not in manifest.get("datasets", {}):
        raise ValueError(f"Dataset id {dataset_id} not found in manifest")

    path = path or ""
    s3_info = manifest["datasets"][dataset_id]
    s3_base_path = S3Path("s3://") / s3_info["bucket"] / s3_info["prefix"]
    return s3_base_path / path
