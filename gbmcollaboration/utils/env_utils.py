"""Utilities to check the environment."""

import os
from pathlib import Path as LocalPath
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Type, Union, Dict
import json
from warnings import warn
from cloudpathlib import S3Path

PathType = Union[LocalPath, S3Path]

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
    manifest_path = LocalPath("/home/sagemaker-user/.config/abstra/project.json")
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


class Path:
    """
    A custom Path class that returns either an S3Path or a LocalPath based on the input.

    This class acts as a factory, creating the appropriate path object depending on whether
    the input path is an S3 URI or a local file system path.

    Examples
    --------
    >>> path = Path("s3://my-bucket/my-file.txt")
    >>> isinstance(path, S3Path)
    True

    >>> path = Path("/home/user/my-file.txt")
    >>> isinstance(path, LocalPath)
    True
    """

    # HACK : the methods of Path are not correctly identified by mypy,
    # see https://github.com/python/mypy/issues/1021
    # list of methods taken from https://github.com/drivendataorg/cloudpathlib/blob/master/README.md
    if TYPE_CHECKING:
        __truediv__: ClassVar[Callable]
        absolute: ClassVar[Callable]
        anchor: ClassVar[Callable]
        as_uri: ClassVar[Callable]
        drive: ClassVar[Callable]
        exists: ClassVar[Callable]
        glob: ClassVar[Callable]
        is_absolute: ClassVar[Callable]
        is_dir: ClassVar[Callable]
        is_file: ClassVar[Callable]
        is_relative_to: ClassVar[Callable]
        iterdir: ClassVar[Callable]
        joinpath: ClassVar[Callable]
        mkdir: ClassVar[Callable]
        name: ClassVar[Callable]
        parent: ClassVar[Callable]
        parents: ClassVar[Callable]
        parts: ClassVar[Callable]
        read_bytes: ClassVar[Callable]
        read_text: ClassVar[Callable]
        relative_to: ClassVar[Callable]
        rename: ClassVar[Callable]
        replace: ClassVar[Callable]
        resolve: ClassVar[Callable]
        rglob: ClassVar[Callable]
        rmdir: ClassVar[Callable]
        samefile: ClassVar[Callable]
        stat: ClassVar[Callable]
        stem: ClassVar[Callable]
        suffix: ClassVar[Callable]
        suffixes: ClassVar[Callable]
        touch: ClassVar[Callable]
        unlink: ClassVar[Callable]
        with_name: ClassVar[Callable]
        with_stem: ClassVar[Callable]
        with_suffix: ClassVar[Callable]
        write_bytes: ClassVar[Callable]
        write_text: ClassVar[Callable]
        as_url: ClassVar[Callable]
        clear_cache: ClassVar[Callable]
        cloud_prefix: ClassVar[Callable]
        copy: ClassVar[Callable]
        copytree: ClassVar[Callable]
        download_to: ClassVar[Callable]
        etag: ClassVar[Callable]
        fspath: ClassVar[Callable]
        is_junction: ClassVar[Callable]
        is_valid_cloudpath: ClassVar[Callable]
        rmtree: ClassVar[Callable]
        upload_from: ClassVar[Callable]
        validate: ClassVar[Callable]
        walk: ClassVar[Callable]
        with_segments: ClassVar[Callable]
        bucket: ClassVar[Callable]
        key: ClassVar[Callable]

    def __new__(  # type: ignore
        cls, path: Union[str, LocalPath, S3Path], *args: Any, **kwargs: Any
    ) -> Union[LocalPath, S3Path]:
        """
        Create and return a new instance of either ``S3Path`` or ``LocalPath``.

        This method determines whether the given path is an S3 path or a local path
        and returns the appropriate instance.

        Parameters
        ----------
        path : Union[str, LocalPath, S3Path]
            The path to be converted into a Path object.
        *args : Any
            Additional positional arguments to be passed to the path constructor.
        **kwargs : Any
            Additional keyword arguments to be passed to the path constructor.

        Returns
        -------
        Union[LocalPath, S3Path]
            An instance of either LocalPath or S3Path, depending on the input path.
        """
        if isinstance(path, (LocalPath, S3Path)):
            return path
        if str(path).startswith("s3://"):
            return S3Path(path, *args, **kwargs)
        else:
            return LocalPath(path)

    @classmethod
    def __class_getitem__(cls, *args: Any, **kwargs: Any) -> Type["Path"]:
        """
        Support for type hinting with subscripts (e.g., Path[str]).

        This method allows the Path class to be used with type subscripts without
        changing its behavior. It's primarily used for static type checking and
        doesn't affect runtime behavior.

        Returns
        -------
        Type[pathlib.Path]
            The Path class itself.
        """
        return cls


ABSTRA_ROOT_PATH = Path("/home/owkin")
SAGEMAKER_ROOT_PATH = Path("/home/sagemaker-user")
SAGEMAKER_JOB_ROOT_PATH = Path("/opt/ml")


def am_i_on_sagemaker() -> bool:
    """
    Check if the user is on SageMaker.

    Returns
    -------
    bool
        True if the user is on SageMaker, False otherwise.
    """
    return SAGEMAKER_ROOT_PATH.exists() or SAGEMAKER_JOB_ROOT_PATH.exists()


def am_i_on_abstra() -> bool:
    """
    Check if the user is on Abstra.

    Returns
    -------
    bool
        True if the user is on Abstra, False otherwise.
    """
    return ABSTRA_ROOT_PATH.exists()


def am_i_on_sagemaker_training_job():
    """Check if the code is running on a SageMaker training job.

    Returns
    -------
    bool
        ``True`` if running on a SageMaker training job, ``False`` otherwise.
    """
    if os.environ.get("SM_NUM_CPUS") is not None:
        return True

    return False


def get_project_root_path() -> Path:
    """Get the project root path.

    Returns
    -------
    pathlib.Path
        The project root path.

    Raises
    ------
    RuntimeError
        If the environment is unknown.
    """
    if am_i_on_abstra():
        return ABSTRA_ROOT_PATH
    elif am_i_on_sagemaker():
        return SAGEMAKER_ROOT_PATH
    else:
        raise RuntimeError("Unknown environment")
