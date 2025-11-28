"""Visium Cohort Dataset."""

from __future__ import annotations

import copy
import warnings
from typing import Any, Literal

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as sd
from joblib import Parallel, delayed
from loguru import logger

from gbmcollaboration.data.dataset_classes.visium_sample import VisiumSample
from gbmcollaboration.utils import (
    check_duplicated_block,
    check_path,
    get_common_list,
    get_object_paths
)
from cloudpathlib import S3Path


class VisiumCohort:
    """Dataset class integrating multiple samples' zarr store's features.

    Args
    ----
    path_dataset : str | list[str | S3Path]
        Paths to the zarr stores or a glob pattern.
        '/path/to/dataset/*.zarr'.
    normalization : str
        The normalization / normalization to use for the dataset.
    radius : int | None
        If not None, compute the neighbors of the spots with the given radius.
    deconvolution_slot_name : str
        The name of the deconvolution wanted when several are available.
        By default, "deconv_cohort_level_2" is used.
    list_sample_ids : list[str] | None
        List of sample ids to load. If None, load all samples.
    ignore_sample_ids : list[str] | None
        List of sample ids to ignore.
        Useful when you want to load all samples but avoid errors due to QC issues.
    gene_threshold : int
        Minimum number of genes to keep a sample.
    gene_combination_method : Literal["intersection", "union"]
        How to match the genes between samples, should be either
            - "intersection" : only keep the genes present in all samples
            - "union" : keep all genes present in any sample
    """

    def __init__(
        self,
        path_dataset: str | list[str | S3Path],
        normalization: str,
        radius: int | None = None,
        deconvolution_slot_name: str = "deconv_cohort_level_2",
        list_sample_ids: list[str] | None = None,
        ignore_sample_ids: list[str] | None = None,
        gene_threshold: int = 0,
        gene_combination_method: Literal["intersection", "union"] = "union",
    ):
        # Check if provided path match any files
        check_path(path_dataset)

        # Fill cohort level attributes
        self._path_dataset = path_dataset
        self._normalization = normalization
        self._radius = radius
        self._deconvolution_slot_name = deconvolution_slot_name
        self._gene_threshold = gene_threshold
        list_path_samples = get_object_paths(self._path_dataset)
        self._list_path_samples = list_path_samples

        if list_sample_ids is None:
            list_sample_ids = [S3Path(path).stem for path in list_path_samples]
        if ignore_sample_ids is not None:
            list_sample_ids = [
                sample_id
                for sample_id in list_sample_ids
                if not any(ign in sample_id for ign in ignore_sample_ids)
            ]

        # Check that there is no duplicated block in the dataset, with sampleX_vis and sampleX_vbu for example
        check_duplicated_block(list_sample_ids)

        # get sample level attributes and compute cohort level attributes.
        list_gene_names = []
        list_celltypes = []
        self._dict_visium_samples: dict[str, VisiumSample] = {}
        self._coords = pd.DataFrame()
        self.neighbors_idx: dict[str, Any] = {}
        self._list_sample_ids = []

        sample_list = Parallel(n_jobs=-1, prefer="threads")(
            delayed(self.initialize_sample)(sample_id) for sample_id in list_sample_ids
        )

        sample_list = [
            sample_id_sample
            for sample_id_sample in sample_list
            if sample_id_sample[1] is not None
        ]

        for sample_id_sample in sample_list:
            sample_id = sample_id_sample[0]
            sample = sample_id_sample[1]
            if len(sample.list_gene_names) < self._gene_threshold:
                logger.warning(
                    f"Skipping sample {sample_id} with {len(sample.list_gene_names)} genes"
                )
                continue
            self._list_sample_ids.append(sample_id)
            self._dict_visium_samples[sample_id] = sample
            list_gene_names.append(sample.list_gene_names)
            list_celltypes.append(sample.list_celltypes)
            self._coords = pd.concat([self._coords, sample._coords])
            if radius is not None:
                self.neighbors_idx[sample_id] = sample.neighbors_idx

        self._common_genes = get_common_list(
            list_gene_names, how=gene_combination_method
        )
        self._common_celltypes = get_common_list(list_celltypes)

        self._loaded_counts = False
        self._loaded_deconvolution = False

    def initialize_sample(self, sample_id: str):
        """Initialize a sample."""
        try:
            path_sample = [
                str(path) for path in self._list_path_samples if sample_id in str(path)
            ][0]
            sample = VisiumSample(
                path_sample_zarr=path_sample,
                sample_id=sample_id,
                normalization=self._normalization,
                radius=self._radius,
                deconvolution_slot_name=self._deconvolution_slot_name,
            )
            return sample_id, sample
        except Exception as e:
            logger.error(
                f"Error while initializing sample {sample_id}: {e}. The sample is removed from the cohort."
            )
            return sample_id, None

    @property
    def common_genes(self) -> list[str]:
        """Return the common genes."""
        return self._common_genes

    @property
    def common_celltypes(self) -> list[str]:
        """Return the common celltypes."""
        return self._common_celltypes

    @property
    def list_sample_ids(self) -> list[str]:
        """Return the list of samples."""
        return self._list_sample_ids

    @property
    def radius(self) -> int | None:
        """Return the radius chosen for the cohort."""
        return self._radius

    @property
    def coords(self) -> pd.DataFrame:
        """Return the coordinates of the spots."""
        return self._coords

    @property
    def loaded_counts(self) -> bool:
        """Return True if the counts are loaded."""
        return self._loaded_counts

    @property
    def loaded_deconvolution(self) -> bool:
        """Return True if the deconvolution is loaded."""
        return self._loaded_deconvolution

    def load_counts_sample(self, sample_id: str) -> tuple[str, np.ndarray]:
        """Load the counts data for a given sample."""
        counts = self._dict_visium_samples[sample_id].load_counts(self._common_genes)
        return sample_id, counts

    def load_counts_cohort(self, reload: bool = False) -> dict[str, Any]:
        """Load the counts data for all samples in the cohort."""
        if self._loaded_counts and not reload:
            counts_dict = {
                sample_id: self._dict_visium_samples[sample_id].counts
                for sample_id in self.list_sample_ids
            }
            return counts_dict

        # Parallel version
        counts_list = Parallel(n_jobs=-1, prefer="threads")(
            delayed(self.load_counts_sample)(sample_id)
            for sample_id in self.list_sample_ids
        )
        # Necessary step in parallel version
        # Otherwise the sample's attribute are not modified
        # including the loaded_counts, list_gene_names and counts attributes
        for sample_id_counts in counts_list:
            id = sample_id_counts[0]
            counts = sample_id_counts[1]
            self._dict_visium_samples[id]._loaded_counts = True
            self._dict_visium_samples[id]._counts = counts
            if self._common_genes is not None:
                self._dict_visium_samples[id]._list_gene_names = self._common_genes
        counts_dict = dict(counts_list)

        self._loaded_counts = True

        return counts_dict

    def load_deconv_sample(
        self, sample_id: str, celltype_mapping_dict, keep_original_celltypes, reload
    ) -> tuple[str, np.ndarray, np.ndarray]:
        """Load the deconvolution data for a given sample."""
        deconv, deconv_frac = self._dict_visium_samples[sample_id].load_deconvolution(
            list_common_cell_types=self._common_celltypes,
            celltype_mapping_dict=celltype_mapping_dict,
            keep_original_celltypes=keep_original_celltypes,
            reload=reload,
        )
        return sample_id, deconv, deconv_frac

    def load_image_sample(self, sample_id: str) -> tuple[str, np.ndarray]:
        """Load the counts data for a given sample."""
        uns_spatial, obsm_spatial = self._dict_visium_samples[sample_id].load_images()
        return sample_id, uns_spatial, obsm_spatial

    def load_images_cohort(self) -> dict[str, Any]:

        images_list = Parallel(n_jobs=-1, prefer="threads")(
            delayed(self.load_image_sample)(sample_id)
            for sample_id in self.list_sample_ids
        )

        uns_spatial = {}
        obsm_spatial = {}
        for (sample_id, uns_, obsm_) in images_list:
            uns_spatial[sample_id] = uns_
            obsm_spatial[sample_id] = obsm_

        return {'uns_spatial':uns_spatial, 'obsm_spatial':obsm_spatial}

    def load_deconvolution_cohort(
        self,
        celltype_mapping_dict=None,
        frac: bool = False,
        keep_original_celltypes: bool = False,
        reload: bool = False,
    ) -> dict[str, np.ndarray]:
        """Load the deconvolution data for all samples in the cohort."""
        if celltype_mapping_dict is not None:
            if not keep_original_celltypes:
                self._common_celltypes = [
                    celltype
                    for celltype in self._common_celltypes
                    if celltype
                    not in [
                        item
                        for sublist in celltype_mapping_dict.values()
                        for item in sublist
                    ]
                ]
            self._common_celltypes += list(celltype_mapping_dict.keys())
            self._common_celltypes = list(np.unique(self._common_celltypes))
        deconvolution_dict: dict[str, Any] = {}

        # Parallel version
        deconv_list = Parallel(n_jobs=-1, prefer="threads")(
            delayed(self.load_deconv_sample)(
                sample_id,
                celltype_mapping_dict=celltype_mapping_dict,
                keep_original_celltypes=keep_original_celltypes,
                reload=reload,
            )
            for sample_id in self.list_sample_ids
        )

        # Necessary step in parallel version
        # Otherwise the sample's attribute are not modified
        # including the list_celltypes, deconvolution, deconvolution_frac, loaded_deconvolution attributes
        for sample_id_deconvs in deconv_list:
            id, deconvolution, deconvolution_frac = sample_id_deconvs
            sample = self._dict_visium_samples[id]
            sample._list_celltypes = self._common_celltypes
            if not sample._loaded_deconvolution or reload:
                if isinstance(deconvolution, pd.DataFrame):
                    sample._deconvolution, sample._deconvolution_frac = (
                        np.array(deconvolution.loc[:, self._common_celltypes].values),
                        np.array(
                            deconvolution_frac.loc[:, self._common_celltypes].values
                        ),
                    )
                else:
                    # assume deconvolution is a numpy array
                    index_common_celltypes = [
                        i
                        for i, celltype in enumerate(sample._list_celltypes)
                        if celltype in self._common_celltypes
                    ]
                    sample._deconvolution, sample._deconvolution_frac = (
                        deconvolution[:, index_common_celltypes],
                        deconvolution_frac[:, index_common_celltypes],
                    )
                sample._loaded_deconvolution = True

            if frac:
                deconvolution_dict[id] = deconvolution_frac
            else:
                deconvolution_dict[id] = deconvolution

        self._loaded_deconvolution = True

        return deconvolution_dict

    def update_neighbors_sample(self, sample_id: str) -> tuple[str, pd.DataFrame]:
        """Load the neighbors' indices for a given sample."""
        assert self._radius is not None, (
            "You must set a radius before updating neighbors"
        )
        neighbors = self._dict_visium_samples[sample_id].update_neighbors(
            radius=self._radius,
            updated_from_cohort=True,
        )
        return sample_id, neighbors

    def update_neighbors_cohort(self, radius: int):
        """Update the neighbors with a new radius for all samples in the cohort."""
        self._radius = radius
        neighbor_dict = {}
        # Parallel version
        neighbors_list = Parallel(n_jobs=-1, prefer="threads")(
            delayed(self.update_neighbors_sample)(sample_id)
            for sample_id in self.list_sample_ids
        )
        # Necessary step in parallel version
        # Otherwise the sample's attribute are not modified
        # including the radius and neighbors_idx attributes
        for sample_id_neighbors in neighbors_list:
            id = sample_id_neighbors[0]
            neighbors = sample_id_neighbors[1]
            self._dict_visium_samples[id]._radius = self._radius
            self._dict_visium_samples[id].neighbors_idx = neighbors
        neighbor_dict = dict(neighbors_list)
        self.neighbors_idx = neighbor_dict
        logger.info(
            f"Neighbors have been updated for the cohort with radius {radius}, adata "
            "and sdata objects should be reloaded"
        )
        return neighbor_dict

    def filter_umi(
        self,
        umi_threshold: int | None = None,
        warn_threshold: float = 0.5,
        remove_threshold: float = 0.25,
        all_genes: bool = True,
        reload: bool = False,
    ):
        """Filter the spots with valid UMI counts."""
        if not self._loaded_counts or reload:
            raise ValueError("You must load the counts before filtering")
        if umi_threshold is None:
            warnings.warn(
                "No umi threshold found in the zarr store, will use a per-sample "
                "threshold if one is found in the zarr stores.",
                stacklevel=2,
            )
        updated_coords = pd.DataFrame()
        updated_neighbors_idx: dict[str, Any] = {}
        for sample_id in self.list_sample_ids:
            try:
                self._dict_visium_samples[sample_id].filter_umi(
                    umi_threshold=umi_threshold,
                    warn_threshold=warn_threshold,
                    remove_threshold=remove_threshold,
                    all_genes=all_genes,
                )
                updated_coords = pd.concat(
                    [updated_coords, self._dict_visium_samples[sample_id]._coords]
                )
                if bool(self.neighbors_idx):
                    updated_neighbors_idx[sample_id] = self._dict_visium_samples[
                        sample_id
                    ].neighbors_idx
            except Exception as e:
                logger.error(f"{e} The sample is removed from the cohort.")
                self._list_sample_ids.remove(sample_id)
                self._dict_visium_samples.pop(sample_id)

        self._coords = updated_coords
        self.neighbors_idx = updated_neighbors_idx

    def get_sdata(self, sample_id: str) -> sd.SpatialData:
        """Return the SpatialData object for a given sample."""
        return self._dict_visium_samples[sample_id].get_sdata()

    def get_adata(
        self, sample_id: str, load_images: bool = False, load_obs: bool = False
    ) -> ad.AnnData:
        """Return the AnnData object for a given sample."""
        return self._dict_visium_samples[sample_id].get_adata(
            load_images=load_images, load_obs=load_obs
        )

    def get_adata_concat(
        self, load_images: bool = False, load_obs: bool = False
    ) -> ad.AnnData:
        """Return the concatenated AnnData object for all samples."""
        adata_list = []
        for sample_id in self.list_sample_ids:
            adata_sample = self._dict_visium_samples[sample_id].get_adata(
                load_images=load_images, load_obs=load_obs
            )
            adata_sample.obs["sample_id"] = sample_id
            adata_sample.obs_names = adata_sample.obs["sample_id"].str.cat(
                adata_sample.obs_names, sep="."
            )
            adata_list.append(adata_sample)
        adata = ad.concat(adata_list)
        adata.obs.index.name = None
        return adata

    def getitem_sample(self, sample_id: str) -> dict[str, Any]:
        """Return the sample object."""
        return self._dict_visium_samples[sample_id].getitem_sample()

    def get_pseudobulk(self) -> pd.DataFrame:
        """Return pseudobulk counts."""
        counts_cohort = self.load_counts_cohort()
        pseudobulk = pd.DataFrame(
            np.stack([counts.mean(0) for counts in counts_cohort.values()]),
            columns=self.common_genes,
            index=self.list_sample_ids,
        )
        return pseudobulk

    def __getitem__(self, index: str | list[str]) -> VisiumCohort | VisiumSample:
        """Return a sample or a sub-cohort.

        Args
        ----
        index : str | list[str]
            Sample id or list of sample ids to return.

        Returns
        -------
        VisiumCohort | VisiumSample
            The sample or the sub-cohort.
        """
        if isinstance(index, str):
            # assert that index an element of self._list_sample_ids
            assert index in self._list_sample_ids, (
                "Required sample id is not in the original cohort"
            )
            return self._dict_visium_samples[index]  # type: ignore

        if isinstance(index, list):
            # assert that index is included in self._list_sample_ids
            assert all(sample_id in self._list_sample_ids for sample_id in index), (
                "Some required sample ids are not in the original cohort"
            )

            sub_cohort = VisiumCohort(
                path_dataset=self._path_dataset,
                normalization=self._normalization,
                radius=self._radius,
                deconvolution_slot_name=self._deconvolution_slot_name,
                list_sample_ids=index,
                ignore_sample_ids=None,
                gene_threshold=self._gene_threshold,
            )

            # the following attributes should be the same as the original cohort
            sub_cohort._common_genes = self._common_genes
            sub_cohort._common_celltypes = self._common_celltypes
            sub_cohort._loaded_counts = self._loaded_counts
            sub_cohort._loaded_deconvolution = self._loaded_deconvolution

            # in case there has been some umi filtering at the cohort level
            for id in sub_cohort._list_sample_ids:
                sub_cohort._dict_visium_samples[
                    id
                ]._list_barcodes = self._dict_visium_samples[id]._list_barcodes
                sub_cohort._dict_visium_samples[id]._coords = self._dict_visium_samples[
                    id
                ]._coords
                if self._dict_visium_samples[id].neighbors_idx is not None:
                    sub_cohort._dict_visium_samples[
                        id
                    ].neighbors_idx = self._dict_visium_samples[id].neighbors_idx

            # avoid reloading the counts and deconvolution if they already were
            if self._loaded_counts:
                for id in sub_cohort._list_sample_ids:
                    sub_cohort._dict_visium_samples[
                        id
                    ]._counts = self._dict_visium_samples[id]._counts
                    sub_cohort._dict_visium_samples[
                        id
                    ]._loaded_counts = self._dict_visium_samples[id]._loaded_counts
                    sub_cohort._dict_visium_samples[
                        id
                    ]._list_gene_names = self._dict_visium_samples[id]._list_gene_names

            if self._loaded_deconvolution:
                for id in sub_cohort._list_sample_ids:
                    sub_cohort._dict_visium_samples[
                        id
                    ]._deconvolution_frac = self._dict_visium_samples[
                        id
                    ]._deconvolution_frac
                    sub_cohort._dict_visium_samples[
                        id
                    ]._deconvolution = self._dict_visium_samples[id]._deconvolution
                    sub_cohort._dict_visium_samples[
                        id
                    ]._loaded_deconvolution = self._dict_visium_samples[
                        id
                    ]._loaded_deconvolution
                    sub_cohort._dict_visium_samples[
                        id
                    ]._list_celltypes = self._dict_visium_samples[id]._list_celltypes

            return sub_cohort

    def __len__(self):
        return len(self._list_sample_ids)

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index < len(self):
            sample_id = self._list_sample_ids[self.index]
            sample = self[sample_id]
            self.index += 1
            return sample_id, sample
        else:
            raise StopIteration

    @classmethod
    def from_visium_samples(
        cls,
        samples: list[VisiumSample],
        gene_combination_method: Literal["intersection", "union"] = "intersection",
        load_counts: bool = False,
        load_deconvolution: bool = False,
        copy_samples: bool = True,
        **kwargs,
    ) -> VisiumCohort:
        """Create a VisiumCohort from a list of VisiumSample objects.

        All samples must have the same normalization and deconvolution_slot_name.

        Parameters
        ----------
        samples : list[VisiumSample]
            List of VisiumSample objects to include in the cohort.
        gene_combination_method : Literal["intersection", "union"]
            How to match the genes between samples.
        load_counts : bool
            Whether to load counts for all samples after creating the cohort.
            If some samples already have counts loaded, they will be reloaded
            to ensure consistency with the common genes.
        load_deconvolution : bool
            Whether to load deconvolution for all samples after creating the cohort.
            If some samples already have deconvolution loaded, they will be reloaded
            to ensure consistency with the common cell types.
        copy_samples : bool
            Whether to copy the VisiumSample objects into the VisiumCohort object.
            If False, the VisiumSample objects will be references to the original
            objects.
        **kwargs
            Additional arguments to pass to the VisiumCohort constructor.

        Returns
        -------
        VisiumCohort
            A new VisiumCohort object containing the provided samples.

        Raises
        ------
        ValueError
            If the samples have different normalization or deconvolution_slot_name.
        ValueError
            If the samples have duplicate sample IDs.
        ValueError
            If the samples have different radius.
        """
        assert len(samples) > 0, "At least one sample must be provided"
        assert all(isinstance(sample, VisiumSample) for sample in samples), (
            "All samples must be VisiumSample objects"
        )

        # Verify all samples have the same normalization and deconvolution_slot_name
        normalization = samples[0].normalization
        deconvolution_slot_name = samples[0]._deconvolution_slot_name
        radius = samples[0]._radius

        for sample in samples[1:]:
            if sample.normalization != normalization:
                raise ValueError(
                    f"All samples must have the same normalization. Found {normalization} and {sample.normalization}"
                )
            if sample._deconvolution_slot_name != deconvolution_slot_name:
                raise ValueError(
                    f"All samples must have the same deconvolution_slot_name. Found {deconvolution_slot_name} and {sample._deconvolution_slot_name}"
                )
            if sample._radius != radius:
                warnings.warn(
                    f"Sample {sample.sample_id} has radius {sample._radius} which differs from {radius}. Using {radius} for the cohort.",
                    UserWarning,
                    stacklevel=2,
                )

        # Check for duplicate sample IDs
        sample_ids = [sample.sample_id for sample in samples]
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError(
                "Duplicate sample IDs found. Each sample must have a unique ID."
            )

        path_dataset = [sample._path_sample_zarr for sample in samples]

        # Create a minimal cohort object
        cohort = cls(
            path_dataset=path_dataset,
            normalization=normalization,
            radius=radius,
            deconvolution_slot_name=deconvolution_slot_name,
            list_sample_ids=sample_ids,
            gene_combination_method=gene_combination_method,
            **kwargs,
        )

        # Replace the automatically created samples with the provided ones
        if copy_samples:
            cohort._dict_visium_samples = {
                sample.sample_id: copy.deepcopy(sample) for sample in samples
            }
        else:
            cohort._dict_visium_samples = {
                sample.sample_id: sample for sample in samples
            }

        # Recalculate cohort-level attributes
        list_gene_names = [sample.list_gene_names for sample in samples]
        list_celltypes = [sample.list_celltypes for sample in samples]

        # Update coordinates
        cohort._coords = pd.DataFrame()
        for sample in samples:
            cohort._coords = pd.concat([cohort._coords, sample._coords])

        # Update neighbors if available
        if radius is not None:
            cohort.neighbors_idx = {}
            for sample in samples:
                if sample.neighbors_idx is not None:
                    cohort.neighbors_idx[sample.sample_id] = sample.neighbors_idx

        # Update common genes and celltypes
        cohort._common_genes = get_common_list(
            list_gene_names, how=gene_combination_method
        )
        cohort._common_celltypes = get_common_list(list_celltypes)

        # Check if all samples have counts and deconvolution loaded
        counts_loaded = all(sample._loaded_counts for sample in samples)
        deconv_loaded = all(sample._loaded_deconvolution for sample in samples)

        # Check for inconsistent loading state
        some_counts_loaded = any(sample._loaded_counts for sample in samples)
        some_deconv_loaded = any(sample._loaded_deconvolution for sample in samples)

        if some_counts_loaded and not counts_loaded and not load_counts:
            warnings.warn(
                "Some samples have counts loaded but others don't. This may cause inconsistent behavior. "
                "Consider using load_counts=True to ensure all samples have counts loaded.",
                UserWarning,
                stacklevel=2,
            )

        if some_deconv_loaded and not deconv_loaded and not load_deconvolution:
            warnings.warn(
                "Some samples have deconvolution loaded but others don't. This may cause inconsistent behavior. "
                "Consider using load_deconvolution=True to ensure all samples have deconvolution loaded.",
                UserWarning,
                stacklevel=2,
            )

        # Set loaded flags based on sample state
        cohort._loaded_counts = counts_loaded
        cohort._loaded_deconvolution = deconv_loaded

        # Load counts if requested
        if load_counts:
            logger.info("Loading counts for all samples to ensure consistency")
            cohort.load_counts_cohort(reload=True)

        # Load deconvolution if requested
        if load_deconvolution:
            logger.info("Loading deconvolution for all samples to ensure consistency")
            cohort.load_deconvolution_cohort(reload=True)

        return cohort
