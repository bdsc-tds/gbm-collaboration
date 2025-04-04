"""Visium sample class."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as sd
import zarr
from anndata._io.specs import read_elem
from loguru import logger
from spatialdata_io.experimental import from_legacy_anndata

from gbmcollaboration.utils import (
    check_path,
    compute_neighbors_as_multiindex,
    extract_counts_dataframe,
    extract_deconvolution_dataframe,
    load_coords_as_multiindex,
    read_elem_with_path,
)


class VisiumSample:
    """Dataset class loading a single sample data from a Zarr store.

    Args
    ----
    path_sample_zarr : str
        path to the zarr store or a glob pattern '/path/to/dataset/*.zarr'.
    normalization : str
        The normalization to use for the dataset.
    sample_id : str | None
        The sample id to load.
    radius : int | None
        If not None, compute the neighbors of the spots with the given radius.
    deconvolution_slot_name : str
        The name of the deconvolution wanted when several are available.
        By default, "deconvolution" is used, as it was the original name in
        Mosaic visium pipeline.
    indices : np.ndarray | None
        If not None, filter spots with the given indices.

    Raises
    ------
    ValueError
        If the deconvolution slot name is not available.
    """

    def __init__(
        self,
        path_sample_zarr: str,
        normalization: str,
        sample_id: str | None = None,
        radius: int | None = None,
        deconvolution_slot_name: str = "deconv_cohort_level_2",
        indices: np.ndarray | None = None,
    ):
        # Check if provided path match any files
        check_path(path_sample_zarr)
        self._path_sample_zarr = path_sample_zarr
        self._zarr_store = zarr.open(path_sample_zarr, "r")
        if sample_id is None:
            sample_id = Path(path_sample_zarr).stem
        self._sample_id = sample_id
        self._normalization = normalization
        self._radius = radius
        self._deconvolution_slot_name = deconvolution_slot_name
        self._counts: np.ndarray | None = None
        self._deconvolution: np.ndarray | None = None
        self._deconvolution_frac: np.ndarray | None = None

        # fill sample level attributes
        self._list_gene_names = read_elem_with_path(
            self._zarr_store, ["tables", "table", "var", "_index"]
        )
        available_normalizations = list(self._zarr_store["tables/table/layers"].keys())
        if normalization not in available_normalizations:
            raise ValueError(
                f"Normalization '{normalization}' is not available. "
                f"Available normalizations are: {available_normalizations}"
            ) from None

        if "umi_threshold" in self._zarr_store["tables/table/uns"]:
            self._umi_threshold = read_elem_with_path(
                self._zarr_store,
                ["tables", "table", "uns", "umi_threshold", "threshold"],
            ).item()
        else:
            self._umi_threshold = None

        # get the deconvolution table, the name is given by the user
        # or is by default "deconvolution"
        available_deconvolution = list(self._zarr_store["tables/table/obsm"].keys())
        # get names ending with weights
        deconvolution_available = [
            name for name in available_deconvolution if name.endswith("_weights")
        ]
        # remove string _frac and _weight for each available name
        if len(deconvolution_available) > 0:
            deconvolution_available = [
                name.replace("_weights", "") for name in deconvolution_available
            ]
            deconvolution_available = list(set(deconvolution_available))
            self._deconvolution_available = deconvolution_available
        try:
            if deconvolution_slot_name == "deconvolution":
                self._list_celltypes = read_elem_with_path(
                    self._zarr_store,
                    ["tables", "table", "obsm", f"{self._deconvolution_slot_name}"],
                ).columns.tolist()
            else:
                self._list_celltypes = read_elem_with_path(
                    self._zarr_store,
                    [
                        "tables",
                        "table",
                        "obsm",
                        f"{self._deconvolution_slot_name}_weights",
                    ],
                ).columns.tolist()
            self._original_list_celltypes = self._list_celltypes
        except KeyError:
            raise ValueError(
                f"Deconvolution slot name '{deconvolution_slot_name}' is not available for sample {sample_id}. "
                f"Available deconvolution slots are: {deconvolution_available}"
                f"If you are loading an actually available slot, please inform #mosaic-data-users, it may be a corrupted zarr store."
            ) from None

        self._coords = load_coords_as_multiindex(self._zarr_store, self._sample_id)
        self._list_barcodes = list(self._coords["barcode"])
        if indices is not None:
            self._list_barcodes = list(np.array(self._list_barcodes)[indices])
            self._coords = self._coords[
                self._coords["barcode"].isin(self._list_barcodes)
            ]

        self.neighbors_idx: pd.DataFrame | None = None
        if self._radius is not None:
            self.neighbors_idx = compute_neighbors_as_multiindex(
                self._coords[["norm_x", "norm_y"]].loc[self._sample_id],  # type: ignore[arg-type]
                radius=self._radius,
                sample_id=self._sample_id,
            )

        self._loaded_counts = False
        self._loaded_deconvolution = False

    @property
    def sample_id(self) -> str:
        """Return the sample id."""
        return self._sample_id

    @property
    def normalization(self) -> str:
        """Return the normalization."""
        return self._normalization

    @property
    def list_gene_names(self) -> list[str]:
        """Return the genes."""
        return self._list_gene_names

    @property
    def list_celltypes(self) -> list[str]:
        """Return the cell types."""
        return self._list_celltypes

    @property
    def original_list_celltypes(self) -> list[str]:
        """Return the original cell types."""
        return self._original_list_celltypes

    @property
    def list_barcodes(self) -> list[str]:
        """Return the barcodes."""
        return self._list_barcodes

    @property
    def coords(self) -> pd.DataFrame:
        """Return the coordinates of the spots."""
        return self._coords

    @property
    def counts(self) -> np.ndarray | None:
        """Return the counts."""
        return self._counts

    @property
    def deconvolution(self) -> np.ndarray | None:
        """Return the deconvolution."""
        return self._deconvolution

    @property
    def deconvolution_frac(self) -> np.ndarray | None:
        """Return the deconvolution fraction."""
        return self._deconvolution_frac

    @property
    def deconvolution_available(self) -> list[str] | None:
        """Return the deconvolution slots available."""
        return self._deconvolution_available

    def __len__(self) -> int:
        return len(self._coords)

    def get_adata(
        self, load_images: bool = False, load_obs: bool = False
    ) -> ad.AnnData:
        """Format the VisiumSample object to an AnnData object."""
        if not self._loaded_counts:
            raise ValueError("You must load the counts before running the analysis")

        if load_images or load_obs:
            # Create a mask of vlaid barcodes. User to filter barcode only if loading image or obs
            all_barcodes = list(self._zarr_store["tables"]["table"]["obs"]["_index"])
            valid_barcodes = self.list_barcodes
            mask_valid_barcodes = [
                barcode in valid_barcodes for barcode in all_barcodes
            ]

        counts = pd.DataFrame(
            self.counts,
            index=self.coords["barcode"],
            columns=self.list_gene_names,
        )

        if load_obs:
            obs = read_elem(self._zarr_store["tables"]["table"]["obs"])
            obs.index.name = "barcode"
            obs = obs.loc[valid_barcodes]
        else:
            obs = self.coords.set_index("barcode")

        adata = ad.AnnData(
            X=counts,
            layers={
                self.normalization: counts,
            },
            obs=obs,
        )

        if self._loaded_deconvolution:
            deconvolution = pd.DataFrame(
                self.deconvolution,
                index=self.coords["barcode"],
                columns=self.list_celltypes,
            )
            deconvolution_frac = pd.DataFrame(
                self.deconvolution_frac,
                index=self.coords["barcode"],
                columns=self.list_celltypes,
            )
            adata.obsm["deconvolution"] = deconvolution
            adata.obsm["deconvolution_frac"] = deconvolution_frac
            adata.obs["DominantCellType"] = np.array(deconvolution.columns)[
                np.argmax(deconvolution_frac, axis=1)
            ]

        if self.neighbors_idx is not None:
            neighbors = pd.DataFrame(
                self.neighbors_idx.values,
                index=self.coords["barcode"],
            )
            adata.obsm["neighbors_index"] = neighbors

        if load_images:
            adata.uns["spatial"] = read_elem(
                self._zarr_store["tables"]["table"]["uns"]["spatial"]
            )
            adata.obsm["spatial"] = read_elem(
                self._zarr_store["tables"]["table"]["obsm"]["spatial"]
            )[mask_valid_barcodes, :]

        return adata

    def get_sdata(self) -> sd.SpatialData:
        """Format the VisiumSample object to a SpatialData object."""
        return from_legacy_anndata(self.get_adata())

    def load_counts(
        self,
        list_common_genes: list | None = None,
        store: bool = True,
        reload: bool = False,
    ):
        """Load the visium counts data for a given sample."""
        if list_common_genes is None:
            list_common_genes = self._list_gene_names
        self._list_gene_names = list_common_genes
        if not self._loaded_counts or reload:
            counts_df = extract_counts_dataframe(
                self._zarr_store,
                self._normalization,
                list_valid_spot=self._list_barcodes,
            )
            missing_genes = [
                gene for gene in list_common_genes if gene not in counts_df.columns
            ]  # useful only when using union of genes. Impute undetected genes with 0
            counts_df.loc[:, missing_genes] = 0
            counts_ = counts_df.loc[:, list_common_genes].values

            if store:
                self._counts = np.array(counts_)
                self._loaded_counts = True
            return counts_
        return self._counts

    def load_deconvolution(
        self,
        list_common_cell_types: list | None = None,
        celltype_mapping_dict: dict | None = None,
        store: bool = True,
        keep_original_celltypes: bool = False,
        reload: bool = False,
    ):
        """Load the deconvolution data for a given sample."""
        if celltype_mapping_dict is not None:
            if not keep_original_celltypes:
                self._list_celltypes = [
                    celltype
                    for celltype in self._list_celltypes
                    if celltype
                    not in [
                        item
                        for sublist in celltype_mapping_dict.values()
                        for item in sublist
                    ]
                ]
            self._list_celltypes += list(celltype_mapping_dict.keys())
            self._list_celltypes = list(np.unique(self._list_celltypes))
            list_common_cell_types = self._list_celltypes
        if list_common_cell_types is None:
            list_common_cell_types = self._list_celltypes
        self._list_celltypes = list_common_cell_types

        if not self._loaded_deconvolution or reload:
            if self._deconvolution_slot_name == "deconvolution":
                to_extract = f"{self._deconvolution_slot_name}"
            else:
                to_extract = f"{self._deconvolution_slot_name}_weights"
            df_deconvolution, df_deconvolution_frac = extract_deconvolution_dataframe(
                self._zarr_store,
                list_valid_spot=self._list_barcodes,
                refactor_dict=celltype_mapping_dict,
                deconvolution_slot_name=to_extract,
                keep_original_celltypes=keep_original_celltypes,
            )

            if store:
                self._deconvolution, self._deconvolution_frac = (
                    np.array(df_deconvolution.loc[:, list_common_cell_types].values),
                    np.array(
                        df_deconvolution_frac.loc[:, list_common_cell_types].values
                    ),
                )
                self._loaded_deconvolution = True
            return (
                df_deconvolution.loc[:, list_common_cell_types].values,
                df_deconvolution_frac.loc[:, list_common_cell_types].values,
            )

        return self._deconvolution, self._deconvolution_frac

    def remove_spots(self, barcodes: str | list[str]):
        """Remove spots based on a list of spot barcodes."""
        if isinstance(barcodes, str):
            barcodes = [barcodes]

        # Remove duplicated barcodes
        barcodes_unique = np.unique(barcodes)

        # Find spot indices of the spots to keep
        to_keep_idx = np.argwhere(
            ~np.isin(self._list_barcodes, barcodes_unique)
        ).flatten()
        if len(to_keep_idx) + len(barcodes_unique) > len(self._list_barcodes):
            wrong_barcodes = list(set(barcodes_unique) - set(self._list_barcodes))
            logger.warning(f"Some barcodes to remove were not found: {wrong_barcodes}")

        # Remove the spots
        self._list_barcodes = [b for b in self._list_barcodes if b not in barcodes]
        self._coords = self._coords.iloc[to_keep_idx]
        if self._loaded_counts:
            self._counts = self._counts[to_keep_idx]  # type: ignore[index]
        if self._loaded_deconvolution:
            self._deconvolution = self._deconvolution[to_keep_idx]  # type: ignore[index]
            self._deconvolution_frac = self._deconvolution_frac[to_keep_idx]  # type: ignore[index]
        if self.neighbors_idx is not None:
            self.neighbors_idx = compute_neighbors_as_multiindex(
                self._coords[["norm_x", "norm_y"]].loc[self._sample_id],  # type: ignore[arg-type]
                radius=self._radius,  # type: ignore[arg-type]
                sample_id=self._sample_id,
            )

    def filter_umi(
        self,
        umi_threshold: int | None = None,
        warn_threshold: float = 0.5,
        remove_threshold: float = 0.25,
        all_genes: bool = True,
    ):
        """Filter the spots with valid UMI counts."""
        if not self._loaded_counts:
            raise ValueError("You must load the counts before filtering")
        if umi_threshold is None:
            if self._umi_threshold is None:
                raise ValueError(
                    "No umi threshold found in the zarr store, please provide a "
                    "threshold"
                )
            umi_threshold = self._umi_threshold
        if self._normalization != "raw_count" and not all_genes:
            total_counts = (
                extract_counts_dataframe(
                    self._zarr_store,
                    "raw_count",
                    list_valid_spot=self._list_barcodes,
                )
                .loc[:, self._list_gene_names]
                .values.sum(axis=1)
            )
        elif all_genes:
            total_counts = extract_counts_dataframe(
                self._zarr_store,
                "raw_count",
                list_valid_spot=self._list_barcodes,
            ).values.sum(axis=1)
        else:
            total_counts = self._counts.sum(axis=1)  # type: ignore[union-attr]
        index_umi_not_valid = total_counts < umi_threshold
        barcodes_to_remove = list(np.array(self._list_barcodes)[index_umi_not_valid])

        # Warn the user in case of low QC sample
        if 1 - len(barcodes_to_remove) / len(self._list_barcodes) < warn_threshold:
            print(
                f"Sample {self._sample_id} has {round((len(barcodes_to_remove) / len(self._list_barcodes)) * 100, 1)}% of spots with less than {umi_threshold} UMI"
            )

        # If we are about to filter too many spots, return to avoid errors
        if 1 - len(barcodes_to_remove) / len(self._list_barcodes) < remove_threshold:
            raise ValueError(
                f"Too many spots have been filtered out by filter_umi threshold for sample {self._sample_id}."
            )

        # Remove the spots with UMI below threshold
        if len(barcodes_to_remove) > 0:
            self.remove_spots(barcodes=barcodes_to_remove)

    def filter_artifact(
        self,
        warn_threshold: float = 0.5,
        umi_threshold_range: float = 0.15,
        type: str = "cassette_misalignment",
    ):
        """Filter the spots affected by artifacts (WIP: only cassette misalignment for now).

        Parameters
        ----------
        warn_threshold: float
            The threshold to warn the user about the quality of the sample.
        umi_threshold_range: float
            The range of the umi_threshold to consider for the filtering.
        type: str
            The type of artifact to consider for the filtering.
            Available options are: "cassette_misalignment", "other_artifacts" or "all".
            Only "cassette_misalignment" is available for now but soon other
            artifacts will be available from the Curator annotations.

        Raises
        ------
        ValueError
            If the counts are not loaded before filtering.
        """
        if not self._loaded_counts:
            raise ValueError("You must load the counts before filtering")

        barcodes_to_remove = []
        umi_range = self._umi_threshold + self._umi_threshold * umi_threshold_range

        if type in ["cassette_misalignment", "all"]:
            barcodes_to_remove.extend(self._filter_cassette_misalignment(umi_range))

        if type in ["other_artifacts", "all"]:
            barcodes_to_remove.extend(
                self._filter_other_artifacts()
            )  # TODO implement this method

        barcodes_to_remove = list(set(barcodes_to_remove))  # ensure uniqueness

        logger.info(f"{len(barcodes_to_remove)} barcodes have been removed")

        # Warn the user in case of low QC sample
        if 1 - len(barcodes_to_remove) / len(self._list_barcodes) < warn_threshold:
            print(
                f"Warning: sample {self._sample_id} has {round((len(barcodes_to_remove) / len(self._list_barcodes)) * 100, 1)}% of spots with less than {umi_range} UMI counts."
            )

        # If we are about to filter out all the spots, return to avoid errors
        if len(self._list_barcodes) == len(barcodes_to_remove):
            raise ValueError(
                "All spots have been filtered out by filter_artifact, check that this threshold is not stringent. If not, if you are using VisiumCohort class, remove this sample using `ignore_sample_ids` attribute. Otherwise if you are using VisiumSample class, don't consider this sample."
            )

        # Remove the spots with UMI below threshold
        if len(barcodes_to_remove) > 0:
            self.remove_spots(barcodes=barcodes_to_remove)

    def _filter_cassette_misalignment(self, umi_range: float):
        """Filter the spots affected by cassette misalignment."""
        # prepare the counts, barcodes and coordinates
        if self._normalization != "raw_count":
            total_counts = extract_counts_dataframe(
                self._zarr_store,
                "raw_count",
                list_valid_spot=self._list_barcodes,
            ).sum(axis=1)
        else:
            total_counts = self._counts.sum(axis=1)  # type: ignore[union-attr]

        counts_df = pd.DataFrame(
            self._counts, columns=self._list_gene_names, index=self._list_barcodes
        )
        total_counts = counts_df.sum(axis=1)
        spatial_coords = self._coords.reset_index()[
            ["barcode", "array_col", "array_row"]
        ]

        # calculate the row and column sums of counts
        row_sums = spatial_coords.groupby("array_row")["barcode"].apply(
            lambda barcodes: total_counts.loc[barcodes].mean()
        )
        col_sums = spatial_coords.groupby("array_col")["barcode"].apply(
            lambda barcodes: total_counts.loc[barcodes].mean()
        )

        # Get the umi_threshold that is sample-specific, use a range in case
        # the umi threshold is not accurate
        logger.info(f"UMI threshold for sample {self._sample_id} is {umi_range}")

        # for col_sums < umi_threshold_range, identify the spots within it
        low_rows = row_sums[row_sums < umi_range].index
        low_cols = col_sums[col_sums < umi_range].index

        misaligned_barcodes = list(
            spatial_coords[
                spatial_coords["array_col"].isin(low_cols)
                | spatial_coords["array_row"].isin(low_rows)
            ]["barcode"]
        )

        logger.info(f"{len(misaligned_barcodes)} spots predicted as misaligned")
        # need to ensure we remove known barcodes
        barcodes_to_remove = [
            barcode for barcode in self._list_barcodes if barcode in misaligned_barcodes
        ]

        # to improve the detect we can look at the neighbours of the misaligned spots
        # and remove those below the umi_threshold
        # first put a radius of 2 because 1 is not enough
        save_radius = self._radius
        self.update_neighbors(radius=2, updated_from_cohort=True)
        coords = self._coords.reset_index()[["barcode", "x", "y"]]
        coords_filtered = coords[coords["barcode"].isin(barcodes_to_remove)].set_index(
            ["x", "y"]
        )
        neigh_idx = (
            self.neighbors_idx.reset_index()  # type: ignore[union-attr]
            .set_index(["x", "y"])
            .drop(columns=["sample_id"])
        )
        neigh_idx = coords_filtered.merge(
            neigh_idx, left_index=True, right_index=True, how="left"
        )
        neigh_idx_values = neigh_idx.drop(columns=["barcode"]).values.flatten()
        neigh_idx_values = np.unique(neigh_idx_values)
        neigh_idx_values = neigh_idx_values[
            ~np.isnan(neigh_idx_values)
        ]  # drop the "nan"
        neigh_barcode = [self._list_barcodes[int(idx)] for idx in neigh_idx_values]
        total_counts = total_counts[neigh_barcode]
        barcodes_below_threshold = total_counts[total_counts < umi_range].index.tolist()
        # remove the barcodes that are already in the misaligned spots
        barcodes_below_threshold = [
            barcode
            for barcode in barcodes_below_threshold
            if barcode not in misaligned_barcodes
        ]
        barcodes_below_threshold = list(set(barcodes_below_threshold))
        logger.info(
            f"{len(barcodes_below_threshold)} neighbooring barcodes have been added to the misaligned spots"
        )
        # reset the radius
        self.update_neighbors(
            radius=save_radius,  # type: ignore[arg-type]
            updated_from_cohort=True,
        )

        barcodes_to_remove += barcodes_below_threshold

        return barcodes_to_remove

    def _filter_other_artifacts(self):
        """Filter other artifacts like hair, bubbles, etc."""
        logger.info("Filtering other artifacts is not implemented yet")
        return []

    def update_neighbors(self, radius: int, updated_from_cohort: bool = False):
        """Update the neighbors with a new radius."""
        self._radius = radius
        self.neighbors_idx = compute_neighbors_as_multiindex(
            self._coords[["norm_x", "norm_y"]].loc[self._sample_id],  # type: ignore[arg-type]
            radius=self._radius,
            sample_id=self._sample_id,
        )
        if not updated_from_cohort:
            logger.info(
                f"Neighbors have been updated with radius {radius} for sample "
                f"{self._sample_id}, adata and sdata objects should be reloaded"
            )
        return self.neighbors_idx

    def __getitem__(
        self,
        index: int,
        list_common_genes: list[str] | None = None,
        list_common_celltypes: list[str] | None = None,
    ) -> dict[str, Any]:
        output: dict[str, Any] = {
            "features": {},
            "coords": {},
        }

        output["features"]["counts"] = self.load_counts(list_common_genes)[index]
        output["features"]["deconvolution"] = self.load_deconvolution(
            list_common_celltypes
        )[0][index]
        output["features"]["deconvolution_frac"] = self.load_deconvolution(
            list_common_celltypes
        )[1][index]
        output["coords"]["self"] = self.coords.loc[self._sample_id].values[index]
        if self.neighbors_idx is not None:
            output["coords"]["neighbors_index"] = self.neighbors_idx.iloc[index].values
        return output

    def getitem_sample(self) -> dict[str, Any]:
        """Get the whole sample dict."""
        list_index = list(range(len(self)))
        out = self[list_index]  # type: ignore[index]
        return out

    def get_subset(self, list_index: list[int] | np.ndarray) -> VisiumSample:
        """Return a copy of the VisiumSample object with a subset of the spots."""
        subset = VisiumSample(
            path_sample_zarr=self._path_sample_zarr,
            sample_id=self._sample_id,
            normalization=self._normalization,
            radius=self._radius,
            indices=np.array(list_index),
        )
        return subset

    def extract_weighted_point_cloud(
        self,
        gene: str,
        xy_source: str = "norm",
        binarise: bool = False,
        binarise_threshold: float | None = 0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Extract weights and locations for given gene.

        Extract x,y coordinates / spatial positions and expression
        values as arrays from a VisiumSample - the dataset class
        from sowce for handling Visium data.

        Parameters
        ----------
        gene: str
            name of a gene whose expression values are required
        xy_source: str
            this arguments allows the user to determine whether
            the norm_x, norm_y values are used as the x,y coordinates
            (default, float) or the array_col, array_row
            values (any value other than norm, integer).
        binarise: bool
            Whether or not to convert the expression levels to binary.
            Defaults to False.
            Some scores require this.
        binarise_threshold: float | None
            Threshold of expression level to use for binarising.
            Defaults to zero.

        Returns
        -------
        a_i: np.ndarray
            weights
        x_i: np.ndarray
            locations
        """
        # get expression of gene in a 1D array
        assert gene in self._list_gene_names, (
            f"Your gene {gene} is not present in the sample"
        )
        # get counts for the gene depending on whether load_counts
        # already executed or not
        if self._loaded_counts:
            counts = (
                pd.DataFrame(
                    self._counts,
                    index=self._list_barcodes,
                    columns=self._list_gene_names,
                )[[gene]]
                .to_numpy()
                .flatten()
            )
        else:
            counts = self.load_counts(list_common_genes=[gene], store=False).flatten()

        # get the x, y coords array positions
        xy_source_cols = (
            ["norm_x", "norm_y"] if xy_source == "norm" else ["array_col", "array_row"]
        )
        vals = pd.concat(
            [
                self._coords[["barcode"] + xy_source_cols].set_index("barcode"),
                pd.Series(counts, name=gene, index=self._list_barcodes),
            ],
            axis=1,
        )
        # whether to binarise
        x_i = vals[xy_source_cols].to_numpy()
        a_i = vals[gene].to_numpy()
        if binarise:
            x_i = x_i[a_i > binarise_threshold]  # type: ignore
            a_i = np.ones(len(x_i))  # type: ignore
        return a_i, x_i
