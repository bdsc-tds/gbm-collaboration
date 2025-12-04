"""Define the DataCenter class."""

from __future__ import annotations

import re
import functools
import gc
import scanpy as sc
import numpy as np
import pandas as pd
from typing import NamedTuple, Optional
from joblib import Parallel, delayed

from gbmcollaboration.data.io.data_structure import DataSource
from gbmcollaboration.data.io.loaders import SingleCellLoader, VisiumCohortLoader
from gbmcollaboration.definitions import (
    BAD_QUALITY_CHUV_VISIUM,
    BAD_QUALITY_UKER_VISIUM,
    BAD_QUALITY_CR_VISIUM,
    PATIENTS_IDH1_MUT
)
from gbmcollaboration.utils.spatial_obj_to_df import (
    celltype_array_to_df,
    genecounts_array_to_df,
)


CENTERS: dict[str, DataCenter] = {}


def register_center(cls):
    """Define decorator to register instances of DataCenter."""

    @functools.wraps(cls)
    def wrapper_register(*args, **kwargs):
        obj = cls(*args, **kwargs)
        CENTERS[obj.name] = obj
        return obj

    return wrapper_register


@register_center
class DataCenter(NamedTuple):
    """Define the structure of a data center."""

    name: str
    sources: dict[str, DataSource]

    def load_tabular(self, exclude_idh1_mut: bool = False):
        """Load tabular data."""

        def exclude_idh1mut(df: pd.DataFrame) -> pd.DataFrame:
            """
            Parameters
            ----------
            exclude_idh1_mut : bool
                If True, removes samples from the data whose indices are in the predefined
                list of IDH1-mutated samples, by default False.
            Returns
            -------
            dict[str, pd.DataFrame]
                Dictionary of dataframes for each source that are tabular.
            """
            # remove the suffix from visium
            mutated_patient = PATIENTS_IDH1_MUT.copy()
            mutated_patient = [s[:-5] for s in mutated_patient]
            # convert mutated_patient to a set for faster lookups
            patients_set = set(mutated_patient)
            # create a pattern to match any of the patient IDs in the set
            pattern = "|".join(map(re.escape, patients_set))

            # check if any column or index contains the pattern
            column_matches = df.columns[df.columns.astype(str).str.contains(pattern)]
            index_matches = df.index[df.index.astype(str).str.contains(pattern)]

            if len(column_matches) > 0:
                print(f"IDH1mut: excluding {len(column_matches)} columns from tabular data.")
                df = df.drop(columns=column_matches)
            elif len(index_matches) > 0:
                print(f"IDH1mut: excluding {len(index_matches)} rows from tabular data.")
                df = df.drop(index=index_matches)
            else:
                print("IDH1 mutation samples not found in the columns or index.")
            return df

        def load_and_filter(file_obj):
            df = file_obj.load()
            if exclude_idh1_mut:
                df = exclude_idh1mut(df)
            return df

        return {
            source_name: {
                file_name: load_and_filter(file_)
                for file_name, file_ in source.files.items()
            } if isinstance(source.files, dict) else load_and_filter(source.files)
            for source_name, source in self.sources.items()
            if source_name in ["clinical", "bulk_rna", "wes"]#, "he"]
        }

    def load_spatial(
        self,
        exclude_low_quality_visium: bool = False,
        exclude_idh1_mut: bool = False,
        list_sample_id: Optional[list[str]] = None,
        normalization: str = "raw_count",
        radius: int = 1,
        deconvolution_slot_name: str = "deconv_cohort_level_4",
        load_gex: bool = True,
        load_deconvolution: bool = True,
        load_images: bool = False,
    ):
        """Load spatial data of the full cohort or a list of samples.

        Parameters
        ----------
        exclude_low_quality_visium : bool
            Exclude samples with low quality.
        exclude_idh1_mut : bool
            If True, removes samples from the data whose indices are in the predefined
            list of IDH1-mutated samples, by default False.
        list_sample_id : list[str]
            List of sample IDs to load. If None, load the full cohort.
            Default is None.
        normalization : str
            Normalization method.
            Default is "raw_count". Supported methods are "raw_count",
            "log_norm_count", "norm_count", "SCTransform".
        radius : int
            Radius for the spatial data.
            Default is 1.
        deconvolution_slot_name : str
            Slot name for the deconvolution data.
            Default is "deconv_GMBatlas_2steps".
        load_gex: bool
            Whether to load gene expression
        load_deconvolution: bool
            Whether to load deconvolutions
        load_images: bool
            Whether to load images

        Returns
        -------
        dict[str, pd.DataFrame]
            Dictionary of dataframes for each source that are spatial.
        """
        source = self.sources.get("spatial")
        # initialize ignore_sample_ids
        ignore_sample_ids = None

        if self.name == "mosaic_uker" and exclude_low_quality_visium:
            len_bad_quality = len(BAD_QUALITY_UKER_VISIUM)
            print(f"Excluding {len_bad_quality} bad quality samples at MOSAIC-UKER")
            ignore_sample_ids = BAD_QUALITY_UKER_VISIUM.copy()
        elif self.name == "mosaic_chuv" and exclude_low_quality_visium:
            len_bad_quality = len(BAD_QUALITY_CHUV_VISIUM)
            print(f"Excluding {len_bad_quality} bad quality samples at MOSAIC-CHUV")
            ignore_sample_ids = BAD_QUALITY_CHUV_VISIUM.copy()
        elif self.name == "mosaic_cr" and exclude_low_quality_visium:
            len_bad_quality = len(BAD_QUALITY_CR_VISIUM)
            print(f"Excluding {len_bad_quality} bad quality samples at MOSAIC-CR")
            ignore_sample_ids = BAD_QUALITY_CR_VISIUM.copy()

        if exclude_idh1_mut:
            if ignore_sample_ids is None:
                ignore_sample_ids = PATIENTS_IDH1_MUT
            else:
                ignore_sample_ids += PATIENTS_IDH1_MUT
            print("Excluding IDH1 mutated samples if any present.")

        if source:
            # instantiate VisiumDataLoader
            print( "Instantiating loader.")
            path_dataset = source.files.get("cohort zarr files").name  # type: ignore[union-attr]
            data_loader = VisiumCohortLoader(
                path_dataset=path_dataset,
                ignore_sample_ids=ignore_sample_ids,
                list_sample_ids=list_sample_id,
                normalization=normalization,
                radius=radius,
                deconvolution_slot_name=deconvolution_slot_name,
            )
            # preparing the dataframes as additional attributes
            spots_info = data_loader.coords
            spots_info.reset_index(inplace=True)
            if load_gex:
                print(
                    "Preparing the gene counts dataframe as attribute (spatial_obj.genecounts_df)."
                )
                data_loader.genecounts_df = genecounts_array_to_df(
                    gene_counts_array=data_loader.load_counts_cohort(),
                    gene_names_list=data_loader.common_genes,
                    spots_metadata=spots_info,
                )
            if load_deconvolution:
                print("Preparing the celltype dataframe (spatial_obj.celltype_df)")
                data_loader.celltype_df = celltype_array_to_df(
                    celltype_array=data_loader.load_deconvolution_cohort(frac=False),
                    celltypes_name_list=data_loader.common_celltypes,
                    spots_metadata=spots_info,
                )
            if load_images:
                print("Preparing the images dict (spatial_obj.images)")
                data_loader.images = data_loader.load_images_cohort()
            return data_loader
        else:
            return {}

    def load_singlecell(self, exclude_idh1_mut: bool = False):
        """Load single-cell data.

        Returns
        -------
        adata: AnnData
            Single-cell adata object.
        exclude_idh1_mut : bool
            If True, removes samples from the data whose indices are in the predefined
            list of IDH1-mutated samples, by default False.
        """
        if "sc_rna" in self.sources:
            file_path = self.sources["sc_rna"].files["Single-cell RNA"].name
            print("Loading the single-cell data.")
            adata = SingleCellLoader(path_data=file_path).load_sc_zarr()
            if exclude_idh1_mut:
                print("Excluding IDH1 mutated samples if any present.")
                adata = adata[~adata.obs["orig.ident"].isin(PATIENTS_IDH1_MUT)]
            return adata
        else:
            raise ValueError("Single-cell data is not available in this data center.")


    def __repr__(self):
        return f"DataCenter <{self.name}> with {len(self.sources)} sources"

class DataCenters:
    """Wrapper that operates on all registered DataCenters automatically."""

    def __init__(self, *data_centers: DataCenter):
        self.data_centers = {dc.name: dc for dc in data_centers}

    def __repr__(self):
        return f"AllDataCenters({', '.join(self.data_centers)})"

    def _apply(self, method_name, *args, mode="threading", n_jobs=3, **kwargs):
        """mode='threading' or 'process'"""
        results = Parallel(n_jobs=n_jobs, backend=mode)(
            delayed(getattr(dc, method_name))(*args, **kwargs)
            for dc in self.data_centers.values()
        )
        return dict(zip(self.data_centers.keys(), results))

    def load_tabular(self, *args, **kwargs):
        return self._apply("load_tabular", *args, **kwargs)

    def load_singlecell(self, *args, **kwargs):
        return sc.concat(self._apply("load_singlecell", *args, **kwargs),label='center')

    def load_spatial(self, *args, **kwargs):
        return self._apply("load_spatial", *args, **kwargs)

    def load_spatial_adata(self,
        exclude_low_quality_visium=True,
        exclude_idh1_mut=True,
        list_sample_id=None,
        normalization="log_norm_count",
        radius=1,
        deconvolution_slot_name="deconv_cohort_level_4",
        load_images = False,
        load_deconvolution = True,
        ):
        """Load entire Visium cohort as anndata."""

        # Load all data
        loaded_objects = self.load_spatial(
            exclude_low_quality_visium=exclude_low_quality_visium,
            exclude_idh1_mut=exclude_idh1_mut,
            list_sample_id=list_sample_id,
            normalization=normalization,
            radius=radius,
            deconvolution_slot_name=deconvolution_slot_name,
            load_images=load_images,
            load_deconvolution=load_deconvolution
        )

        # Concatenate the dataframes
        genes = [obj.genecounts_df.columns for obj in loaded_objects.values()]
        common_genes = list(set(genes[0]).intersection(*genes[1:]))
        counts_df = pd.concat([obj.genecounts_df[common_genes] for obj in loaded_objects.values()])

        # Flatten the MultiIndex for AnnData compatibility
        flat_index = counts_df.index.map(lambda idx: f"{idx[0]}_{idx[1]}")

        if load_images:
            images_uns = {k:v.images['uns_spatial'] for k,v in loaded_objects.items()}
            images_obsm = {k:v.images['obsm_spatial'] for k,v in loaded_objects.items()}


        if load_deconvolution:
            obs_df = pd.concat([obj.celltype_df for obj in loaded_objects.values()])

            # Align obs to counts and create the 'sample' column
            obs_df = obs_df.loc[counts_df.index]
            obs_df['sample'] = obs_df.index.get_level_values(0)
        else:
            obs_df = pd.DataFrame(index=flat_index)

        counts_df.index = flat_index
        obs_df.index = flat_index

        del loaded_objects
        gc.collect()

        # Create the final AnnData object
        adata = sc.AnnData(
            obs=obs_df,
            var=pd.DataFrame(index=counts_df.columns),
            layers={normalization:counts_df.values},
        )

        if load_images:
            adata.uns['images_uns'] = images_uns
            adata.uns['images_obsm'] = images_obsm


        del counts_df, obs_df
        gc.collect()

        return adata
