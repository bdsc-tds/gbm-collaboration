"""Define the DataCenter class."""

from __future__ import annotations

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
    PATIENT_EXCLUSION,
    RELAPSE_EXCLUSION,
    PATIENTS_IDH1_MUT,
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

    def load_tabular(
        self,
        exclude_idh1_mut: bool = False,
        exclude_patient_ids: bool = False,
        exclude_relapse: bool = False,):
        """Load tabular data.

        Parameters
        ----------
        exclude_idh1_mut : bool
            If True, removes samples from the data whose indices are in the predefined
            list of IDH1-mutated samples, by default False.
        exclude_patient_ids : bool
            If True, removes samples from the data whose indices are in the predefined
            list of excluded patients, by default False.
        exclude_relapse : bool
            If True, removes samples from the data whose indices are in the predefined
            list of relapse samples, by default False.

        Returns
        -------
        dict[str, pd.DataFrame]
            Dictionary of dataframes for each source that are tabular.
        """

        def apply_exclusions(
            df: pd.DataFrame,
            exclude_idh1_mut: bool,
            exclude_patient_ids: bool,
            exclude_relapse: bool
        ):
            """Apply all exclusion rules to a dataframe."""

            # --- Exclude IDH1-mut patients ---
            if exclude_idh1_mut:
                patients_set = set(PATIENTS_IDH1_MUT)
                pattern = "|".join(patients_set)

                col_matches = df.columns[df.columns.str.contains(pattern)]
                idx_matches = df.index[df.index.str.contains(pattern)]

                if not col_matches.empty:
                    print(f"IDH1mut: excluding {len(col_matches)} columns.")
                    df = df.drop(columns=col_matches)
                if not idx_matches.empty:
                    print(f"IDH1mut: excluding {len(idx_matches)} rows.")
                    df = df.drop(index=idx_matches)

            # --- Exclude specific patients ---
            if exclude_patient_ids:
                patients_set = set(PATIENT_EXCLUSION)
                pattern = "|".join(patients_set)

                col_matches = df.columns[df.columns.str.contains(pattern)]
                idx_matches = df.index[df.index.str.contains(pattern)]

                if not col_matches.empty:
                    print(f"Patients: excluding {len(col_matches)} columns.")
                    df = df.drop(columns=col_matches)
                if not idx_matches.empty:
                    print(f"Patients: excluding {len(idx_matches)} rows.")
                    df = df.drop(index=idx_matches)

            return df

            # --- Exclude relapse sample ---
            if exclude_relapse:
                patients_set = set(RELAPSE_EXCLUSION)
                pattern = "|".join(patients_set)

                col_matches = df.columns[df.columns.str.contains(pattern)]
                idx_matches = df.index[df.index.str.contains(pattern)]

                if not col_matches.empty:
                    print(f"Relapse samples: excluding {len(col_matches)} columns.")
                    df = df.drop(columns=col_matches)
                if not idx_matches.empty:
                    print(f"Relapse samples: excluding {len(idx_matches)} rows.")
                    df = df.drop(index=idx_matches)

            return df

        return {
            source_name: {
                file_name: apply_exclusions(
                    file_.load(),
                    exclude_idh1_mut=exclude_idh1_mut,
                    exclude_patient_ids=exclude_patient_ids,
                    exclude_relapse=exclude_relapse,
                )
                for file_name, file_ in source.files.items()
            }
            if isinstance(source.files, dict)
            else apply_exclusions(
                source.files.load(),
                exclude_idh1_mut=exclude_idh1_mut,
                exclude_patient_ids=exclude_patient_ids,
                exclude_relapse=exclude_relapse,
            )
            for source_name, source in self.sources.items()
            if source_name in ["clinical", "bulk_rna", "wes"]
        }


    def load_spatial(
        self,
        exclude_low_quality_visium: bool = False,
        exclude_idh1_mut: bool = False,
        exclude_patient_ids: bool = False,
        exclude_relapse: bool = False,
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
        exclude_patient_ids : bool
            If True, removes samples from the data whose indices are in the predefined
            list of excluded patients, by default False.
        exclude_relapse : bool
            If True, removes samples from the data whose indices are in the predefined
            list of relapse samples, by default False.
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

        if exclude_patient_ids:
            if ignore_sample_ids is None:
                ignore_sample_ids = PATIENT_EXCLUSION
            else:
                ignore_sample_ids += PATIENT_EXCLUSION
            print("Excluding predefined patients if any present.")

        if exclude_relapse:
            if ignore_sample_ids is None:
                ignore_sample_ids = RELAPSE_EXCLUSION
            else:
                ignore_sample_ids += RELAPSE_EXCLUSION
            print("Excluding predefined patients if any present.")

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

    def load_singlecell(
        self, exclude_idh1_mut: bool = False,
        exclude_patients: bool = False,
        exclude_relapse: bool = False):
        """Load single-cell data.

        Returns
        -------
        adata: AnnData
            Single-cell adata object.
        exclude_idh1_mut : bool
            If True, removes samples from the data whose indices are in the predefined
            list of IDH1-mutated samples, by default False.
        exclude_patients : bool
            If True, removes samples from the data whose indices are in the predefined
            list of excluded patients, by default False.
        exclude_relapse : bool
            If True, removes samples from the data whose indices are in the predefined
            list of relapse samples, by default False.
        """
        if "sc_rna" in self.sources:
            file_path = self.sources["sc_rna"].files["Single-cell RNA"].name
            print("Loading the single-cell data.")
            adata = SingleCellLoader(path_data=file_path).load_sc_zarr()
            if exclude_idh1_mut:
                print("Excluding IDH1 mutated samples if any present.")
                adata = adata[~adata.obs["orig.ident"].isin(PATIENTS_IDH1_MUT)]
            if exclude_relapse:
                print("Excluding relapse samples if any present.")
                adata = adata[~adata.obs["orig.ident"].isin(RELAPSE_EXCLUSION)]
            if exclude_patients:
                # need to search for the samples that match the patient ids, since
                # the sample ids have suffixes
                PATIENT_EXCLUSION_MATCHES = []
                for patient_id in PATIENT_EXCLUSION:
                    matches = [s for s in adata.obs["orig.ident"].unique() if patient_id in s]
                    PATIENT_EXCLUSION_MATCHES.extend(matches)
                print(f"Excluding {len(PATIENT_EXCLUSION_MATCHES)} samples from predefined patients ids.")
                adata = adata[~adata.obs["orig.ident"].isin(PATIENT_EXCLUSION_MATCHES)]
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
