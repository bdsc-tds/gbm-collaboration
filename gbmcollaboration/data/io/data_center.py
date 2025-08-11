"""Define the DataCenter class."""

from __future__ import annotations

import functools
from typing import NamedTuple, Optional

from gbmcollaboration.data.io.data_structure import DataSource
from gbmcollaboration.data.io.loaders import SingleCellLoader, VisiumCohortLoader
from gbmcollaboration.definitions import (
    BAD_QUALITY_CHUV_VISIUM,
    BAD_QUALITY_UKER_VISIUM,
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

    def load_tabular(self, exclude_idh1_mut: bool = False):
        """Load tabular data.

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

        def exclude_idh1mut(df):
            # remove the suffix from visium
            mutated_patient = PATIENTS_IDH1_MUT.copy()
            mutated_patient = [s[:-5] for s in mutated_patient]
            # convert PATIENTS_IDH1_MUT to a set for faster lookups
            patients_set = set(mutated_patient)
            # create a pattern to match any of the patient IDs in the set
            pattern = "|".join(patients_set)

            # check if any column or index contains the pattern
            column_matches = df.columns[df.columns.str.contains(pattern)]
            index_matches = df.index[df.index.str.contains(pattern)]

            # drop matched columns or index rows
            if not column_matches.empty:
                print(
                    f"IDH1mut: excluding {len(column_matches)} columns from tabular data."
                )
                df = df.drop(columns=column_matches)
            elif not index_matches.empty:
                print(
                    f"IDH1mut: excluding {len(index_matches)} rows from tabular data."
                )
                df = df.drop(index=index_matches)
            else:
                print("IDH1 mutation samples not found in the columns or index.")

            return df

        return {
            source_name: {
                file_name: (
                    exclude_idh1mut(file_.load()) if exclude_idh1_mut else file_.load()
                )
                for file_name, file_ in source.files.items()  # Only call items() if source.files is a dict
            }
            if isinstance(source.files, dict)
            else (
                exclude_idh1mut(source.files.load())
                if exclude_idh1_mut
                else source.files.load()
            )
            # Load only the sources that are clinical, bulk_rna
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
        deconvolution_slot_name: str = "deconv_cohort_level_4_2steps",
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

        Returns
        -------
        dict[str, pd.DataFrame]
            Dictionary of dataframes for each source that are spatial.
        """
        source = self.sources.get("spatial")

        # initialize ignore_sample_ids
        ignore_sample_ids = []

        if self.name == "mosaic_uker" and exclude_low_quality_visium:
            len_bad_quality = len(BAD_QUALITY_UKER_VISIUM)
            print(f"Excluding {len_bad_quality} bad quality samples at MOSAIC-UKER")
            ignore_sample_ids = BAD_QUALITY_UKER_VISIUM.copy()
        elif self.name == "mosaic_chuv" and exclude_low_quality_visium:
            len_bad_quality = len(BAD_QUALITY_CHUV_VISIUM)
            print(f"Excluding {len_bad_quality} bad quality samples at MOSAIC-CHUV")
            ignore_sample_ids = BAD_QUALITY_CHUV_VISIUM.copy()

        if exclude_idh1_mut:
            ignore_sample_ids += PATIENTS_IDH1_MUT
            print("Excluding IDH1 mutated samples if any present.")

        if source:
            # instantiate VisiumDataLoader
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
            print(
                "Preparing the gene counts dataframe as attribute (spatial_obj.genecounts_df)."
            )
            data_loader.genecounts_df = genecounts_array_to_df(
                gene_counts_array=data_loader.load_counts_cohort(),
                gene_names_list=data_loader.common_genes,
                spots_metadata=spots_info,
            )
            # print("Preparing the celltype dataframe (spatial_obj.celltype_df)")
            data_loader.celltype_df = celltype_array_to_df(
                celltype_array=data_loader.load_deconvolution_cohort(frac=False),
                celltypes_name_list=data_loader.common_celltypes,
                spots_metadata=spots_info,
            )
            return data_loader
        else:
            return {}

    def load_singlecell(self):
        """Load single-cell data.

        Returns
        -------
        adata: AnnData
            Single-cell adata object.
        """
        if "sc_rna" in self.sources:
            file_path = self.sources["sc_rna"].files["Single-cell RNA"].name
            print("Loading the single-cell data.")
            adata = SingleCellLoader(path_data=file_path).load_sc_zarr()
            return adata
        else:
            raise ValueError("Single-cell data is not available in this data center.")

    def __repr__(self):
        return f"DataCenter <{self.name}> with {len(self.sources)} sources"
