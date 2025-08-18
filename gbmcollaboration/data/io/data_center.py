"""Define the DataCenter class."""

from __future__ import annotations

import re
import functools
from typing import NamedTuple, Optional

from gbmcollaboration.data.io.data_structure import DataSource
from gbmcollaboration.data.io.loaders import SingleCellLoader, VisiumCohortLoader
from gbmcollaboration.definitions import (
    BAD_QUALITY_CHUV_VISIUM,
    BAD_QUALITY_UKER_VISIUM,
    PATIENTS_IDH1_MUT,
    UKER_MAX_ID,
    CHUV_MAX_ID
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

        def exclude_over_max_id(df: pd.DataFrame) -> pd.DataFrame:
            def keep_label(label: str) -> bool:
                # grab the number at the end
                m = re.search(r"(\d+)$", label)
                if not m:
                    return True  # keep if no number found
                sid = int(m.group(1))
                if "ER" in label:
                    return sid <= UKER_MAX_ID
                else:  # assume CH
                    return sid <= CHUV_MAX_ID

            # apply filter to index
            keep_index = [keep_label(str(idx)) for idx in df.index]
            df = df.loc[keep_index]

            # also check columns in case samples are columns
            keep_cols = [keep_label(str(col)) for col in df.columns]
            df = df.loc[:, keep_cols]

            return df

        def exclude_idh1mut(df: pd.DataFrame) -> pd.DataFrame:
            mutated_patient = [s[:-5] for s in PATIENTS_IDH1_MUT]
            patients_set = set(mutated_patient)
            pattern = "|".join(map(re.escape, patients_set))

            column_matches = df.columns[df.columns.astype(str).str.contains(pattern)]
            index_matches = df.index[df.index.astype(str).str.contains(pattern)]

            if len(column_matches) > 0:
                print(f"IDH1mut: excluding {len(column_matches)} columns from tabular data.")
                df = df.drop(columns=column_matches)
            elif len(index_matches) > 0:
                print(f"IDH1mut: excluding {len(index_matches)} rows from tabular data.")
                df = df.drop(index=index_matches)
            return df

        def load_and_filter(file_obj):
            df = file_obj.load()
            df = exclude_over_max_id(df)
            if exclude_idh1_mut:
                df = exclude_idh1mut(df)
            return df

        return {
            source_name: {
                file_name: load_and_filter(file_)
                for file_name, file_ in source.files.items()
            } if isinstance(source.files, dict) else load_and_filter(source.files)
            for source_name, source in self.sources.items()
            if source_name in ["clinical", "bulk_rna", "wes"]
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
