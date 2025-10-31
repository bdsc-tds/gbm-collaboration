"""A module for building the data centers loaders."""

from __future__ import annotations

from attrs import frozen

from gbmcollaboration.data.io.data_center import DataCenter, DataCenters
from gbmcollaboration.data.io.data_structure import (
    DataFile,
    DataSource,
)
from gbmcollaboration.data.io.loaders import (
    CSVDataLoader,
    PathsLoader,
    XLSXDataLoader,
)
from gbmcollaboration.data.io.preprocessor_tabular import (
    MOSAICBulkRNADataCleaner,
    MOSAICClinicalDataCleaner,
    MOSAICEndpointsDataCleaner,
    MOSAICHEDataCleaner,
    MOSAICHEFeaturesDataCleaner,
    MOSAICKeyEventsDataCleaner,
    MOSAICLightClinicalDataCleaner,
    MOSAICWESDataCleaner,
    MOSAICWESTMBDataCleaner
)
from gbmcollaboration.definitions import (
    DATA_PATH_CHUV_FILES,
    DATA_PATH_UKER_FILES,
    DATA_PATH_CR_FILES,
    MOSAIC_DATA_DICT_FILENAME,
)


@frozen
class SourceNames:
    """Define the names of the sources.

    This class is used to avoid typos when accessing the sources of a center.
    """

    clinical: str = "clinical"
    data_dict: str = "data_dict"
    wes: str = "wes"
    bulk_rna: str = "bulk_rna"
    sc_rna: str = "sc_rna"
    spatial: str = "spatial"
    he: str = "he"

    # force default values by forbidding arguments during initialization.
    def __init__(self):
        self.__attrs_init__()  # type: ignore # pylint: disable=no-member


SOURCES = SourceNames()

# CHUV data
CHUVDataCenter = DataCenter(
    name="mosaic_chuv",
    sources={
        SOURCES.data_dict: DataSource(
            files={
                "Data dict": DataFile(
                    name=MOSAIC_DATA_DICT_FILENAME,
                    loader=XLSXDataLoader(),
                ),
            }
        ),
        SOURCES.clinical: DataSource(
            files={
                "original clinical": DataFile(
                    name=DATA_PATH_CHUV_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICLightClinicalDataCleaner(),
                    ),
                ),
                "processed gbm clinical": DataFile(
                    name=DATA_PATH_CHUV_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICClinicalDataCleaner(),
                    ),
                ),
                "endpoints": DataFile(
                    name=DATA_PATH_CHUV_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICEndpointsDataCleaner(),
                    ),
                ),
                "key events": DataFile(
                    name=DATA_PATH_CHUV_FILES["Key events"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICKeyEventsDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.bulk_rna: DataSource(
            files={
                "raw counts": DataFile(
                    name=DATA_PATH_CHUV_FILES["Bulk RNA raw counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "TPM": DataFile(
                    name=DATA_PATH_CHUV_FILES["Bulk RNA TPM"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "normalized counts": DataFile(
                    name=DATA_PATH_CHUV_FILES["Bulk RNA normalized counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
            }
        ),
        SOURCES.spatial: DataSource(
            files={
                "cohort zarr files": DataFile(
                    name=DATA_PATH_CHUV_FILES["Spatial Zarr Files"],
                    loader=None,
                ),
            }
        ),
        SOURCES.sc_rna: DataSource(
            files={
                "Single-cell RNA": DataFile(
                    name=DATA_PATH_CHUV_FILES["Single-cell RNA"],
                    loader=None,
                ),
            }
        ),
        SOURCES.wes: DataSource(
            files={
                "WES CNV deletion": DataFile(
                    name=DATA_PATH_CHUV_FILES["WES CNV deletion"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES CNV amplification": DataFile(
                    name=DATA_PATH_CHUV_FILES["WES CNV amplification"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES mutations": DataFile(
                    name=DATA_PATH_CHUV_FILES["WES mutations"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES TMB": DataFile(
                    name=DATA_PATH_CHUV_FILES["WES TMB"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESTMBDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.he: DataSource(
            files={
                "HE images": DataFile(
                    name=DATA_PATH_CHUV_FILES["HE_images"],
                    loader=PathsLoader(
                        pattern_files="*_eHnE.svs",
                        data_transformer=MOSAICHEDataCleaner(),
                    ),
                ),
                "H0 features": DataFile(
                    name=DATA_PATH_CHUV_FILES["HE_features_H0"],
                    loader=PathsLoader(
                        pattern_files="*/features.npy",
                        data_transformer=MOSAICHEFeaturesDataCleaner(),
                    ),
                ),
            },
        ),
    },
)

# UKER data
UKERDataCenter = DataCenter(
    name="mosaic_uker",
    sources={
        SOURCES.data_dict: DataSource(
            files={
                "Data dict": DataFile(
                    name=MOSAIC_DATA_DICT_FILENAME,
                    loader=XLSXDataLoader(),
                ),
            }
        ),
        SOURCES.clinical: DataSource(
            files={
                "original clinical": DataFile(
                    name=DATA_PATH_UKER_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICLightClinicalDataCleaner(),
                    ),
                ),
                "processed gbm clinical": DataFile(
                    name=DATA_PATH_UKER_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICClinicalDataCleaner(),
                    ),
                ),
                "endpoints": DataFile(
                    name=DATA_PATH_UKER_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICEndpointsDataCleaner(),
                    ),
                ),
                "key events": DataFile(
                    name=DATA_PATH_UKER_FILES["Key events"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICKeyEventsDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.bulk_rna: DataSource(
            files={
                "raw counts": DataFile(
                    name=DATA_PATH_UKER_FILES["Bulk RNA raw counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "TPM": DataFile(
                    name=DATA_PATH_UKER_FILES["Bulk RNA TPM"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "normalized counts": DataFile(
                    name=DATA_PATH_UKER_FILES["Bulk RNA normalized counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
            }
        ),
        SOURCES.spatial: DataSource(
            files={
                "cohort zarr files": DataFile(
                    name=DATA_PATH_UKER_FILES["Spatial Zarr Files"],
                    loader=None,
                ),
            }
        ),
        SOURCES.sc_rna: DataSource(
            files={
                "Single-cell RNA": DataFile(
                    name=DATA_PATH_UKER_FILES["Single-cell RNA"],
                    loader=None,
                ),
            }
        ),
        SOURCES.wes: DataSource(
            files={
                "WES CNV deletion": DataFile(
                    name=DATA_PATH_UKER_FILES["WES CNV deletion"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES CNV amplification": DataFile(
                    name=DATA_PATH_UKER_FILES["WES CNV amplification"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES mutations": DataFile(
                    name=DATA_PATH_UKER_FILES["WES mutations"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES TMB": DataFile(
                    name=DATA_PATH_UKER_FILES["WES TMB"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESTMBDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.he: DataSource(
            files={
                "HE images": DataFile(
                    name=DATA_PATH_UKER_FILES["HE_images"],
                    loader=PathsLoader(
                        pattern_files="*_eHnE.svs",
                        data_transformer=MOSAICHEDataCleaner(),
                    ),
                ),
                "H0 features": DataFile(
                    name=DATA_PATH_UKER_FILES["HE_features_H0"],
                    loader=PathsLoader(
                        pattern_files="*/features.npy",
                        data_transformer=MOSAICHEFeaturesDataCleaner(),
                    ),
                ),
            },
        ),
    },
)


# CR data
CRDataCenter = DataCenter(
    name="mosaic_cr",
    sources={
        SOURCES.data_dict: DataSource(
            files={
                "Data dict": DataFile(
                    name=MOSAIC_DATA_DICT_FILENAME,
                    loader=XLSXDataLoader(),
                ),
            }
        ),
        SOURCES.clinical: DataSource(
            files={
                "original clinical": DataFile(
                    name=DATA_PATH_CR_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICLightClinicalDataCleaner(),
                    ),
                ),
                "processed gbm clinical": DataFile(
                    name=DATA_PATH_CR_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICClinicalDataCleaner(),
                    ),
                ),
                "endpoints": DataFile(
                    name=DATA_PATH_CR_FILES["Processed clinical"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICEndpointsDataCleaner(),
                    ),
                ),
                "key events": DataFile(
                    name=DATA_PATH_CR_FILES["Key events"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICKeyEventsDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.bulk_rna: DataSource(
            files={
                "raw counts": DataFile(
                    name=DATA_PATH_CR_FILES["Bulk RNA raw counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "TPM": DataFile(
                    name=DATA_PATH_CR_FILES["Bulk RNA TPM"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
                "normalized counts": DataFile(
                    name=DATA_PATH_CR_FILES["Bulk RNA normalized counts"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICBulkRNADataCleaner(), sep="\t"
                    ),
                ),
            }
        ),
        SOURCES.spatial: DataSource(
            files={
                "cohort zarr files": DataFile(
                    name=DATA_PATH_CR_FILES["Spatial Zarr Files"],
                    loader=None,
                ),
            }
        ),
        SOURCES.sc_rna: DataSource(
            files={
                "Single-cell RNA": DataFile(
                    name=DATA_PATH_CR_FILES["Single-cell RNA"],
                    loader=None,
                ),
            }
        ),
        SOURCES.wes: DataSource(
            files={
                "WES CNV deletion": DataFile(
                    name=DATA_PATH_CR_FILES["WES CNV deletion"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES CNV amplification": DataFile(
                    name=DATA_PATH_CR_FILES["WES CNV amplification"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES mutations": DataFile(
                    name=DATA_PATH_CR_FILES["WES mutations"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESDataCleaner(),
                    ),
                ),
                "WES TMB": DataFile(
                    name=DATA_PATH_CR_FILES["WES TMB"],
                    loader=CSVDataLoader(
                        data_transformer=MOSAICWESTMBDataCleaner(),
                    ),
                ),
            }
        ),
        SOURCES.he: DataSource(
            files={
                "HE images": DataFile(
                    name=DATA_PATH_CR_FILES["HE_images"],
                    loader=PathsLoader(
                        pattern_files="*_eHnE.svs",
                        data_transformer=MOSAICHEDataCleaner(),
                    ),
                ),
                "H0 features": DataFile(
                    name=DATA_PATH_CR_FILES["HE_features_H0"],
                    loader=PathsLoader(
                        pattern_files="*/features.npy",
                        data_transformer=MOSAICHEFeaturesDataCleaner(),
                    ),
                ),
            },
        ),
    },
)


AllDataCenters = DataCenters(CHUVDataCenter, UKERDataCenter, CRDataCenter)
