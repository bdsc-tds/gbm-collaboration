"""Pre-processor of tabular data."""

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd

from gbmcollaboration.definitions import COL_DATES, KEYWORDS_COL_TO_REMOVE
from gbmcollaboration.utils.dates import convert_to_datetime, time_to_event


# Opt-in to the future behavior
pd.set_option("future.no_silent_downcasting", True)


class BaseDataCleaner(ABC):
    """Clean Tabular Data. Class to remove typos, preprocess data columns.

    Attributes
    ----------
    dataset_name : str

    Methods
    -------
    transform: Transform the input data.
    """

    @abstractmethod
    def transform(self, x):
        """Transform the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        pass


class MOSAICLightClinicalDataCleaner(BaseDataCleaner):
    """Add Subject Id as index and remove columns with only nan values.

    Useful only for the original clinical data from MOSAIC.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.set_index("Subject Id", inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICClinicalDataCleaner(BaseDataCleaner):
    """Clean processed clinical data from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.set_index("Subject Id", inplace=True)

        # Remove columns specific to other TAs
        # Remove columns if the column name contains a keyword in KEYWORDS_COL_TO_REMOVE
        for keyword in KEYWORDS_COL_TO_REMOVE:
            x = x.loc[:, ~x.columns.str.contains(keyword, case=False)]

        # Replace Unknown values by nan
        x.replace("Unknown", np.nan, inplace=True)

        # Convert date columns to datetimes format
        for col in COL_DATES:
            if col in x.columns:
                convert_to_datetime(x, col, dayfirst=False)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        return x


class MOSAICEndpointsDataCleaner(MOSAICClinicalDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x = super().transform(x)

        x = x[
            [
                "Date of cancer diagnosis",
                "Date of death",
                "Date of last follow-up before inclusion",
                "Vital status at inclusion",
            ]
        ]
        x = x.dropna(how="all")

        # Remove duplicated rows
        x = x[~x.index.duplicated(keep="first")]

        # Postprocess survival dates and status
        x["is_alive"] = x["Vital status at inclusion"] == "Alive"

        # Compute the OS
        x["OS"] = time_to_event(
            x,
            col_start="Date of cancer diagnosis",
            col_end="Date of death",
            col_censored="is_alive",
            col_end_censored="Date of last follow-up before inclusion",
            unit="years",
        )

        x = x.drop(
            columns=[
                "Date of death",
                "Date of last follow-up before inclusion",
                "Vital status at inclusion",
                "is_alive",
            ]
        )

        return x


class MOSAICKeyEventsDataCleaner(BaseDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.set_index("Subject Id", inplace=True)

        # Replace Unknown values by nan
        x.replace("Unknown", np.nan, inplace=True)

        # Remove columns with only nan values
        x.dropna(axis=1, how="all", inplace=True)

        # Convert date columns to datetime
        convert_to_datetime(x, col="date", dayfirst=False)

        return x


class MOSAICBulkRNADataCleaner(BaseDataCleaner):
    """Clean Key Events data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.rename(columns={"gene": "EnsemblID"}, inplace=True)
        x.set_index("EnsemblID", inplace=True)

        # If weird ensembl IDs exist, drop them
        # Appears in the deseq normalized counts for instance
        x = x[~x.index.str.contains("g@-ext")]  # corresponds to
        # immunoglobulin super-loci, mainly used to find fusion genes
        # and therefore not necessary for DGEA

        return x


class MOSAICWESDataCleaner(BaseDataCleaner):
    """Clean WES data file from MOSAIC.

    Remove typos, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x = x.dropna(subset=["gene_name"])
        x = x.set_index("gene_name")
        x = x.drop(columns=["ensembl_gene_id"])
        x = x.T  # Transpose the data to have patient id in index
        x.index.name = "Subject Id"

        return x


class MOSAICHEDataCleaner(BaseDataCleaner):
    """Clean HE data file from MOSAIC.

    Add patient ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "Subject Id"] = x["path"].apply(lambda v: v.name[:9])
        x.loc[:, "patient id"] = x["Subject Id"].apply(lambda v: v[:-1])
        x = x.set_index("Subject Id")
        return x


class MOSAICHEFeaturesDataCleaner(BaseDataCleaner):
    """Clean HE features data file from MOSAIC.

    Add patient ID, preprocess data columns.

    Methods
    -------
    transform: Transform the input data.
    """

    def transform(self, x):
        """Transform the input data.

        Note that this method acts **inplace** on the input data.

        Parameters
        ----------
        x : pd.DataFrame
            Input data

        Returns
        -------
        pd.DataFrame
            x with no typos, parsed dates
        """
        x.loc[:, "Subject Id"] = x["path"].apply(lambda v: v.parent.name[:9])
        x.loc[:, "patient id"] = x["Subject Id"].apply(lambda v: v[:-1])
        x.loc[:, "mask_path"] = x["path"].apply(lambda v: v.parent / "mask.png")
        x.loc[:, "info_path"] = x["path"].apply(
            lambda v: v.parent / "tiling_tool_format.json"
        )
        x = x.set_index("Subject Id")
        return x
