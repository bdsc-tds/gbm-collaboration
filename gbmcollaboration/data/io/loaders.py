"""Module for loading data."""

import anndata as ad
import pandas as pd
from gbmcollaboration.data.dataset_classes.visium_cohort import VisiumCohort
from warnings import warn
import s3fs


class BaseDataLoader:
    """Base class for loader."""

    def __init__(self, data_transformer=None, **load_data_kwargs):
        """Initialize base class.

        Parameters
        ----------
        data_transformer : optional
            An instance of DataCleaner or its subclass for data
            preprocessing, by default None
        **load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        self.data_transformer = data_transformer
        self.load_data_kwargs = load_data_kwargs

    def load_data(self, data_path):
        """Load and transform the input data.

        Parameters
        ----------
        data_path : str
            Path to the data file to be read.

        Returns
        -------
        data : DataFrame
            Loaded and transformed data.
        """
        pass


class XLSXDataLoader(BaseDataLoader):
    """A class for loading Excel data.

    This transformer loads and preprocesses data using pandas.
    """

    def __init__(self, data_transformer=None, **load_data_kwargs):
        """Initialize the CSVDataLoader object.

        Parameters
        ----------
        data_transformer : optional
            An instance of DataCleaner for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        # if data_transformer is None the data will be returned in the same
        # format as loaded

        super().__init__(data_transformer=data_transformer, **load_data_kwargs)

    def load_data(self, data_path):
        """Load and transform the input data.

        Parameters
        ----------
         data_path : str
            Path to the Excel file to be read.

        Returns
        -------
        X_csv : DataFrame
            Loaded and transformed data not including the target columns
        y : DataFrame
            Data corresponding to target columns with index as patient ids.
        """
        data = pd.read_excel(data_path, **self.load_data_kwargs)
        if self.data_transformer is not None:
            data = self.data_transformer.transform(data)
        return data


class CSVDataLoader(BaseDataLoader):
    """A class for loading CSV data.

    This transformer loads and preprocesses data using pandas.
    """

    def __init__(self, data_transformer=None, **load_data_kwargs):
        """Initialize the CSVDataLoader object.

        Parameters
        ----------
        dataset_name : str, optional
            Dataset name,
        data_transformer : optional
            An instance of DataCleaner or its subclass for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        # if data_transformer is None the data will be returned in the same
        # format as loaded
        super().__init__(data_transformer=data_transformer, **load_data_kwargs)

    def load_data(
        self,
        data_path,
    ):
        """Transform the input data using pandas.

        Returns
        -------
        X_csv : DataFrame
            Loaded and transformed data not including the target columns
        y : DataFrame
            Data corresponding to target columns with index as patient ids.
        """
        data = pd.read_csv(filepath_or_buffer=data_path, **self.load_data_kwargs)
        if self.data_transformer is not None:
            data = self.data_transformer.transform(data)
        return data


class PathsLoader(BaseDataLoader):
    """Class to load a data frame with the path to the paths of files in a directory."""

    def __init__(
        self, pattern_files: str = "*", data_transformer=None, **load_data_kwargs
    ):
        """Initialize the PathsLoader object.

        Parameters
        ----------
        pattern_files : str
            Pattern to match the files in the directory.
        data_transformer : optional
            An instance of DataCleaner for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        super().__init__(data_transformer=data_transformer, **load_data_kwargs)
        self.pattern_files = pattern_files

    def load_data(self, data_path):
        """Load and transform the input data.

        Parameters
        ----------
        data_path : str
            Path to the folder containing the files to be included.

        Returns
        -------
        data : DataFrame
            Loaded and transformed data.
        """
        paths = list(data_path.glob(self.pattern_files))
        data = pd.DataFrame({"path": paths})
        if self.data_transformer is not None:
            data = self.data_transformer.transform(data)
        return data


class VisiumCohortLoader(VisiumCohort):
    """Class to load Visium data."""

    def __init__(
        self,
        path_dataset,
        ignore_sample_ids,
        list_sample_ids,
        normalization,
        radius,
        deconvolution_slot_name,
    ):
        """Initialize the VisiumCohortLoader object.

        Parameters
        ----------
        path_dataset : list[str]
            Paths to the data.
        ignore_sample_ids : list[str] | None
            List of sample ids to ignore in case of loading the full cohort.
            Useful when you want to load all samples but avoid errors due to QC issues.
            By default for CHUV, we will exclude samples because of QC issues.
        list_sample_id : list[str]
            List of sample IDs to load. If None, load the full cohort.
            Default is None.
        normalization : str
            Normalization method. Available methods are "raw_count",
            "log_norm_count", "norm_count", "SCTransform".
        radius : int
            Radius for the spatial data.
        deconvolution_slot_name : str
            Slot name for the deconvolution data.
        """
        super().__init__(
            path_dataset=path_dataset,
            ignore_sample_ids=ignore_sample_ids,
            list_sample_ids=list_sample_ids,
            normalization=normalization,
            radius=radius,
            deconvolution_slot_name=deconvolution_slot_name,
        )

    def load_data(self):
        """Load the counts for the cohort.

        Returns
        -------
        The VisiumCohortLoader object.
        To access the counts, add self.load_counts_cohort()
        """
        return self


class SingleCellLoader:
    """Class to load the single cell anndata."""

    def __init__(self, path_data):
        """Initialize the SingleCellLoader object."""
        self.path_data = path_data

    def load_sc_zarr(self):
        try:
            # the use of S3FileSystem improve caching
            s3 = s3fs.S3FileSystem(anon=False)  # Or anon=True if public
            store = s3fs.S3Map(root=self.path_data, s3=s3)
            adata = ad.read_zarr(store)
            return adata
        except ValueError as e:
            warn(
                f"Failed to load data from {self.path_data}: {str(e)}. "
                "Returning empty AnnData object.",
                UserWarning,
            )
