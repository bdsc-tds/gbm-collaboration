"""Module for loading data."""

import anndata as ad
import pandas as pd
from gbmcollaboration.data.dataset_classes.visium_cohort import VisiumCohort
from warnings import warn
import s3fs
from collections.abc import Iterable


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
        data_transformer : optional
            An instance of DataCleaner or its subclass for data
            preprocessing, default is None.
        load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        # if data_transformer is None the data will be returned in the same
        # format as loaded
        super().__init__(data_transformer=data_transformer, **load_data_kwargs)

    def load_data(self, data_path):
        """
        Load and optionally transform data.

        Parameters
        ----------
        data_path : str or list of str
            Path or list of paths to CSV file(s).

        Returns
        -------
        data : DataFrame or transformed output
            Loaded and (optionally) transformed data.
        """
        # If data_path is a list, pass directly to transformer
        if isinstance(data_path, Iterable) and not isinstance(data_path, (str, bytes)):
            if self.data_transformer is not None:
                return self.data_transformer.transform(data_path)
            else:
                raise ValueError("A data_transformer must be provided for list of paths.")

        # Else load a single CSV
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


class Cell2LocationLoader(BaseDataLoader):
    """Class to load cell2location results from h5ad files.
    
    Loads cell2location AnnData objects and extracts cell count estimates.
    Specifically extracts q50_cell_abundance_w_sf values and computes total cells per spot.
    """
    
    def __init__(self, level: int = 2, data_transformer=None, **load_data_kwargs):
        """Initialize the Cell2LocationLoader object.
        
        Parameters
        ----------
        level : int
            Cell type annotation level to load (2, 4, etc.). Default: 2
        data_transformer : optional
            An instance of DataCleaner for data preprocessing, default is None.
        **load_data_kwargs : dict
            Additional keyword arguments for the load_data method.
        """
        super().__init__(data_transformer=data_transformer, **load_data_kwargs)
        self.level = level
    
    def load_data(self, data_path):
        """Load cell2location h5ad file and extract cell abundance estimates.
        
        Parameters
        ----------
        data_path : str or Path
            Path to directory containing cell2location h5ad files.
            Expected structure: {data_path}/{sample_id}/cohort_level_{level}.h5ad
            
        Returns
        -------
        pd.DataFrame
            DataFrame with cell abundance estimates per spot.
            Columns: cell types, Index: spot barcodes
            Also includes 'total_cells' column with sum over all celltypes per spot.
        """
        from pathlib import Path
        from cloudpathlib import S3Path
        import s3fs
        
        # Convert to Path-like object
        if isinstance(data_path, str):
            if data_path.startswith('s3://'):
                base_path = S3Path(data_path)
            else:
                base_path = Path(data_path)
        else:
            base_path = data_path
        
        # Find all sample directories
        sample_dirs = [d for d in base_path.iterdir() if d.is_dir()]
        
        all_results = []
        
        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            h5ad_file = sample_dir / f"cohort_level_{self.level}.h5ad"
            
            if not h5ad_file.exists():
                warn(
                    f"Cell2location file not found for {sample_id}: {h5ad_file}",
                    UserWarning,
                )
                continue
            
            # Load h5ad file
            try:
                if isinstance(h5ad_file, S3Path):
                    # Use s3fs for S3 paths
                    s3 = s3fs.S3FileSystem(anon=False)
                    # S3Path needs to be converted to string for s3fs
                    h5ad_path_str = str(h5ad_file).replace('s3://', '')
                    with s3.open(h5ad_path_str, 'rb') as f:
                        adata = ad.read_h5ad(f)
                else:
                    adata = ad.read_h5ad(h5ad_file)
                
                # Extract q50_cell_abundance_w_sf from obsm
                if 'q50_cell_abundance_w_sf' not in adata.obsm:
                    warn(
                        f"q50_cell_abundance_w_sf not found in {sample_id} AnnData object. "
                        f"Available keys: {list(adata.obsm.keys())}",
                        UserWarning,
                    )
                    continue
                
                # Get cell abundance DataFrame (spots × celltypes)
                # Try to get cell type names from various possible locations
                celltype_names = None
                if 'mod' in adata.uns and 'factor_names' in adata.uns['mod']:
                    celltype_names = adata.uns['mod']['factor_names']
                elif 'factor_names' in adata.uns:
                    celltype_names = adata.uns['factor_names']
                elif 'cell_type_names' in adata.uns:
                    celltype_names = adata.uns['cell_type_names']
                
                if celltype_names is None:
                    # Fallback: create generic names
                    n_celltypes = adata.obsm['q50_cell_abundance_w_sf'].shape[1]
                    celltype_names = [f'celltype_{i}' for i in range(n_celltypes)]
                
                cell_abundance = pd.DataFrame(
                    adata.obsm['q50_cell_abundance_w_sf'],
                    index=adata.obs_names,
                    columns=celltype_names
                )
                
                # Add sample ID to index if not already present
                if 'sample' not in cell_abundance.index.names:
                    cell_abundance.index = pd.MultiIndex.from_tuples(
                        [(sample_id, idx) for idx in cell_abundance.index],
                        names=['sample', 'barcode']
                    )
                
                # Compute total cells per spot (sum over all celltypes)
                cell_abundance['total_cells'] = cell_abundance.sum(axis=1)
                
                all_results.append(cell_abundance)
                
            except Exception as e:
                warn(
                    f"Failed to load cell2location data from {h5ad_file}: {str(e)}",
                    UserWarning,
                )
                continue
        
        if not all_results:
            warn(
                "No cell2location data loaded. Returning empty DataFrame.",
                UserWarning,
            )
            return pd.DataFrame()
        
        # Concatenate all samples
        result_df = pd.concat(all_results, axis=0)
        
        # Apply transformer if provided
        if self.data_transformer is not None:
            result_df = self.data_transformer.transform(result_df)
        
        return result_df
