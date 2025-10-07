"""Script that allows to derive a pandas DataFrame from a spatial object."""

# import necessary packages

from typing import Union

import numpy as np
import pandas as pd


def celltype_array_to_df(
    celltype_array: Union[dict, np.ndarray],
    celltypes_name_list: list,
    spots_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Transform array of cell abundance or fraction to a DataFrame.

    Part of the API to load the data.

    Parameters
    ----------
    celltype_array : dict or np.ndarray
        Dictionary or numpy array containing the celltype array for each sample. Can be
        the abundance or the fraction of each celltype in each spot for each sample.
    celltypes_name_list : list
        List containing the celltypes names, which is the layer of the spatial object
        (eg. spatial_object.common_celltypes)
    spots_metadata : pd.DataFrame
        DataFrame containing the metadata of each spot, which is the layer
        "coords" of the spatial object (eg. spatial_object.coords)

    Returns
    -------
    cellabun_coord_df : pd.DataFrame
        DataFrame containing the cell abundance or fraction for each celltype in each spot
        for each sample, with the coordinates of the spots.
    """
    # set the index of the spots_metadata DataFrame
    print("SPOTS METADATA COLUMNS", spots_metadata.columns)
    spots_metadata = spots_metadata.rename(columns={"sample_id": "sample"})
    spots_metadata = spots_metadata.set_index("sample")
    # if celltype_array is a dictionary, then it is the full cohort
    if isinstance(celltype_array, dict):
        # loop over the celltype_array dictionary to create a DataFrame
        for k, v in celltype_array.items():
            celltype_array[k] = pd.DataFrame(
                v, columns=celltypes_name_list, index=spots_metadata.loc[k, "barcode"]
            )
        celltype_df = pd.concat(
            celltype_array.values(),
            keys=celltype_array.keys(),
            names=["sample", "barcode"],
        )
    # if celltype_array is a numpy array, it means only a sample is loaded
    else:
        celltype_df = pd.DataFrame(
            celltype_array, columns=celltypes_name_list, index=spots_metadata["barcode"]
        )

    spots_metadata = spots_metadata.reset_index()
    spots_metadata = spots_metadata.set_index(["sample", "barcode"])
    coordinates = spots_metadata[["imagerow", "imagecol"]]

    # merge the celltype_df and the coordinates DataFrame
    cellabun_coord_df = celltype_df.merge(
        coordinates, left_index=True, right_index=True
    )

    return cellabun_coord_df


def genecounts_array_to_df(
    gene_counts_array: Union[dict, np.ndarray],
    gene_names_list: list,
    spots_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """
    Transform array of gene counts to a DataFrame.

    Part of the API to load the data.

    Parameters
    ----------
    gene_counts_array : dict or np.ndarray
        Dictionary or numpy array containing the gene counts for each sample.
    gene_names_list : list
        List containing the gene names, which is the layer common_genes of the
        spatial object (eg. spatial_object.common_genes)
    spots_metadata : pd.DataFrame
        DataFrame containing the metadata of each spot, which is the layer
        "coords" of the spatial object (eg. spatial_object.coords)

    Returns
    -------
    counts_coord_df : pd.DataFrame
        DataFrame containing the gene counts for each gene in each spot for each sample,
        with the coordinates of the spots.
    """
    # set the index of the spots_metadata DataFrame
    spots_metadata = spots_metadata.rename(columns={"sample_id": "sample"})
    spots_metadata = spots_metadata.set_index("sample")
    # if gene_counts_array is a dictionary, then it is the full cohort
    if isinstance(gene_counts_array, dict):
        # loop over the gene_counts_array dictionary to create a DataFrame
        for k, v in gene_counts_array.items():
            gene_counts_array[k] = pd.DataFrame(
                v, columns=gene_names_list, index=spots_metadata.loc[k, "barcode"]
            )
        counts_df = pd.concat(
            gene_counts_array.values(),
            keys=gene_counts_array.keys(),
            names=["sample", "barcode"],
        )
    # if gene_counts_array is a numpy array, it means only a sample is loaded
    else:
        counts_df = pd.DataFrame(
            gene_counts_array, columns=gene_names_list, index=spots_metadata["barcode"]
        )

    spots_metadata = spots_metadata.reset_index()
    spots_metadata = spots_metadata.set_index(["sample", "barcode"])
    coordinates = spots_metadata[["imagerow", "imagecol"]]

    # merge the counts_df and the coordinates DataFrame
    counts_coord_df = counts_df.merge(coordinates, left_index=True, right_index=True)

    return counts_coord_df
