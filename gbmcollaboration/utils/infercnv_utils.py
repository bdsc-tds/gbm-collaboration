"""Script with utility functions to use with inferCNVpy."""


def add_reference_column(
    single_cell_obj,
    cell_type_annotation_key,
    celltypes_to_reference,
    threshold_n_genes=400,
    threshold_n_counts=10000,
):
    """Add the reference column to the object.

    The reference column will be used in the reference_key argument of the inferCNVpy function.

    Parameters
    ----------
    single_cell_obj : anndata object
        Annotated data object.
    cell_type_annotation_key : str
        Key in the obs to use as the cell type annotation for mapping.
    celltypes_to_reference : dict
        Dictionary mapping the cell types to the reference.
    threshold_n_genes : int, optional
        Threshold to consider a cell as low_genes, by default 400.
    threshold_n_counts : int, optional
        Threshold to consider a cell as high_counts, by default 10000.

    Returns
    -------
    anndata object
        Annotated data object with the reference column.
    """
    # create a "reference" column in the obs to filter the healthy reference cells
    # use the mapping provided by the user
    single_cell_obj.obs["reference"] = single_cell_obj.obs[
        cell_type_annotation_key
    ].map(celltypes_to_reference)

    # calculate n_genes and n_counts to further remove cells from the reference
    single_cell_obj.obs["n_genes"] = (
        single_cell_obj.layers["unnormalized_counts"] > 0
    ).sum(axis=1)
    single_cell_obj.obs["n_counts"] = single_cell_obj.layers["unnormalized_counts"].sum(
        axis=1
    )

    # now annotate the cells accordingly if the conditions are met
    # therefore they will not be used as the "healthy" reference
    single_cell_obj.obs["reference"][
        (single_cell_obj.obs["n_genes"] < threshold_n_genes)
    ] = "low_genes"
    single_cell_obj.obs["reference"][
        (single_cell_obj.obs["n_counts"] > threshold_n_counts)
    ] = "high_counts"

    return single_cell_obj
