"""Script to plot various variable on Visium image."""

from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_celltype_abun_visium(
    df: pd.DataFrame,
    sample_id: str,
    cell_type: str,
    colorbar: bool = True,
    ax: Optional[matplotlib.axes.Axes] = None,
):
    """Plot the cell type abundance on Visium for a given patient.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the Visium data and metadata with the imagecol
        and imagerow columns from the coordinates layers of the spatial object
        and the cell type abundance columns.
    sample_id : str
        The sample ID of the patient to plot.
    cell_type : str
        The cell type to plot. (eg. "Mural cell")
    colorbar: bool
        If True, add a colorbar legend to the heatmap plot.
    ax: matplotlib.axes.Axes
        (Optional) If given, use this matplotlib axis to plot the confusion matrix.
        Otherwise create one and return the matplotlib figure associated with this ax.
    """
    df_filt = df.copy()
    patient_target = df_filt.loc[sample_id]
    if ax is None:
        _, ax = plt.subplots()
    ax.scatter(
        patient_target["imagecol"],
        patient_target["imagerow"],
        c=patient_target[cell_type],
        s=10,
    )
    ax.set_title(f"{cell_type} abundance in {sample_id}")
    # Display figure with colorbar
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(ax.collections[0], cax=cax, orientation="vertical")


def plot_gene_expression_visium(
    df,
    sample_id,
    gene,
    colorbar: bool = True,
    ax: Optional[matplotlib.axes.Axes] = None,
):
    """Plot the gene expression on Visium for a given patient.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the Visium data and metadata with the imagecol
        and imagerow columns from the coordinates layers of the spatial object
        and the gene expression columns.
    sample_id : str
        The sample ID of the patient to plot.
    gene : str
        The gene to plot for expression. (eg. EGFR, CD34)
    colorbar: bool
        If True, add a colorbar legend to the heatmap plot.
    ax: matplotlib.axes.Axes
        (Optional) If given, use this matplotlib axis to plot the confusion matrix.
        Otherwise create one and return the matplotlib figure associated with this ax.
    """
    df_filt = df.copy()
    patient_target = df_filt.loc[sample_id]
    if ax is None:
        _, ax = plt.subplots()
    ax.scatter(
        patient_target["imagecol"],
        patient_target["imagerow"],
        c=patient_target[gene],
        s=10,
    )
    ax.set_title(f"{gene} expression in {sample_id}")
    # Display figure with colorbar
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(ax.collections[0], cax=cax, orientation="vertical")
