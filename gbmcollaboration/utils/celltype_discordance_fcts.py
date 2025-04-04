"""Functions to process spatial data and calculate discordance between deconvolution methods."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def process_spatial_data(
    spatial_object_step1,
    spatial_object_step2,
    threshold_umi: int = 200,
    sort_by: str = "AC-like",
):
    """Process spatial data and calculate discordance between deconvolution methods.

    It first loads the counts and deconvolution data for two deconvolution methods
    (step 1 and step 2). It then filters the spots based on the UMI threshold.
    Then, it concatenates the data to create dataframes.
    Finally, the function calculates the discordance between the two deconvolution methods
    and returns several dataframes.

    Parameters
    ----------
    spatial_object_step1 : SpatialObject
        Spatial object with the data for the first deconvolution method.
    spatial_object_step2 : SpatialObject
        Spatial object with the data for the second deconvolution method.
    threshold_umi : int, optional
        Threshold to filter spots based on the UMI count. Default is 200.
    sort_by : str, optional
        Column to sort the discordant_samples DataFrame. Default is "AC-like".

    Returns
    -------
    discordant_samples : pd.DataFrame
        DataFrame with cell types as columns and samples as rows.
        Values are the average discordance between the two deconvolution methods.
    dominant_celltype_freq : pd.DataFrame
        DataFrame with samples as rows and cell types as columns.
        Values are the % of the dominant cell type in each sample.
    concat_abun_step1 : pd.DataFrame
        DataFrame with cell types as columns and samples and spots as rows.
        Values are the cell type abundance for the first deconvolution method.
    concat_abun_step2 : pd.DataFrame
        DataFrame with cell types as columns and samples and spots as rows.
        Values are the cell type abundance for the second deconvolution method.
    """
    # Load counts and deconvolution data for step 1 and step 2
    spatial_object_step1.load_counts_cohort()
    spatial_object_step1.filter_umi(threshold_umi)
    celltype_abund_step1 = spatial_object_step1.load_deconvolution_cohort()

    spatial_object_step2.load_counts_cohort()
    spatial_object_step2.filter_umi(threshold_umi)
    celltype_abund_step2 = spatial_object_step2.load_deconvolution_cohort()

    # Load common cell types
    common_celltypes_step1 = spatial_object_step1.common_celltypes
    common_celltypes_step2 = spatial_object_step2.common_celltypes

    # Concatenate the dataframes
    concat_abun_step2 = pd.concat(
        [
            pd.DataFrame(v, columns=common_celltypes_step2)
            for v in celltype_abund_step2.values()
        ],
        keys=celltype_abund_step2.keys(),
        names=["sample", "spot"],
    )
    concat_abun_step1 = pd.concat(
        [
            pd.DataFrame(v, columns=common_celltypes_step1)
            for v in celltype_abund_step1.values()
        ],
        keys=celltype_abund_step1.keys(),
        names=["sample", "spot"],
    )

    # Calculate discordance between the two deconvolution methods
    discordant_samples = discordant_deconvolution(concat_abun_step1, concat_abun_step2)
    discordant_samples = discordant_samples.sort_values(by=sort_by, ascending=False)

    # Calculate dominant cell type fractions for the second deconvolution method
    celltype_frac = spatial_object_step2.load_deconvolution_cohort(frac=True)

    df_frac = pd.concat(
        [
            pd.DataFrame(v, columns=common_celltypes_step2)
            for v in celltype_frac.values()
        ],
        keys=celltype_frac.keys(),
        names=["sample", "spot"],
    )

    # Calculate dominant cell type frequency
    dominant_celltype_freq = dominant_celltype_frac(
        df_frac, common_celltypes_step2, freq=True
    )

    return (
        discordant_samples,
        dominant_celltype_freq,
        concat_abun_step1,
        concat_abun_step2,
    )


def discordant_deconvolution(df_deconv1, df_deconv2):
    """Find the samples with most discordance for deconvolution methods.

    Parameters
    ----------
    df_deconv1 : pd.DataFrame
        DataFrame with cell types as columns and samples as rows.
        Values should be the cell type abundance or fraction.
    df_deconv2 : pd.DataFrame
        DataFrame with cell types as columns and samples as rows.
        Values should be the cell type abundance or fraction.

    Returns
    -------
    df_discordance : pd.DataFrame
        DataFrame with cell types as columns and samples as rows.
        Values are the average discordance between the two deconvolution methods.
    """
    # First assess the cell types intersection
    common_celltypes = list(set(df_deconv1.columns) & set(df_deconv2.columns))
    if len(common_celltypes) == 0:
        raise ValueError("No common cell types between the two deconvolution methods.")
    else:
        df_deconv1 = df_deconv1[common_celltypes]
        df_deconv2 = df_deconv2[common_celltypes]
    # Calculate the absolute difference between the two deconvolution methods
    df_diff = (df_deconv1 - df_deconv2).abs()
    df_diff.reset_index(inplace=True)
    df_diff.drop(columns="spot", inplace=True)
    df_diff.set_index("sample", inplace=True)
    df_discordance = pd.DataFrame(
        index=list(df_diff.index.unique()), columns=df_diff.columns
    )
    for sample in list(df_diff.index.unique()):
        df_discordance.loc[sample] = df_diff.loc[sample].mean(axis=0)
    return df_discordance


def dominant_celltype_frac(df_frac, common_celltypes, freq=False):
    """Calculate the % of the dominant cell types in each sample.

    Parameters
    ----------
    df_frac : pd.DataFrame
        DataFrame with cell types as columns and samples and spots as rows.
        Values should be the cell type fraction.
    common_celltypes : list
        List of common cell types that in the DataFrame.
    freq : bool, optional
        If True, the output will be the frequency of the dominant cell type in each sample.
        If False, the output will be the relative abundance of each cell type in each spot.

    Returns
    -------
    df_dominant : pd.DataFrame
        DataFrame with samples as rows and cell types as columns.
        Values are the % of the dominant cell type in each sample.
    """
    # calculate dominant cell type in the frac data
    df_frac["dominant_celltype"] = df_frac.loc[:, common_celltypes].idxmax(axis=1)

    if freq:
        # for each sample calculate the frequency of the dominant cell type
        df_frac = (
            df_frac.groupby("sample")["dominant_celltype"]
            .value_counts(normalize=True)
            .unstack()
            .fillna(0)
        )
        # convert to frequency (rounded, no comma)
        df_frac = (df_frac * 100).round(0).astype(int)
        df_frac = df_frac.rename(columns={col: f"{col} %" for col in df_frac.columns})

    return df_frac


def violin_plot_celltype_abundance(df_abun, sample):
    """Create a violin plot of the cell type abundance for each sample."""
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=df_abun.loc[sample], density_norm="width", inner="point")
    plt.title("Violin Plot of Cell Type Abundance")
    plt.ylabel("Abundance")
    plt.xticks(rotation=45)
    plt.show()


def pdf_violin_plot_across_cohort(df_abun, filename):
    """Create a PDF with violin plots of the cell type abundance for each sample."""
    plt.switch_backend("agg")
    from matplotlib.backends.backend_pdf import PdfPages

    # create a pdf file
    pdf = PdfPages(filename)
    # increase margin
    plt.rcParams["figure.autolayout"] = True

    # create a copy
    df_abun_copy = df_abun.copy()
    df_abun_copy.reset_index(inplace=True)
    df_abun_copy.drop(columns="spot", inplace=True)
    df_abun_copy.set_index("sample", inplace=True)

    # for each sample, create a violin plot of the cell type abundance
    for sample in list(df_abun_copy.index.unique()):
        plt.figure(
            # reduce the height of the plot
            figsize=(8, 4)
        )
        sns.violinplot(
            data=df_abun_copy.loc[sample], density_norm="width", inner="point"
        )
        plt.title(f"Violin Plot of Cell Type Abundance for Sample {sample}")
        plt.ylabel("Abundance")
        plt.xticks(rotation=45)
        pdf.savefig()
        plt.close()

    # close the pdf file
    pdf.close()
