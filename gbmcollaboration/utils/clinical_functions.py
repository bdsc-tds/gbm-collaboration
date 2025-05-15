
"""Useful functions for clinical data analysis."""

import itertools
from typing import Dict, List, Optional

import lifelines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.preprocessing import PowerTransformer
from statsmodels.stats.multitest import fdrcorrection


# Useful functions
def plot_correlations(df: pd.DataFrame, values: str, var_list: List[str]):  # noqa: PLR0912, PLR0915
    """Plot a matrix of correlations or p-values.

    Parameters
    ----------
    df: pd.DataFrame
    values: str
        What should be plotted. Either 'FDR' or 'correlations'
    """
    missing_vars = set(var_list) - set(df.columns)
    if missing_vars:
        raise ValueError(
            f"The following variables are missing from the DataFrame: {missing_vars}"
        )
    df_subset = df[var_list].copy()
    nvars = len(var_list)

    # Initialize matrices
    r_arr = df_subset.corr(method="pearson").values
    p_arr = np.zeros((nvars, nvars))

    # Compute p-values
    for i, j in itertools.combinations(range(nvars), 2):
        tmp = df_subset.iloc[:, [i, j]].dropna()
        if tmp.empty:
            p = np.nan
        else:
            _, p = pearsonr(tmp.iloc[:, 0], tmp.iloc[:, 1])
        p_arr[i, j] = p
        p_arr[j, i] = p  # Symmetric matrix

        # Create mask
        mask = np.triu(np.ones_like(r_arr, dtype=bool))

        # Extract the traingular inferior matrix of df
        p_vect = p_arr[np.tril_indices(nvars, k=-1)]
        p_vect = p_vect[~np.isnan(p_vect)]

    # Apply FDR correction to the lower triangle
    _, p_corr_vect = fdrcorrection(p_arr[np.tril_indices(nvars, k=-1)])

    # Initialize arrays
    sign_arr = [["" for _ in var_list] for _ in var_list]
    sign_arr_dis = np.zeros((nvars, nvars))
    p_corr_arr = np.full((nvars, nvars), np.nan)

    # Assign corrected p-values and significance
    ind = 0
    for i in range(nvars):
        for j in range(nvars):
            if j < i:
                p_corr = p_corr_vect[ind]
                p_corr_arr[i, j] = p_corr
                if p_corr < 0.001:
                    sign = "***"
                    sign_arr_dis[i, j] = 3
                elif p_corr < 0.01:
                    sign = "**"
                    sign_arr_dis[i, j] = 2
                elif p_corr < 0.05:
                    sign = "*"
                    sign_arr_dis[i, j] = 1
                else:
                    sign = "ns"
                    sign_arr_dis[i, j] = 0
                sign_arr[i][j] = sign
                ind += 1
            else:
                sign_arr[i][j] = ""

    # Plot
    f, ax = plt.subplots(figsize=(nvars, nvars))

    if values == "FDR":
        sns.heatmap(
            sign_arr_dis,
            mask=mask,
            square=True,
            linewidths=0.5,
            cmap="coolwarm",
            cbar_kws={"shrink": 0.4, "ticks": [-1, -0.5, 0, 0.5, 1]},
            vmin=0,
            vmax=3,
            fmt="",
            annot=sign_arr,
            annot_kws={"size": 12},
            cbar=False,
        )

        # add the column names as labels
        ax.set_yticklabels(df_subset.columns, rotation=0)
        ax.set_xticklabels(df_subset.columns, rotation=90)

        sns.set_style({"xtick.bottom": True}, {"ytick.left": True})

        plt.title("False discovery rate (FDR)")

    elif values == "correlations":
        sns.heatmap(
            r_arr,
            mask=mask,
            square=True,
            linewidths=0.5,
            cmap="coolwarm",
            cbar_kws={"shrink": 0.4, "ticks": [-1, -0.5, 0, 0.5, 1]},
            vmin=-1,
            vmax=1,
            fmt="",
            annot=np.round(r_arr, 1),
            annot_kws={"size": 12},
        )

        # add the column names as labels
        ax.set_yticklabels(df_subset.columns, rotation=0)
        ax.set_xticklabels(df_subset.columns, rotation=90)

        sns.set_style({"xtick.bottom": True}, {"ytick.left": True})

        plt.title("Pearson's correlation coefficient (r)")

    else:
        raise ValueError(
            "Unknown 'values' argument. Supported values are 'FDR' and 'correlations'."
        )


def univariate(
    df: pd.DataFrame, endpoint: str, var_names_all: Dict[str, str]
) -> pd.DataFrame:
    """Perform univariate analysis.

    Parameters
    ----------
    df: DataFrame
        Clinical data.
    endpoint: str
        Name of the column containing the endpoint. Ex: 'label' or 'OS'.
    var_names_all: Dict[str, str]
        Names of the clinical variables and associated type ('categorical' or 'continuous')

    Returns
    -------
    df_all : DataFrame
        Univariate analysis in the same format as a cph.summary dataframe.
    """
    event_col = f"{endpoint}.event"
    time_col = f"{endpoint}.time"

    scaler = PowerTransformer()
    is_first_var = True

    # todo Convert the time-to-event label into sklearn format
    df.loc[:, time_col] = df[endpoint].abs()
    df.loc[:, event_col] = df[endpoint] > 0

    # Compute univariate analysis for each variable
    for var, var_type in var_names_all.items():
        data_train = df[[time_col, event_col, var]]

        # Remove patients for which the variable value is unknown
        data_train = data_train.dropna()

        # Normalize values of continuous variables
        if var_type == "continuous":
            data_train[var] = scaler.fit_transform(
                np.asarray(data_train[var]).reshape(-1, 1)
            )

        cph = lifelines.CoxPHFitter()
        cph.fit(data_train, time_col, event_col=event_col, show_progress=False)

        cph_summary = cph.summary

        cph_summary["c_index"] = lifelines.utils.concordance_index(
            data_train[time_col], data_train[var], data_train[event_col]
        )

        # For the first variable, initiate the dataframe. Then, append information for
        # the other variables.
        if is_first_var:
            df_all = cph_summary
            is_first_var = False
        else:
            df_all = pd.concat([df_all, cph_summary])

    # Add false discovery rate
    _, df_all["fdr"] = fdrcorrection(df_all["p"])

    return df_all


def plot_forest(
    df_all: pd.DataFrame,
    distances: Optional[List[float]] = None,
    figsize: Optional[tuple[float, float]] = (10, 10),
    outcome: Optional[str] = None,
    significance: str = "FDR",
):
    """Generate a forest plot.

    Parameters
    ----------
    df_all : DataFrame
        Univariate analysis in the same format as a cph.summary dataframe.
    distances: Optional[List[float]]
        Distance between the name of the genes and the bars, between the bars
        and the HRs, and between the HRs and the p-values/FDR.
        Default is [0.2, 0.2, 0.2]
    figsize: Optional[Tuple[int]]
        Figure size
    outcome: Optional[str]
        'OS', 'DSS', 'PFS', or 'PFI'
    significance: str
        Significance measure to be used (p-value or FDR)
    """

    def get_limits(var, df):
        x0 = df[df["covariate"] == var]["exp(coef) lower 95%"].values[0]
        x1 = df[df["covariate"] == var]["exp(coef) upper 95%"].values[0]
        x = df[df["covariate"] == var]["exp(coef)"].values[0]
        return x0, x1, x

    def get_color(significance_val):
        if significance_val < 0.05:
            return "black"
        return "gray"

    def plot_significance_text(significance_val, ind, dist3):
        if significance_val >= 0.001:
            plt.text(dist3, ind, f"{np.round(significance_val, 3)}")
        else:
            plt.text(dist3, ind, "<0.001")

    space = 10
    var_names_all = df_all["covariate"]
    plt.figure(figsize=figsize)
    ax = plt.axes(frameon=False)
    ax.axes.get_yaxis().set_visible(False)  # type: ignore[union-attr]

    mini = df_all["exp(coef) lower 95%"].min()
    maxi = df_all["exp(coef) upper 95%"].max()

    if distances is None:
        distances = [0.2, 0.2, 0.2]

    dist1, dist2, dist3 = (
        mini - (maxi - mini) * distances[0],
        maxi + (maxi - mini) * distances[1],
        maxi + (maxi - mini) * distances[2],
    )

    plt.text(dist2, space, "HR (95% CI)", weight="bold")
    plt.text(dist3, space, significance, weight="bold")

    for ind, var in enumerate(var_names_all):
        x0, x1, x = get_limits(var, df_all)
        significance_val = df_all[df_all["covariate"] == var][
            significance.lower()
        ].values[0]
        c = get_color(significance_val)

        plt.plot([x0, x1], [ind * -space, ind * -space], color=c)
        plt.scatter(x, ind * -space, marker="o", color=c)
        plt.text(dist1, ind * -space, var)
        plt.text(
            dist2,
            ind * -space,
            f"{np.round(x, 2)} ({np.round(x0, 2)}-{np.round(x1, 2)})",
        )

        plot_significance_text(significance_val, ind * -space, dist3)

    if outcome:
        plt.xlabel(f"{outcome} Hazard Ratio")

    plt.axvline(x=1, color="black", linestyle="--")

    return plt
