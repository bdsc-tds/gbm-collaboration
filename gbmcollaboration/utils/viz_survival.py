"""Distribution plots."""

from collections.abc import Iterable
from typing import Any, Optional
from loguru import logger

import matplotlib
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import (
    KaplanMeierFitter,
    logrank_test,
    multivariate_logrank_test,
)
from matplotlib import pyplot as plt
from pandas.api.types import is_integer_dtype

AVG_NUM_DAYS_PER_MONTH = 30.437
DEFAULT_FIG_SIZE = (3, 4)


def check_column(df: pd.DataFrame, col_name: str) -> None:
    """Check that the column is in the dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe to check
    col_name : str
        column name to check

    Raises
    ------
    ValueError
        If the column is not in the dataframe.
    """
    if col_name not in df.columns:
        logger.error(
            f"Column {col_name} not found in dataframe, available columns: {df.columns}"
        )
        raise ValueError(f"Column {col_name} not found in dataframe")


def check_hue_order(df: pd.DataFrame, hue: str, hue_order: Iterable):
    """Check that hue_order contains all and only the values in df[hue]."""
    if not set(df[hue].unique()).issubset(set(hue_order)):
        logger.error(
            f"Values in {hue} not found in hue_order:"
            f" {set(df[hue].unique()) - set(hue_order)}"
        )
        raise ValueError(
            f"Values in {hue} not found in hue_order:"
            f" {set(df[hue].unique()) - set(hue_order)}"
        )
    if not set(hue_order).issubset(set(df[hue].unique())):
        logger.error(
            f"Values in hue_order not found in {hue}:"
            f" {set(hue_order) - set(df[hue].unique())}"
        )
        raise ValueError(
            f"Values in hue_order not found in {hue}:"
            f" {set(hue_order) - set(df[hue].unique())}"
        )


def _checks(data: pd.DataFrame, x: str, hue: Optional[str] = None, hue_order=None):
    """Check a dataframe to be used with survival plots."""
    check_column(data, x)
    # Check that values in data[x] are integers
    if not is_integer_dtype(data[x]):
        raise ValueError(f"Column {x} must be of integer type. Got {data[x].dtype}")
    if hue is not None:
        check_column(data, hue)
    if hue_order is not None and hue is None:
        raise ValueError("Cannot use hue_order without hue")
    if hue is not None and hue_order is not None:
        check_hue_order(data, hue, hue_order)


def _compute_hazard_ratio(df: pd.DataFrame, x: str, hue: str) -> tuple[Any, list[Any]]:
    """Compute the hazard ratio for the data.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to use.
    x : str
        Column name for the duration to plot.
    hue : str
        Column name to use for the hue.

    Returns
    -------
    Any
        Hazard ratio.
    list[Any]
        Confidence interval for the hazard ratio.

    Raises
    ------
    ValueError
        If hue is None.
    ValueError
        If there are less than 2 values in hue.
    """
    if hue is None:
        raise ValueError("Cannot compute hazard ratio without hue")
    hue_values = df[hue].unique()
    if len(hue_values) < 2:
        raise ValueError(
            f"Cannot compute hazard ratio with less than 2 values in {hue}"
        )
    if len(hue_values) > 2:
        raise ValueError("Cannot compute hazard ratio with more than 2 values in hue")

    cph = CoxPHFitter()
    # Convert hue to a integer variable
    df_tmp = df.copy()[[x, "event", hue]]
    df_tmp[hue] = df_tmp[hue].apply(lambda x: 0 if x == hue_values[0] else 1)
    df_tmp[hue] = df_tmp[hue].astype(int)
    cph.fit(df_tmp, duration_col=x, event_col="event")
    hr = cph.summary["exp(coef)"][hue]
    hr_ci = [
        cph.summary["exp(coef) lower 95%"][hue],
        cph.summary["exp(coef) upper 95%"][hue],
    ]
    return hr, hr_ci


def _compute_logrank_pvalue(df: pd.DataFrame, x: str, hue: str) -> float:
    """Compute the logrank p-value for the data.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to use.
    x : str
        Column name for the duration to plot.
    hue : str
        Column name to use for the hue.

    Returns
    -------
    float
        Logrank p-value.

    Raises
    ------
    ValueError
        If hue is None.
    ValueError
        If there are less than 2 values in hue.
    """
    if hue is None:
        raise ValueError("Cannot compute logrank p-value without hue")
    hue_values = df[hue].unique()
    if len(hue_values) < 2:
        raise ValueError(
            f"Cannot compute logrank p-value with less than 2 values in {hue}"
        )

    if len(hue_values) == 2:
        mask1 = df[hue] == hue_values[0]
        mask2 = df[hue] == hue_values[1]
        results = logrank_test(
            df[mask1][x], df[mask2][x], df[mask1]["event"], df[mask2]["event"]
        )
    else:
        results = multivariate_logrank_test(df[x], df[hue], df["event"])

    return results.p_value


def kmplot(
    data: pd.DataFrame,
    x: str,
    hue: Optional[str] = None,
    hue_order: Optional[Iterable] = None,
    convert_days_to_months: bool = True,
    add_hazard_ratio: bool = False,
    add_logrank_pvalue: bool = False,
    at_risk_counts: bool = False,
    ax: Optional[matplotlib.axes.Axes] = None,
    title: str = None
) -> Optional[matplotlib.figure.Figure]:
    """Plot a Kaplan Meier curve of the data.

    It is assumed that censorship is encoded with a negative value in data[x].
    It is also assumed that the data is given in days.

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to plot.
    x : str
        Column name for the duration to plot.
        data[x] must represent durations and censorship encoded in one value:
        positive values are (known) times to event of interest (e.g. death)
        negative values are censored data with absolute value being the time
        of the last event.
    hue : Optional[str], optional
        Column name to use for the hue, by default None.
        If given multiple curves will be plotted corresponding to
         each possible value of hue in data.
    hue_order : Optional[Iterable], optional
        Order for the hue levels, by default None
    convert_days_to_months: bool
        if True: assumes the inputs are given in days and convert them into months.
    add_hazard_ratio: bool
        if True, add the hazard ratio to the plot.
        Only supported if hue is used and has exactly two possible values.
    add_logrank_pvalue: bool
        if True, add the logrank p-value to the plot.
    ax : Optional[matplotlib.axes.Axes], optional
        Axes to plot on, by default None. Not figure is returned if ax is not None.
        In this case the figure will be added inplace to the axis.

    Returns
    -------
    Optional[matplotlib.figure.Figure]
        Figure with the plot. Return None of ax is not None.
        In this case the figure will be added inplace to the axis.
    """
    df = data.copy()
    df = df[df[x].notna()]  # Remove unknown durations
    df[x] = df[x].astype(int)  # Ensure that durations are integers
    _checks(df, x, hue, hue_order)

    # Prepare data for plot
    unit = "days"
    if convert_days_to_months:
        df[x] /= AVG_NUM_DAYS_PER_MONTH
        unit = "months"
    df["event"] = df[x] > 0  # True iff event was observed
    df[x] = df[x].abs()  # Remove the sign of the duration

    if ax is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=DEFAULT_FIG_SIZE)
    else:
        fig = None

    kmf_list = []

    kmf = KaplanMeierFitter()
    if hue is None:
        kmf.fit(df[x])
        kmf.plot(ax=ax)
        kmf_list.append(kmf)

        if title is None:
            title = x
    else:
        if hue_order is None:
            hue_order = df[hue].unique()
        for hue_value in hue_order:
            mask = df[hue] == hue_value
            n_total = mask.sum()
            n_events = df.loc[mask, "event"].sum()
            n_censored = n_total - n_events

            kmf.fit(df[mask][x], df[mask]["event"])
            kmf.plot(ax=ax, label=f"{hue_value} (n={n_total}, censored={int(n_censored)})")
            kmf_list.append(kmf)

        ax.legend(title=hue, loc="best")
        if title is None:
            title = hue
        if add_logrank_pvalue:
            x_pos = 0.01 * ax.get_xlim()[1]
            y_pos = 0.01
            if (p_value := _compute_logrank_pvalue(df, x, hue)) < 0.001:
                p_val_str = "p<0.001"
            else:
                p_val_str = f"p={p_value:.3f}"
            ax.text(x_pos, y_pos, p_val_str, fontsize=10)
        if add_hazard_ratio:
            hazard_ratio, ci = _compute_hazard_ratio(df, x, hue)
            x_pos = 0.01 * ax.get_xlim()[1]
            y_pos = 0.05
            text = f"HR={hazard_ratio:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]"
            ax.text(x_pos, y_pos, text, fontsize=10)

    # Add at-risk counts
    if at_risk_counts:
        add_at_risk_counts(*kmf_list, ax=ax)

    ax.set_title(title)
    ax.set_xlabel(f"Time ({unit})")
    ax.set_ylabel("Survival probability")

    return fig
