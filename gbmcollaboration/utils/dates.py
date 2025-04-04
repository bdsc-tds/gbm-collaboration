"""Utils functions for dates processing."""

from typing import Any, Optional

import numpy as np
import pandas as pd
from loguru import logger

from gbmcollaboration.definitions import AVG_NUM_DAYS_PER_MONTH, AVG_NUM_DAYS_PER_YEAR


def convert_to_datetime(
    df: pd.DataFrame,
    col: str,
    dayfirst: bool,
    format: Optional[str] = None,
    verbose: bool = False,
) -> None:
    """
    Convert a column of a dataframe to datetime in place.

    Parameters
    ----------
    df: pd.DataFrame
        dataframe with the column to convert.
    col: str
        name of the column to convert.
    dayfirst: bool
        if True, the first value of the date is the day.
    format: str
        format of the date.
    verbose: bool
        if True, print warning for values that are not time and were not NaNs.

    Raises
    ------
    ValueError
        if the column cannot be converted to datetime.
    """
    try:
        if format is not None:
            new_col_values = pd.to_datetime(df[col], format=format, exact=True)
        else:
            new_col_values = pd.to_datetime(df[col], dayfirst=dayfirst)
    except Exception as e:
        error_msg = "Could not convert column '%s' to datetime. See:\n%s\n%s" % (
            col,
            str(df.loc[~df[col].isna(), col].head()),
            str(e),
        )
        raise ValueError(error_msg) from e

    if verbose:
        # report warning for values that are not time and were not NaNs
        nat_index = new_col_values.isnull()
        nan_index = df[col].notna()
        wrong_index = nat_index & nan_index
        if wrong_index.sum() > 0:
            logger.warning("Found values not convertible to time in column '%s':" % col)
            logger.warning(df.loc[wrong_index, col])

    df[col] = new_col_values


def time_to_event(
    df: pd.DataFrame,
    col_start: str,
    col_end: str,
    col_censored: str,
    col_end_censored: str,
    unit: str = "years",
) -> np.ndarray[Any, np.dtype[Any]]:
    """
    Compute the time difference between two dates.

    Parameters
    ----------
    df: pd.DataFrame
        dataframe with the dates.
    col_start: str
        name of the first date column.
    col_end: str
        name of the second date column.
    col_censored: str
        name of the column indicating if the event is censored. Should contain boolean.
    col_end_censored: str
        name of the column containing the censored date.
    unit: str
        unit of the time difference output.

    Returns
    -------
    pd.Series
        the time difference between the two dates. Censored dates are encoded as negative values.
    """
    # Check the time columns data types
    if df[col_start].dtype != "datetime64[ns]":
        raise ValueError("Column %s is not a datetime column." % col_start)
    if df[col_end].dtype != "datetime64[ns]":
        raise ValueError("Column %s is not a datetime column." % col_end)
    if df[col_end_censored].dtype != "datetime64[ns]":
        raise ValueError("Column %s is not a datetime column." % col_end_censored)
    if df[col_censored].dtype != "bool":
        raise ValueError("Column %s is not a boolean column." % col_censored)

    delta = np.where(
        df[col_censored],
        (-1)
        * (df[col_end_censored] - df[col_start]).dt.days,  # if censored (e.g. alive)
        (df[col_end] - df[col_start]).dt.days,  # if not censored (e.g. dead)
    )
    if unit == "months":
        delta = delta / AVG_NUM_DAYS_PER_MONTH
    elif unit == "years":
        delta = delta / AVG_NUM_DAYS_PER_YEAR
    elif unit != "days":
        raise ValueError("Unit %s not supported." % unit)

    return delta
