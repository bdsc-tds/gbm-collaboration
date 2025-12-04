"""Script to produce the demographic tables across centers."""

import os

import numpy as np
import pandas as pd

import os
import numpy as np
import pandas as pd


def _coerce_to_boolean(series: pd.Series) -> pd.Series:
    """Safely coerce common boolean representations to pandas BooleanDtype (nullable)."""
    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")  # nullable bool

    # Map common truthy/falsey forms; leave others as NA
    mapping = {
        "true": True, "false": False,
        "yes": True, "no": False,
        "y": True, "n": False,
        "1": True, "0": False,
        1: True, 0: False,
        True: True, False: False,
    }
    s = series.copy()

    # Handle strings safely
    is_str = s.map(lambda x: isinstance(x, str))
    s.loc[is_str] = s.loc[is_str].str.strip().str.lower()

    # Map values and keep NaNs where mapping not found
    s = s.map(mapping).astype("boolean")
    return s


def create_demographic_table(df: pd.DataFrame, path_to_file: str, verbose: bool = True):
    """Create demographic table.

    Parameters

    df : pd.DataFrame
        Dataframe containing all features at baseline.
    path_to_file : str
        Path to the folder containing the template demo to fill.
    verbose : bool, optional
        Whether to print information about the process, by default True

    Returns
    -------
    pd.DataFrame
        Demographic table.
    """
    # load template
    template_demo = pd.read_excel(
        os.path.join(path_to_file, "Mosaic Cohort Description Table.xlsx"),
        sheet_name="Demographic Table",
        header=0,
        usecols=["ID", "Category", "Subcategory", "Variable", "Type", "Unit"],
    )

    # initiate demographic table
    demographic_table = template_demo.copy()
    demographic_table["Observed"] = np.nan
    demographic_table["Missing"] = np.nan
    demographic_table["Percentage"] = np.nan
    demographic_table["Min"] = np.nan
    demographic_table["Max"] = np.nan
    demographic_table["Mean"] = np.nan
    demographic_table["SD"] = np.nan
    # your new columns
    demographic_table["Count"] = np.nan
    demographic_table["Categories"] = np.nan

    # Pre-cast df columns based on template "Type"
    for _, row in template_demo.iterrows():
        var = row["Variable"]
        typ = str(row["Type"]).strip().lower()
        if var not in df.columns:
            continue
        if typ == "boolean":
            df[var] = _coerce_to_boolean(df[var])
        elif typ == "categorical":
            # ensure object for categorical summaries
            df[var] = df[var].astype("object")

    # Fill table
    for col in template_demo.Variable:
        if verbose:
            print(f"Looking for variable {col}")
        if col in df.columns:
            expected_type = str(
                template_demo.loc[template_demo["Variable"] == col, "Type"].values[0]
            ).strip().lower()

            # calculate missing and observed values
            observed = df[col].notna().sum()
            missing = df[col].isna().sum()
            percentage = (missing / (observed + missing) * 100) if (observed + missing) else np.nan

            demographic_table.loc[demographic_table["Variable"] == col, "Observed"] = np.round(observed)
            demographic_table.loc[demographic_table["Variable"] == col, "Missing"] = np.round(missing)
            demographic_table.loc[demographic_table["Variable"] == col, "Percentage"] = np.round(percentage, 2)

            if expected_type == "boolean":
                # Leave min/max/mean/sd empty (NaN)
                # Count = number of "Yes" -> interpreted as True after coercion
                yes_count = (df[col] == True).sum(skipna=True)  # noqa: E712
                demographic_table.loc[demographic_table["Variable"] == col, "Count"] = np.round(yes_count)
                # Categories left empty unless you also want breakdown; not requested
                if verbose:
                    print(f"Information about variable {col} as boolean added (Count=True).")

            elif expected_type == "categorical":
                # Categorical summary:
                # Categories = dict of category -> count (exclude NaN)
                counts = df[col].value_counts(dropna=True).to_dict()
                demographic_table.loc[demographic_table["Variable"] == col, "Min"] = np.nan
                demographic_table.loc[demographic_table["Variable"] == col, "Max"] = np.nan
                demographic_table.loc[demographic_table["Variable"] == col, "Mean"] = np.nan
                demographic_table.loc[demographic_table["Variable"] == col, "SD"] = np.nan
                demographic_table.loc[demographic_table["Variable"] == col, "Categories"] = str(counts)
                # Count left as NaN (rule only asked for boolean)
                if verbose:
                    print(f"Information about variable {col} as categorical added (Categories dict).")

            else:
                # numeric (float) summary as before
                if pd.api.types.is_numeric_dtype(df[col]) and not pd.api.types.is_bool_dtype(df[col]):
                    minimum = df[col].min()
                    maximum = df[col].max()
                    mean = df[col].mean()
                    std = df[col].std()
                    demographic_table.loc[demographic_table["Variable"] == col, "Min"] = np.round(minimum, 3)
                    demographic_table.loc[demographic_table["Variable"] == col, "Max"] = np.round(maximum, 3)
                    demographic_table.loc[demographic_table["Variable"] == col, "Mean"] = np.round(mean, 3)
                    demographic_table.loc[demographic_table["Variable"] == col, "SD"] = np.round(std, 3)
                    if verbose:
                        print(f"Information about variable {col} as numeric added.")
                else:
                    if verbose:
                        print(f"Variable {col} has unaccepted or mixed type for numeric summary.")
        elif verbose:
            print(f"Variable {col} not found in this center.")

    # Final formatting
    # Convert Observed, Missing, Count to integers (drop .0)
    for col_int in ["Observed", "Missing", "Count"]:
        if col_int in demographic_table.columns:
            demographic_table[col_int] = demographic_table[col_int].apply(
                lambda x: int(x) if pd.notnull(x) else np.nan
            )

    # Round percentage, Min, Max, Mean, SD to 2 decimals
    for col_float in ["Percentage", "Min", "Max", "Mean", "SD"]:
        if col_float in demographic_table.columns:
            demographic_table[col_float] = demographic_table[col_float].apply(
                lambda x: round(x, 2) if pd.notnull(x) else np.nan
            )

    return demographic_table
