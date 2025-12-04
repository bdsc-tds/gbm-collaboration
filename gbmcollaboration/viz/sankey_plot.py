"""Utils function to plot a Sankey plot comparing two annotations."""

# Packages
import pandas as pd
import plotly.graph_objects as go
from typing import Optional

# Function
def plot_sankey_plot(
    df: pd.DataFrame,
    source_col: str,
    target_col: str,
    title: str = "Sankey plot",
    output_html: Optional[str] = None,
    font_size: int = 10,
    pad: int = 15,
    thickness: int = 20
) -> None:
    """
    Plot a Sankey diagram comparing two categorical annotations in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing the annotations to compare.
    source_col : str
        Column name in `df` for the source annotation (e.g., previous label).
    target_col : str
        Column name in `df` for the target annotation (e.g., new label).
    title : str, optional
        Title of the Sankey plot. Default is "Sankey Plot".
    output_html : str, optional
        If provided, saves the plot as an interactive HTML file to this path.
    font_size : int, optional
        Font size used in the plot. Default is 10.
    pad : int, optional
        Padding between Sankey nodes. Default is 15.
    thickness : int, optional
        Thickness of the Sankey nodes. Default is 20.

    Returns
    -------
    None
        Displays the Sankey plot and optionally saves it to an HTML file.
    """
    # Count transitions between the two annotations
    transition_counts = df.groupby([source_col, target_col]).size().reset_index(name="Count")

    # Create labels for the Sankey nodes
    labels = list(pd.unique(transition_counts[[source_col, target_col]].values.ravel()))
    label_index = {label: i for i, label in enumerate(labels)}

    # Map annotations to node indices
    sources = transition_counts[source_col].map(label_index)
    targets = transition_counts[target_col].map(label_index)
    values = transition_counts["Count"]

    # Create Sankey diagram
    fig = go.Figure(go.Sankey(
        node=dict(
            pad=pad,
            thickness=thickness,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values
        )
    ))

    fig.update_layout(title_text=title, font_size=font_size)

    if output_html:
        fig.write_html(output_html)
    fig.show()
