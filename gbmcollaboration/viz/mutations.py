"""Visualization utilities for mutation analysis."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import ListedColormap, to_rgba
from typing import Tuple


def create_oncoplot(patient_gene_matrix: pd.DataFrame, figsize: Tuple[int, int] = (20, 10)) -> None:
    """
    Create an oncoplot-style heatmap showing alterations per patient-gene pair.
    
    Parameters
    ----------
    patient_gene_matrix : pd.DataFrame
        Patient X gene matrix with columns like {gene}_mut, {gene}_amp, {gene}_lof
        Rows should be samples/patients, columns should be alteration types
    figsize : tuple
        Figure size (width, height)
        
    Returns
    -------
    None
        Displays the plot
    """
    # Extract unique genes from column names
    genes = set()
    for col in patient_gene_matrix.columns:
        if '_' in col:
            gene = col.rsplit('_', 1)[0]
            genes.add(gene)
        else:
            # Summary column (just gene name)
            genes.add(col)
    
    genes = sorted(genes)
    samples = patient_gene_matrix.index.tolist()
    
    # Create matrix: 0 = no alteration, 1 = mut, 2 = lof, 3 = amp
    # Priority: mut > lof > amp (if multiple, show mut)
    plot_matrix = []
    gene_labels = []
    
    for gene in genes:
        row = []
        mut_col = f"{gene}_mut" if f"{gene}_mut" in patient_gene_matrix.columns else None
        lof_col = f"{gene}_lof" if f"{gene}_lof" in patient_gene_matrix.columns else None
        amp_col = f"{gene}_amp" if f"{gene}_amp" in patient_gene_matrix.columns else None
        
        for sample in samples:
            value = 0  # No alteration (transparent)
            if mut_col and patient_gene_matrix.loc[sample, mut_col]:
                value = 1  # Mutation (green)
            elif lof_col and patient_gene_matrix.loc[sample, lof_col]:
                value = 2  # Loss-of-function (red)
            elif amp_col and patient_gene_matrix.loc[sample, amp_col]:
                value = 3  # Amplification (blue)
            
            row.append(value)
        
        # Only add row if gene has at least one alteration
        if any(v > 0 for v in row):
            plot_matrix.append(row)
            gene_labels.append(gene)
    
    if not plot_matrix:
        print("No alterations found to plot.")
        return
    
    plot_matrix = np.array(plot_matrix)
    
    # Create custom colormap: transparent, green, red, blue
    colors = [
        to_rgba('white', alpha=0.0),  # 0: transparent (no alteration)
        to_rgba('green', alpha=0.8),  # 1: mutation
        to_rgba('red', alpha=0.8),    # 2: loss-of-function
        to_rgba('blue', alpha=0.8),   # 3: amplification
    ]
    cmap = ListedColormap(colors)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    sns.heatmap(
        plot_matrix,
        cmap=cmap,
        cbar=False,
        xticklabels=samples,
        yticklabels=gene_labels,
        linewidths=0.5,
        linecolor='lightgray',
        ax=ax
    )
    
    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=90, ha='right', fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=10)
    
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Driver Genes', fontsize=12)
    ax.set_title('Driver Alterations: Mutations (Green), Loss-of-Function (Red), Amplifications (Blue)', 
                 fontsize=14, pad=20)
    
    # Create legend
    legend_elements = [
        mpatches.Patch(facecolor='green', alpha=0.8, label='Mutation'),
        mpatches.Patch(facecolor='red', alpha=0.8, label='Loss-of-Function'),
        mpatches.Patch(facecolor='blue', alpha=0.8, label='Amplification'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1))
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nPlotted {len(gene_labels)} genes with alterations across {len(samples)} samples")

