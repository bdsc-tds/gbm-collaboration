"""Visualization utilities for mutation analysis."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
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
    
    samples = patient_gene_matrix.index.tolist()
    
    # Calculate proportions for ALL genes (including those with no alterations)
    gene_proportions = {}
    
    for gene in genes:
        mut_col = f"{gene}_mut" if f"{gene}_mut" in patient_gene_matrix.columns else None
        lof_col = f"{gene}_lof" if f"{gene}_lof" in patient_gene_matrix.columns else None
        amp_col = f"{gene}_amp" if f"{gene}_amp" in patient_gene_matrix.columns else None
        
        # Count samples with any alteration for this gene
        samples_with_alteration = 0
        for sample in samples:
            has_alteration = False
            if mut_col and patient_gene_matrix.loc[sample, mut_col]:
                has_alteration = True
            elif lof_col and patient_gene_matrix.loc[sample, lof_col]:
                has_alteration = True
            elif amp_col and patient_gene_matrix.loc[sample, amp_col]:
                has_alteration = True
            
            if has_alteration:
                samples_with_alteration += 1
        
        proportion = samples_with_alteration / len(samples) if samples else 0
        gene_proportions[gene] = proportion
    
    if not genes:
        print("No genes found to plot.")
        return
    
    # Sort genes by proportion (descending - most frequent first, genes with 0% at the end)
    genes_sorted = sorted(genes, key=lambda g: gene_proportions.get(g, 0), reverse=True)
    
    # Create matrix: store multiple alteration types per cell
    # Each cell is a list of alteration types: [] = no alteration, [1] = mut, [2] = lof, [3] = amp
    # Or combinations like [1, 3] = mut + amp
    plot_matrix = []
    gene_labels = []
    
    for gene in genes_sorted:
        row = []
        mut_col = f"{gene}_mut" if f"{gene}_mut" in patient_gene_matrix.columns else None
        lof_col = f"{gene}_lof" if f"{gene}_lof" in patient_gene_matrix.columns else None
        amp_col = f"{gene}_amp" if f"{gene}_amp" in patient_gene_matrix.columns else None
        
        for sample in samples:
            alterations = []
            if mut_col and patient_gene_matrix.loc[sample, mut_col]:
                alterations.append(1)  # Mutation (green)
            if lof_col and patient_gene_matrix.loc[sample, lof_col]:
                alterations.append(2)  # Loss-of-function (red)
            if amp_col and patient_gene_matrix.loc[sample, amp_col]:
                alterations.append(3)  # Amplification (blue)
            
            row.append(alterations)
        
        plot_matrix.append(row)
        # Add percentage to gene label (even if 0%)
        proportion = gene_proportions.get(gene, 0)
        percentage = proportion * 100
        gene_labels.append(f"{gene} ({percentage:.1f}%)")
    
    n_genes = len(plot_matrix)
    n_samples = len(samples) if plot_matrix else 0
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Color mapping
    color_map = {
        1: 'green',   # Mutation
        2: 'red',     # Loss-of-function
        3: 'blue',    # Amplification
    }
    
    # Draw cells manually to support split cells
    for i, row in enumerate(plot_matrix):
        for j, alterations in enumerate(row):
            if not alterations:
                # No alteration - white/transparent
                rect = Rectangle((j, n_genes - 1 - i), 1, 1, 
                               facecolor='white', edgecolor='lightgray', linewidth=0.5)
                ax.add_patch(rect)
            elif len(alterations) == 1:
                # Single alteration - fill entire cell
                color = color_map[alterations[0]]
                rect = Rectangle((j, n_genes - 1 - i), 1, 1,
                               facecolor=color, alpha=0.8, edgecolor='lightgray', linewidth=0.5)
                ax.add_patch(rect)
            else:
                # Multiple alterations - split cell vertically (top/bottom)
                n_alt = len(alterations)
                if n_alt == 2:
                    # Split vertically: top half and bottom half
                    # Top half: first alteration
                    color1 = color_map[alterations[0]]
                    rect1 = Rectangle((j, n_genes - 1 - i + 0.5), 1, 0.5,
                                    facecolor=color1, alpha=0.8, edgecolor='lightgray', linewidth=0.5)
                    ax.add_patch(rect1)
                    # Bottom half: second alteration
                    color2 = color_map[alterations[1]]
                    rect2 = Rectangle((j, n_genes - 1 - i), 1, 0.5,
                                    facecolor=color2, alpha=0.8, edgecolor='lightgray', linewidth=0.5)
                    ax.add_patch(rect2)
                elif n_alt == 3:
                    # Split vertically into three parts: top, middle, bottom
                    color1 = color_map[alterations[0]]
                    color2 = color_map[alterations[1]]
                    color3 = color_map[alterations[2]]
                    # Top third
                    rect1 = Rectangle((j, n_genes - 1 - i + 2/3), 1, 1/3,
                                    facecolor=color1, alpha=0.8, edgecolor='lightgray', linewidth=0.5)
                    ax.add_patch(rect1)
                    # Middle third
                    rect2 = Rectangle((j, n_genes - 1 - i + 1/3), 1, 1/3,
                                    facecolor=color2, alpha=0.8, edgecolor='lightgray', linewidth=0.5)
                    ax.add_patch(rect2)
                    # Bottom third
                    rect3 = Rectangle((j, n_genes - 1 - i), 1, 1/3,
                                    facecolor=color3, alpha=0.8, edgecolor='lightgray', linewidth=0.5)
                    ax.add_patch(rect3)
    
    # Set axis limits and labels
    ax.set_xlim(0, n_samples)
    ax.set_ylim(0, n_genes)
    
    # Set y-axis labels (gene names only)
    gene_names_only = [g.split(' (')[0] for g in gene_labels]
    ax.set_yticks([i + 0.5 for i in range(n_genes)])
    ax.set_yticklabels(reversed(gene_names_only), fontsize=10)
    
    # Remove x-axis ticks
    ax.set_xticks([])
    
    # Add percentage annotations on the right side
    for i, (gene, label) in enumerate(zip(genes_sorted, gene_labels)):
        proportion = gene_proportions.get(gene, 0)
        percentage = proportion * 100
        # Position annotation at the right edge of the plot
        ax.text(
            n_samples + 0.5,  # Right edge + small offset
            n_genes - 1 - i + 0.5,  # Center of the row
            f'{percentage:.1f}%',
            ha='left',
            va='center',
            fontsize=9
        )
    
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Driver Genes', fontsize=12)
    ax.set_title('Driver Alterations: Mutations (Green), Loss-of-Function (Red), Amplifications (Blue)', 
                 fontsize=14, pad=20)
    
    # Create legend - move further out
    legend_elements = [
        mpatches.Patch(facecolor='green', alpha=0.8, label='Mutation'),
        mpatches.Patch(facecolor='red', alpha=0.8, label='Loss-of-Function'),
        mpatches.Patch(facecolor='blue', alpha=0.8, label='Amplification'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nPlotted {len(gene_labels)} genes with alterations across {len(samples)} samples")

