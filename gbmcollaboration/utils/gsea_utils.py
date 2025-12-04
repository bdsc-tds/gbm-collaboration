import gc
import os
import time
import warnings

import decoupler as dc
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from typing import Dict, List, Set, Tuple, Optional
from pathlib import Path


def calculate_n_spots_per_sample(adata: sc.AnnData, layer=None) -> pd.DataFrame:
    """Calculate the number of non-zero spots for each gene within each sample.

    Args:
        adata (sc.AnnData): AnnData object containing single-cell/spatial data.
            Requires:
                - adata.obs['sample']: sample ID per cell/spot.
                - adata.layers[layer] or adata.X: expression matrix.
        layer (str or None): Which layer to use for counts (default uses adata.X).

    Returns:
        pd.DataFrame: Index = genes, columns = samples,
                      values = number of spots/cells with expression > 0.
    """

    samples = adata.obs['sample'].unique()
    n_spots_samples = pd.DataFrame(index=adata.var_names)

    for sample in samples:
        sample_idx = adata.obs['sample'] == sample
        sample_adata = adata[sample_idx]

        if layer is None:
            sample_counts = sample_adata.X
        else:
            sample_counts = sample_adata.layers[layer]

        # Handle sparse matrix
        with warnings.catch_warnings(action="ignore"):
            if scipy.sparse.issparse(sample_counts):
                n_spots_samples[sample] = (sample_counts > 0).sum(axis=0).A1
            else:
                n_spots_samples[sample] = (sample_counts > 0).sum(axis=0)


    return n_spots_samples


def gene_filter_by_expression(
    adata_case: Optional[sc.AnnData] = None,
    adata_ctrl: Optional[sc.AnnData] = None,
    layer: str = "log_norm_count",
    min_spots: int = 5,
    min_samples: int = 2,
) -> Tuple[Dict[str, pd.DataFrame], pd.Index]:
    """
    Filter genes based on expression presence across samples for one or both groups.

    If both adata_case and adata_ctrl are provided, keeps genes that meet the
    criteria in *both* groups (intersection). If only one is provided,
    it keeps genes that meet the criteria in that single group.

    Args:
        adata_case (Optional[sc.AnnData]): Case samples AnnData.
        adata_ctrl (Optional[sc.AnnData]): Control samples AnnData.
        layer (str): Which layer to use for counts.
        min_spots (int): Minimum expressing spots/cells per sample for a gene.
        min_samples (int): Minimum number of samples meeting the min_spots threshold.

    Returns:
        A tuple containing:
        - n_spots_dfs: A dictionary mapping group name ('case', 'ctrl') to its
          n_spots_per_sample DataFrame.
        - keep_genes: A pandas Index of the final filtered gene names.
    """
    if adata_case is None and adata_ctrl is None:
        raise ValueError("At least one of adata_case or adata_ctrl must be provided.")

    n_spots_dfs = {}
    pass_masks = []

    # Process case group if provided
    if adata_case is not None:
        n_spots_case = calculate_n_spots_per_sample(adata_case, layer=layer)
        n_spots_dfs['case'] = n_spots_case
        case_pass = (n_spots_case >= min_spots).sum(axis=1) >= min_samples
        pass_masks.append(case_pass)

    # Process control group if provided
    if adata_ctrl is not None:
        n_spots_ctrl = calculate_n_spots_per_sample(adata_ctrl, layer=layer)
        n_spots_dfs['ctrl'] = n_spots_ctrl
        ctrl_pass = (n_spots_ctrl >= min_spots).sum(axis=1) >= min_samples
        pass_masks.append(ctrl_pass)

    # Combine all available boolean masks with a logical AND
    # pd.concat aligns by index (genes) and .all(axis=1) checks for True in all columns
    final_pass_mask = pd.concat(pass_masks, axis=1).all(axis=1)

    keep_genes = final_pass_mask.index[final_pass_mask]

    return n_spots_dfs, keep_genes


def plot_gene_filtering_thresholds(
    n_spots_dfs: Dict[str, pd.DataFrame],
    max_min_spots: int = 20,
    max_min_samples: int = None
):
    """
    Plot number of genes retained as a function of min_spots and min_samples.

    Works for one or more groups provided in the n_spots_dfs dictionary.

    Args:
        n_spots_dfs (Dict[str, pd.DataFrame]): Dictionary from group name to its
            n_spots_per_sample DataFrame, as returned by gene_filter_by_expression.
        max_min_spots (int): Max min_spots threshold to test (1 to max_min_spots).
        max_min_samples (int): Max min_samples threshold to test. If None, uses
            the minimum number of samples across all provided groups.
    """
    if not n_spots_dfs:
        print("Warning: n_spots_dfs dictionary is empty. Nothing to plot.")
        return

    if max_min_samples is None:
        # Generalize to find the minimum number of samples across all provided groups
        max_min_samples = min(df.shape[1] for df in n_spots_dfs.values())

    min_spots_range = range(1, max_min_spots + 1)
    min_samples_range = range(1, max_min_samples + 1)
    results = []

    for min_spots in min_spots_range:
        for min_samples in min_samples_range:
            pass_masks = []
            # Calculate the passing mask for each group
            for n_spots_df in n_spots_dfs.values():
                group_pass = (n_spots_df >= min_spots).sum(axis=1) >= min_samples
                pass_masks.append(group_pass)

            # Combine masks and count the genes that pass in all groups
            final_mask = pd.concat(pass_masks, axis=1).all(axis=1)
            genes_retained = final_mask.sum()

            results.append({
                "min_spots": min_spots,
                "min_samples": min_samples,
                "genes_retained": genes_retained
            })

    df_results = pd.DataFrame(results)
    pivot_table = df_results.pivot(index="min_samples", columns="min_spots", values="genes_retained")

    plt.figure(figsize=(5, 4))
    c = plt.pcolor(pivot_table, cmap="viridis", edgecolors='k', linewidths=1)
    plt.colorbar(c, label='Number of genes retained')
    plt.xlabel("min_spots threshold")
    plt.ylabel("min_samples threshold")

    group_names = ", ".join(n_spots_dfs.keys())
    plt.title(f"Genes retained across group(s): {group_names}")

    plt.xticks(np.arange(0.5, len(pivot_table.columns), 1), pivot_table.columns)
    plt.yticks(np.arange(0.5, len(pivot_table.index), 1), pivot_table.index)
    plt.show()

def load_gene_sets(
    gene_sets_urls: Optional[Dict[str, str]] = None,
    gene_sets_files: Optional[Dict[str, str]] = None,
    gene_sets_dir: str = "gene_sets",
    valid_genes: Optional[Set[str]] = None,
    min_gene_set_size: int = 5,
    max_gene_set_size: int = 500,
    plot: bool = False,
    col_pathway = 'pathway_name',
    col_genes = 'gene',
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Downloads, loads, and filters gene sets from GMT (via URL) or CSV (local file)
    sources into pandas DataFrames.

    The filtered DataFrame contains lists of genes in a column.

    Args:
        gene_sets_urls: A dictionary mapping a short name to a GMT file URL.
        gene_sets_files: A dictionary mapping a short name to a local CSV file path.
                         The CSV must contain col_pathway and col_genes columns.
        gene_sets_dir: The directory to store downloaded GMT files.
        valid_genes: A set of gene symbols to keep in each gene set.
        min_gene_set_size: The minimum number of genes a set must have after filtering.
        max_gene_set_size: The maximum number of genes a set can have after filtering.
        plot: If True, displays a histogram of the filtered gene set sizes.

    Returns:
        A tuple containing:
        - filtered_df: A DataFrame of accepted gene sets with columns:
          ['source_fname', col_pathway, col_genes (list of strings)].
        - excluded_df: A DataFrame of excluded gene sets with columns:
          ['source_fname', col_pathway, 'filtered_size'].
    """
    if not gene_sets_urls and not gene_sets_files:
        raise ValueError("Either 'gene_sets_urls' or 'gene_sets_files' must be provided.")

    gene_sets_dir = Path(gene_sets_dir)
    gene_sets_dir.mkdir(exist_ok=True, parents=True)

    # 1. Create a unified list of all sources to process
    all_sources = []
    if gene_sets_urls:
        for name, url in gene_sets_urls.items():
            all_sources.append({'name': name, 'location': url, 'type': 'url'})
    if gene_sets_files:
        for name, path in gene_sets_files.items():
            all_sources.append({'name': name, 'location': path, 'type': 'file'})

    filtered_data = []
    excluded_data = []
    seen_pathways = set()

    # 2. Iterate over the unified list of sources
    for source in all_sources:
        fname = source['name']
        gene_sets_dict = {}

        # --- Loading Logic ---
        if source['type'] == 'url':
            path = gene_sets_dir / f'{fname}.gmt'
            if not path.exists():
                print(f"Downloading {fname} from {source['location']}...")
                subprocess.run(
                    ["wget", "--progress=bar:force", "-O", str(path), source['location']],
                    check=True,
                    capture_output=True
                )
            gene_sets_dict = gp.parser.read_gmt(path.as_posix())

        elif source['type'] == 'file':
            path = Path(source['location'])

            print(f"Loading local CSV {fname} from {path}...")
            # Assume CSV has col_pathway and col_genes columns
            df = pd.read_csv(path)
            gene_sets_dict = df.groupby(col_pathway)[col_genes].apply(list).to_dict()


        # --- Common Filtering Logic ---
        print(f"Loaded {len(gene_sets_dict)} gene sets from {fname}, filtering...")
        for name, genes in gene_sets_dict.items():
            if name in seen_pathways:
                print(f"Warning: Duplicate pathway '{name}' already seen. Keeping first instance.")
                continue
            seen_pathways.add(name)

            if valid_genes is not None:
                filtered_genes = [g for g in genes if g in valid_genes]
            else:
                filtered_genes = genes
            filtered_size = len(filtered_genes)

            if min_gene_set_size <= filtered_size <= max_gene_set_size:
                filtered_data.append([fname, name, filtered_genes])
            else:
                # Bug Fix: Appending the size, not the gene list, for excluded data
                excluded_data.append([fname, name, filtered_size])

    # Create the final DataFrames from the lists
    filtered_df = pd.DataFrame(
        filtered_data,
        columns=['source_fname', col_pathway, col_genes]
    )
    excluded_df = pd.DataFrame(
        excluded_data,
        columns=['source_fname', col_pathway, 'filtered_size']
    )

    print(f"\nKeeping {len(filtered_df)} gene sets.")
    print(f"Excluded {len(excluded_df)} gene sets after filtering.")

    if plot and not filtered_df.empty:
        sizes = filtered_df[col_genes].apply(len)
        sns.histplot(data=sizes, log_scale=True)
        plt.title('Filtered Gene Set Sizes')
        plt.xlabel('Number of Genes')
        plt.ylabel('Count')
        plt.show()

    return filtered_df, excluded_df


def score_gene_sets(
    adata: sc.AnnData,
    gene_sets: dict,
    min_size: int = 5,
    max_size: int = 500,
    layer: str = None,
    method: str = 'ulm',
    verbose=False,
    out_file=None,
    decoupler_storage='long',
    **kwargs
) -> pd.DataFrame:
    """Calculate scores using a specified method on adata.X. Assume expression data in adata.X is pre-transformed.

    Args:
        adata: An AnnData object with spots/cells as obs and genes as var.
        gene_sets: Gene sets as a dictionary {pathway_name: [gene_list]}.
        min_size: Minimum number of genes in a pathway/regulon.
        max_size: Maximum number of genes in a pathway/regulon.
        method: Scoring method to use 'ssgsea' or decoupler methods (e.g., 'ulm', 'mlm', 'viper', 'wmean',...).
        verbose: Print info or not
        out_file: Optional output file
        decoupler_storage (str): if using decoupler, store results in long or wide format
        **kwargs: Additional arguments passed to the scoring function.

    Returns:
        A pandas DataFrame with spots/cells as rows and pathway scores as columns.
    """
    if out_file is not None:
        out_file = Path(out_file)

    t_start = time.time()

    # 1. Run the specified method
    if method == 'ssgsea':
        print("Preparing data...")

        data = adata.to_df(layer=layer).T

        print(f"Running {method}...")
        ssgsea_results = gp.ssgsea(
            data=data,
            gene_sets=gene_sets,
            min_size=min_size,
            max_size=max_size,
            threads=1,
            verbose=verbose,
            **kwargs
        )
        df = ssgsea_results.res2d.set_index('Name')

    else:
        print("Preparing data...")
        # For decoupler, pre-filter gene_sets by size as it lacks a max_size argument
        filtered_gene_sets = {
            pathway: genes
            for pathway, genes in gene_sets.items()
            if min_size <= len(genes) <= max_size
        }

        if not filtered_gene_sets:
            print("  -> Skipping, no gene sets remained after size filtering.")
            return pd.DataFrame()

        # Convert the filtered signature dictionary to the required DataFrame format
        net_df = pd.Series(filtered_gene_sets).explode().reset_index()
        net_df.columns = ['source', 'target']
        net_df['weight'] = 1.0

        # Prepare anndata object with layer as .X
        X = adata.X if layer is None else adata.layers[layer]
        ad_dc = sc.AnnData(X)
        ad_dc.obs_names = adata.obs_names
        ad_dc.var_names = adata.var_names

        # Run decoupler method. min_n is still useful for post-filtering
        # based on genes actually present in the data.
        run_method_func = getattr(dc, f'run_{method}')
        run_method_func(
            ad_dc,
            net=net_df,
            source='source',
            target='target',
            weight='weight',
            use_raw=False,
            min_n=min_size, #tmin in newer decoupler versions
            verbose=verbose,
            **kwargs
        )
        obsm = ad_dc.obsm
        keys = list(obsm.keys())

        if decoupler_storage == 'long':
            print("Reshaping output...")
            # melt decoupler results into long dataframe

            df = obsm[keys[0]].melt(ignore_index=False, var_name='Term', value_name=keys[0])

            # if pvalues matrix also there, add as extra column
            if len(keys) > 1:
                assert obsm[keys[0]].shape == obsm[keys[1]].shape # assume same shape (--> same melting index order)
                df[keys[1]] = obsm[keys[1]].melt(ignore_index=False, var_name='Term', value_name=keys[1])[keys[1]]
            df.index.name ='Name'

    if out_file:
        print("Saving output...")
        if decoupler_storage =='wide' and method != 'ssgsea':
            for k in keys:
                df = ad_dc.obsm[k]
                out_file_ = out_file.with_stem(f'{out_file.stem}_{k}')
                df.to_parquet(out_file_)
        else:
            df.to_parquet(out_file)

    t_end = time.time()
    print(f"  -> Finished in {t_end - t_start:.2f} sec.")

    return df


def load_visium_adata(
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
    list_sample_id=None,
    normalization="log_norm_count",
    radius=1,
    deconvolution_slot_name="deconv_cohort_level_4",
    ):
    """Load entire Visium cohort as anndata."""
    from gbmcollaboration import CHUVDataCenter, UKERDataCenter

    # Define data centers and common loading parameters
    data_centers = [CHUVDataCenter, UKERDataCenter]

    # Load all data
    loaded_objects = []
    for data_center in data_centers:
        print('Loading',data_center)

        loaded_objects.append(data_center.load_spatial(
            exclude_low_quality_visium=exclude_low_quality_visium,
            exclude_idh1_mut=exclude_idh1_mut,
            list_sample_id=list_sample_id,
            normalization=normalization,
            radius=radius,
            deconvolution_slot_name=deconvolution_slot_name,
        ))

    # Concatenate the dataframes
    counts_df = pd.concat([obj.genecounts_df for obj in loaded_objects])
    obs_df = pd.concat([obj.celltype_df for obj in loaded_objects])

    del loaded_objects
    gc.collect()

    # Align obs to counts and create the 'sample' column
    obs_df = obs_df.loc[counts_df.index]
    obs_df['sample'] = obs_df.index.get_level_values(0)

    # Flatten the MultiIndex for AnnData compatibility
    flat_index = counts_df.index.map(lambda idx: f"{idx[0]}_{idx[1]}")
    obs_df.index = flat_index
    counts_df.index = flat_index

    counts_df.sum(axis=0).sort_values()

    # Create the final AnnData object
    adata = sc.AnnData(
        obs=obs_df,
        var=pd.DataFrame(index=counts_df.columns),
        layers={normalization:counts_df.values}
    )

    del counts_df, obs_df
    gc.collect()

    # --- Output ---
    print(" Final adata created:")
    print(" • shape:", adata.shape)
    print(" • layers:", list(adata.layers.keys()))
    print(" • obs columns:", adata.obs.columns.tolist())
    print(" • samples:", adata.obs['sample'].unique().tolist())

    return adata
