
import pandas as pd
import numpy as np
import squidpy as sq
import scipy.stats as st
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from pygam import LinearGAM, s, l, te
from joblib import Parallel, delayed
from tqdm.auto import tqdm

def agg_neighbor_values(
    A,
    features,
    weighted=False,
    dummy_na=False,
):
    """
    Aggregate neighbor values using a spatial adjacency matrix.

    Parameters
    ----------
    A : scipy sparse matrix (n_cells x n_cells)
        Spatial adjacency matrix.
    features : pd.Series or pd.DataFrame
        Cell-level features to aggregate.
        - Series (numeric): summed
        - Series (categorical): one-hot then summed
        - DataFrame: summed column-wise
    weighted : bool
        Whether to use edge weights in A.
    dummy_na : bool
        Whether to include NaN as category (for categorical Series).

    Returns
    -------
    pd.DataFrame (n_cells x n_features)
        Neighbor-aggregated features.
    """
    A = A.copy()
    if not weighted:
        A.data = np.ones_like(A.data)

    # --- Prepare feature matrix L ---
    if isinstance(features, pd.Series):
        if pd.api.types.is_numeric_dtype(features):
            L = features.to_frame()
        else:
            L = pd.get_dummies(features, dummy_na=dummy_na)

    elif isinstance(features, pd.DataFrame):
        L = features

    else:
        raise TypeError("features must be a pandas Series or DataFrame")

    # --- Matrix multiplication ---
    out = A @ L.values

    return pd.DataFrame(
        out,
        index=L.index,
        columns=L.columns
    )

def find_spatial_clusters(
    adata,
    score_key,
    sample_key='sample',
    coord_type='generic',
    obsp="spatial_distances",
    radius=20,
    min_cluster_size=30,
    score_threshold=1,
    score_mode='sum',
    score_cats=None,
    prefix='TLS',
    recompute_neighbors=False,
    set_diag=True,
):
    score_col = f'{prefix}_{score_mode}'
    label_col = f'{prefix}_label'
    split_col = f'{prefix}_split'
    filtered_col = f'{prefix}_filtered'

    # Get score based on neighbor counts or scores
    if recompute_neighbors or obsp not in adata.obsp:
        adata.obs[sample_key] = adata.obs[sample_key].astype('category')
        sq.gr.spatial_neighbors(adata, radius=radius, library_key=sample_key, coord_type=coord_type, set_diag=True)
        A = adata.obsp[obsp].copy()
    else:
        # Remove entries > radius
        A = adata.obsp[obsp].copy()
        A.data[A.data > radius] = 0.
        A.eliminate_zeros()

    if set_diag:
        # Add back diagonal to distances (set_diag=True only adds it to connectivities)
        A.setdiag(0.)

    nbh_values = agg_neighbor_values(A, adata.obs[score_key])
    if not pd.api.types.is_numeric_dtype(adata.obs[score_key]):
        if score_cats is None:
            score_cats = nbh_values.columns
        adata.obs[score_col] = nbh_values[score_cats].agg(score_mode,axis=1)
    else:
        adata.obs[score_col] = nbh_values

    # Binarize based on threshold
    if isinstance(score_threshold,str):
        score_threshold = adata.obs[score_col].agg(score_threshold)
    mask = (adata.obs[score_col] >= score_threshold).values

    adata.obs[label_col] = np.where(mask, prefix, f"non-{prefix}")
    adata.obs[split_col] = f"non-{prefix}"
    adata.obs[filtered_col] = f"non-{prefix}"

    # Split based on spatial connectivity
    if mask.any():
        sub_adj = A[mask][:,mask]
        n_components, labels = connected_components(sub_adj, directed=False)
        adata.obs.loc[mask,split_col] = labels.astype(str)

        # Filter small clusters
        unique_ids, counts = np.unique(labels, return_counts=True)
        big_cluster_ids = unique_ids[counts >= min_cluster_size]

        final_labels = np.where(np.isin(labels, big_cluster_ids), labels.astype(str), f"non-{prefix}")
        adata.obs.loc[mask, filtered_col] = final_labels


class SpatialTissueAnalyzer:
    def __init__(self, radius=100, n_perms=50, n_jobs=-1):
        self.radius = radius
        self.n_perms = n_perms
        self.n_jobs = n_jobs

    def _prepare_data(self, df_sample, source_col='is_TIL', target_col='ECM_score', rank=True, smooth_target=True, adj_matrix=None):
        """
        One-pass spatial aggregation.
        Stores raw counts for contrast and ranked values for regression.
        """
        coords = df_sample[['x', 'y']].values

        # Use precomputed adjacency matrix if provided
        if adj_matrix is not None:
            if not isinstance(adj_matrix, csr_matrix):
                raise ValueError("adj_matrix must be a scipy.sparse.csr_matrix")
            adj_indices = [list(adj_matrix[i].indices) for i in range(adj_matrix.shape[0])]
        else:
            # Compute adjacency using KDTree
            tree_all = cKDTree(coords)
            adj_indices = tree_all.query_ball_tree(tree_all, self.radius)

        # Aggregate Predictor
        src_values = df_sample[source_col].values
        src_raw = np.array([np.sum(src_values[idx]) for idx in adj_indices])

        # Aggregate Target
        tgt_values = df_sample[target_col].values
        if smooth_target:
            tgt_raw = np.array([np.mean(tgt_values[idx]) for idx in adj_indices])
        else:
            tgt_raw = tgt_values

        # Store
        df_sample['source_raw'] = src_raw
        df_sample['target_raw'] = tgt_raw

        if rank:
            df_sample['source_final'] = df_sample['source_raw'].rank(pct=True)
            df_sample['target_final'] = df_sample['target_raw'].rank(pct=True)
        else:
            df_sample['source_final'] = df_sample['source_raw']
            df_sample['target_final'] = df_sample['target_raw']

        return df_sample, adj_indices

    @staticmethod
    def _spearman_perm(seed_i, raw_labels, adj_indices, y):
        rs = np.random.default_rng(seed_i)
        perm_labels = rs.permutation(raw_labels)
        fake_counts = np.array([np.sum(perm_labels[idx]) for idx in adj_indices])
        rho, _ = st.spearmanr(fake_counts, y)
        return rho

    def analyze_sample(self, df_sample, source_col='is_TIL', target_col='ECM_score',
                       method="gam_contrast", n_splines=25, rank=True, smooth_target=True, fit_null=False, adj_matrix=None, seed=0):

        # Prepare data and retrieve the neighbor index list
        if method == 'spearman_spatial_perm':
            smooth_target = False

        df, adj_indices = self._prepare_data(
            df_sample.copy(), source_col, target_col,
            rank=rank, smooth_target=smooth_target,adj_matrix=adj_matrix)

        N = len(df)
        y = df['target_final'].values
        coords = df[['x', 'y']].values

        if method == "spearman_spatial_perm":
            source_raw = df['source_raw'].values
            target_raw = df['target_raw'].values

            obs_rho, _ = st.spearmanr(source_raw, target_raw ,nan_policy='omit')
            null_rhos = []

            # Use raw labels for permutation to preserve total count

            print(f"Running {self.n_perms} permutations...")
            # Generate independent seeds for reproducibility
            base_rng = np.random.default_rng(seed)
            perm_seeds = base_rng.integers(0, 2**32 - 1, size=self.n_perms)

            if self.n_jobs == 1:
                null_rhos = [
                    self._spearman_perm(s, source_raw, adj_indices, target_raw)
                    for s in tqdm(perm_seeds)
                ]
            else:
                null_rhos = Parallel(n_jobs=self.n_jobs)(
                    delayed(self._spearman_perm)(s, source_raw, adj_indices, target_raw)
                    for s in perm_seeds
                )

            null_rhos = np.array(null_rhos)

            p_val = (np.sum(np.abs(null_rhos) >= np.abs(obs_rho)) + 1) / (self.n_perms + 1)
            return_ = {'beta': np.arctanh(obs_rho), 'rho': obs_rho, 'se': 1/np.sqrt(N-3), 'p': p_val, 'n': N}

        elif method == "gam_spatial":
            X_mat = np.column_stack([df['source_final'].values, coords])
            gam = LinearGAM(l(0) + te(1, 2, n_splines=n_splines)).fit(X_mat, y)

            if fit_null:
                gam_null = LinearGAM(te(1, 2, n_splines=n_splines)).fit(X_mat, y)

            return_ = {
                'beta': gam.coef_[0],
                'se': np.sqrt(gam.statistics_['cov'][0, 0]),
                'p': gam.statistics_['p_values'][0],
                'n': N,
                'gam': gam,
            }


        elif method == "gam_contrast":
            # Direct use of the pre-computed raw counts
            raw_counts = df['source_raw'].values

            is_none = (raw_counts == 0).astype(int)
            is_high = (raw_counts > 1).astype(int)

            if len(np.unique(is_high)) < 2: return None

            X_mat = np.column_stack([is_none, is_high, coords])
            # Intercept is exactly 1 TIL. l(1) is >1 TIL vs 1 TIL.
            gam = LinearGAM(l(0) + l(1) + te(2, 3, n_splines=n_splines)).fit(X_mat, y)
            return_ =  {
                'beta': gam.coef_[1],
                'se': np.sqrt(gam.statistics_['cov'][1, 1]),
                'p': gam.statistics_['p_values'][1],
                'n': N,
                'gam': gam,
                }

            if fit_null:
                gam_null = LinearGAM(te(2, 3, n_splines=n_splines)).fit(X_mat, y)

        if fit_null:
            # Calculate Deviance Explained
            r2_null = gam_null.score(X_mat, y)
            r2_full = gam.score(X_mat, y)
            r2_diff = r2_full - r2_spatial

            # Effect Multiplier (Relative Strength)
            strength_ratio = (r2_full - r2_spatial) / r2_spatial if r2_spatial > 0 else 0
            return_.update(
                {
                'gam_null':gam_null,
                'r2_null':r2_spatial,
                'r2_full':r2_full,
                'r2_diff':r2_diff,
                'relative_strength':strength_ratio
                }
            )

        return return_

    def run_meta(self, df, sample_col='sample', source_col='is_TIL', target_col='ECM_score', method='gam_contrast', n_splines=25,
     rank=True, smooth_target=True, fit_null=False, adj_matrix=None, seed=0):
        self.method = method

        sample_results = []
        for sid in np.unique(df[sample_col]):
            print(f"\nAnalyzing {sid}")

            idx_sample = (df[sample_col]==sid).values
            group = df[idx_sample]

            if adj_matrix is not None:
                sadj_matrix = adj_matrix[idx_sample][:,idx_sample]
            else:
                sadj_matrix= None

            res = self.analyze_sample(group, source_col=source_col, target_col=target_col,
                                      method=method, n_splines=n_splines,rank=rank, smooth_target=smooth_target,
                                      fit_null=fit_null,seed=seed, adj_matrix=sadj_matrix)
            if res:
                res[sample_col] = sid
                sample_results.append(res)

        res_df = pd.DataFrame(sample_results)
        meta_dict = self._meta_pool(res_df)
        return res_df, meta_dict

    def _meta_pool(self,res_df):
        # DerSimonian-Laird Random Effects Meta-Analysis
        y = res_df['beta'].values
        v = res_df['se'].values**2
        w = 1.0 / v

        # 1. Fixed Effect Pooled
        weighted_avg = np.sum(w * y) / np.sum(w)

        # 2. Between-sample variance (Tau^2)
        Q = np.sum(w * (y - weighted_avg)**2)
        df = len(y) - 1
        denom = np.sum(w) - (np.sum(w**2) / np.sum(w))
        tau2 = max(0, (Q - df) / denom)

        # 3. Random Effect weights
        w_star = 1.0 / (v + tau2)
        pooled_beta = np.sum(w_star * y) / np.sum(w_star)
        pooled_se = np.sqrt(1.0 / np.sum(w_star))
        p_val = 2 * (1 - st.norm.cdf(abs(pooled_beta / pooled_se)))

        return_ = {'beta': pooled_beta, 'p': p_val, 'tau2': tau2}
        if self.method == 'spearman_spatial_perm':
            return_['rho'] = np.tanh(pooled_beta)

        return return_



def plot_cross_platform_comparison(
    df,
    x='logFC_xenium',
    y='logFC_visium',
    pval_x='pvals_adj',
    pval_y='FDR',
    name_col='names',
    alpha=0.05,
    top_n = 25,
    title='Xenium vs Visium'
):
    df = df.copy()

    # 1. Create the 4-way classification
    sig_x = df[pval_x] < alpha
    sig_y = df[pval_y] < alpha

    conditions = [
        (sig_x & sig_y),   # Both
        (sig_x & ~sig_y),  # Xenium Only
        (~sig_x & sig_y),  # Visium Only
        (~sig_x & ~sig_y)  # Neither
    ]
    choices = ['Both', 'Xenium Only', 'Visium Only', 'Neither']
    df['Status'] = np.select(conditions, choices, default='Neither')

    # 2. Sizing: Use the best p-value (minimum) for sizing
    # -log10 of the smallest p-value makes the most significant points largest
    min_p = df[[pval_x, pval_y]].min(axis=1)
    df['minus_log10_p'] = -np.log10(min_p + 1e-300)

    # 3. Define the Plot
    plt.figure(figsize=(10, 8))

    # Define custom palette for clarity
    color_map = {
        'Both': '#8e44ad',        # Purple
        'Xenium Only': '#2ecc71', # Green
        'Visium Only': '#3498db', # Blue
        'Neither': '#bdc3c7'      # Grey
    }

    ax = sns.scatterplot(
        data=df,
        x=x,
        y=y,
        hue='Status',
        size='minus_log10_p',
        hue_order=['Both', 'Xenium Only', 'Visium Only', 'Neither'],
        palette=color_map,
        sizes=(40, 400),
        alpha=0.8,
        edgecolor='white',
        linewidth=0.5
    )

    # 4. Add labels for top_n per quadrant
    texts = []

    # Define which categories we want to label
    label_groups = ['Both', 'Xenium Only', 'Visium Only']

    for group in label_groups:
        # Filter for this specific quadrant
        group_df = df[df['Status'] == group].copy()

        if group_df.empty:
            continue

        # Sort by importance (e.g., most significant p-value)
        # We use the -log10_p we calculated earlier for sizing
        group_df = group_df.sort_values('minus_log10_p', ascending=False)

        # Take the top N
        df_top = group_df.head(top_n)

        for i in range(df_top.shape[0]):
            texts.append(
                plt.text(
                    df_top[x].iloc[i],
                    df_top[y].iloc[i],
                    df_top[name_col].iloc[i],
                    fontsize=6,
                    weight='bold' if group == 'Both' else 'normal',
                    # Optional: color text to match the point
                    color='black'
                )
            )

    if texts:
        # expand_text=(1.2, 1.2) helps push labels a bit further from points
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
            expand_text=(1.2, 1.2)
        )

    # 5. Styling
    plt.axhline(0, color='black', linestyle='--', alpha=0.3)
    plt.axvline(0, color='black', linestyle='--', alpha=0.3)

    plt.title(title, fontsize=15)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Significance")
    plt.tight_layout()
    plt.show()
