# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import h5py
import scipy
import scanpy as sc
from cloudpathlib import S3Path
from scipy.sparse import issparse
from joblib import Parallel, delayed
from tqdm.auto import tqdm
from cell2location import run_colocation
from gbmcollaboration.utils.io_utils import s3_read_h5ad_partial
from gbmcollaboration.definitions import DATA_PATH_CHUV_FILES, DATA_PATH_UKER_FILES, DATA_PATH_CR_FILES
from gbmcollaboration.definitions import (
    PATIENTS_IDH1_MUT,
    BAD_QUALITY_CHUV_VISIUM,
    BAD_QUALITY_CR_VISIUM,
    BAD_QUALITY_UKER_VISIUM,
    PATIENT_EXCLUSION,
    RELAPSE_EXCLUSION,
    SURVIVAL_EXCLUSION
)


def merge_cell2location_dicts(dict_list):
    # 1. Take intersection of gene names and factor names
    # Note: np.intersect1d automatically sorts the resulting array alphabetically
    common_vars = dict_list[0]['var_names']
    common_facts = dict_list[0]['factor_names']

    for d in dict_list[1:]:
        common_vars = np.intersect1d(common_vars, d['var_names'])
        common_facts = np.intersect1d(common_facts, d['factor_names'])

    # 2. Build the top-level keys
    result = {
        'date': dict_list[0].get('date'),
        'model_name': dict_list[0].get('model_name'),
        'factor_filter': dict_list[0].get('factor_filter'),
        'var_names': common_vars,
        'factor_names': common_facts,
        # Concatenate cell barcodes
        'obs_names': np.concatenate([d['obs_names'] for d in dict_list], axis=0)
    }

    # 3. Explicitly process the known dictionary groups
    post_sample_keys =[
        'post_sample_means', 'post_sample_q05',
        'post_sample_q95', 'post_sample_stds'
    ]

    for psk in post_sample_keys:
        if psk not in dict_list[0]:
            continue

        result[psk] = {}
        inner_keys = dict_list[0][psk].keys()

        for inner_k in inner_keys:
            arrays_to_concat =[]

            for d in dict_list:
                arr = d[psk][inner_k]
                orig_vars = d['var_names']
                orig_facts = d['factor_names']

                # Check if this is a matrix that needs column subsetting
                if arr.ndim > 1:
                    # If columns match the original number of genes
                    if arr.shape[1] == len(orig_vars):
                        # Get indices of shared genes in the original array
                        _, orig_idx, _ = np.intersect1d(orig_vars, common_vars, return_indices=True)
                        arr = arr[:, orig_idx] # Subset columns

                    # If columns match the original number of factors
                    elif arr.shape[1] == len(orig_facts):
                        # Get indices of shared factors in the original array
                        _, orig_idx, _ = np.intersect1d(orig_facts, common_facts, return_indices=True)
                        arr = arr[:, orig_idx] # Subset columns

                arrays_to_concat.append(arr)

            # Concatenate rows together
            result[psk][inner_k] = np.concatenate(arrays_to_concat, axis=0)

    return result

def load_sample(center_name, sample_id, h5ad_file):
    """Load and process a single sample using Cell2LocationLoader approach."""
    adata = s3_read_h5ad_partial(h5ad_file,uns=True,obs=True,var=True,obsm=True)

    if 'means_cell_abundance_w_sf' not in adata.obsm:
        return None
    adata.obs['center'] = center_name
    adata.obs['sample'] = sample_id
    return adata

# %%
EXCLUDE_SAMPLES = (
    PATIENTS_IDH1_MUT+
    BAD_QUALITY_CHUV_VISIUM+
    BAD_QUALITY_CR_VISIUM+
    BAD_QUALITY_UKER_VISIUM+
    PATIENT_EXCLUSION+
    RELAPSE_EXCLUSION+
    SURVIVAL_EXCLUSION
)
path_efs = '/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/'
path_ecm_annot = path_efs+'Analysis_Repo/Notebooks/ecm_neighborhood_annotations_v4.csv'
df_ecm_annot = pd.read_csv(path_ecm_annot,index_col='spot_id')
u_samples = df_ecm_annot['sample'].unique()

# %%
# Load all centers, all samples per center in parallel
centers = {
    'CHUV': DATA_PATH_CHUV_FILES["Cell2location"],
    'UKER': DATA_PATH_UKER_FILES["Cell2location"],
    'CR': DATA_PATH_CR_FILES["Cell2location"],
}

all_results = []

for center_name, cell2loc_path in centers.items():
    print(f"\nProcessing {center_name}...")
    base_path = S3Path(str(cell2loc_path))
    sample_dirs = [d.name for d in base_path.iterdir() if d.is_dir()]

    # Prepare sample paths
    samples_to_process = [
        (center_name, sample_id, base_path / sample_id / "cohort_level_4.h5ad")
        for sample_id in sample_dirs
        if sample_id in u_samples
        # if not any(ign in sample_id for ign in EXCLUDE_SAMPLES)
    ]

    # Process samples in parallel with tqdm for this center
    # center_results = []
    # for center_name, sample_id, h5ad_file in tqdm(samples_to_process, desc=f"  {center_name}"):
    #     center_results.append(load_sample(center_name, sample_id, h5ad_file))
    center_results = Parallel(n_jobs=20)(
        delayed(load_sample)(center_name, sample_id, h5ad_file)
        for center_name, sample_id, h5ad_file in tqdm(samples_to_process, desc=f"  {center_name}")
    )

    # Filter out None results and add center info
    center_results = [ad_ for ad_ in center_results if ad_ is not None]
    all_results.extend(center_results)

# %%
adata_vis = sc.concat(all_results,join='inner')
adata_vis.obs_names = adata_vis.obs['sample'] + '_' + adata_vis.obs_names
adata_vis.uns['mod'] = merge_cell2location_dicts([ad_.uns['mod'] for ad_ in all_results])
assert set(adata_vis.obs_names) == set(df_ecm_annot.index)

# %%
run_name = 'means_cell_abundance_w_sf'
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(5, 30), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3, # number of training restarts
      "cell_abundance_obsm": "means_cell_abundance_w_sf",

    },
    # the hyperparameters of NMF can be also adjusted:
    model_kwargs={'alpha': 0.01, 'init': 'random', "nmf_kwd_args": {"tol": 0.000001}},
    export_args={'path': path_efs / f'results/visium/cell2loc/{run_name}/CoLocatedComb/'}
)
# %%
import cell2location.models as c2l
# create model class
n_fact = int(5)
# X_data = adata_vis.uns['mod']['post_sample_q05']['spot_factors']
X_data = adata_vis.obsm['means_cell_abundance_w_sf'].values
mod_sk = c2l.CoLocatedGroupsSklearnNMF(n_fact, X_data,
        n_iter = 10000,
        verbose = True,
        var_names=adata_vis.uns['mod']['factor_names'],
        obs_names=adata_vis.obs_names,
        fact_names=['fact_' + str(i) for i in range(n_fact)],
        sample_id=adata_vis.obs['sample'],
        init='random', random_state=0,
        nmf_kwd_args={'tol':0.0001})

# train 5 times to evaluate stability
mod_sk.fit(n=5, n_type='restart')
