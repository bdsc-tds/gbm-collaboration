# %%
# Autoreload extension for IPython/Jupyter
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('load_ext', 'autoreload')
        ipython.run_line_magic('autoreload', '2')
except (ImportError, NameError):
    pass

# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import s3fs
from cloudpathlib import S3Path
from scipy.sparse import issparse
from joblib import Parallel, delayed
from tqdm import tqdm
from gbmcollaboration.definitions import DATA_PATH_CHUV_FILES, DATA_PATH_UKER_FILES, DATA_PATH_CR_FILES

# %%
# Function to load one sample
def load_sample(center_name, sample_id, h5ad_file):
    """Load and process a single sample using Cell2LocationLoader approach."""
    try:
        # Load h5ad file
        if isinstance(h5ad_file, S3Path) or str(h5ad_file).startswith('s3://'):
            s3 = s3fs.S3FileSystem(anon=False)
            h5ad_path_str = str(h5ad_file).replace('s3://', '')
            with s3.open(h5ad_path_str, 'rb') as f:
                adata = ad.read_h5ad(f)
        else:
            adata = ad.read_h5ad(h5ad_file)
        
        # Extract means_cell_abundance_w_sf from obsm
        if 'means_cell_abundance_w_sf' not in adata.obsm:
            return None
        
        # Get cell type names
        celltype_names = None
        if 'mod' in adata.uns and 'factor_names' in adata.uns['mod']:
            celltype_names = adata.uns['mod']['factor_names']
        elif 'factor_names' in adata.uns:
            celltype_names = adata.uns['factor_names']
        elif 'cell_type_names' in adata.uns:
            celltype_names = adata.uns['cell_type_names']
        
        if celltype_names is None:
            n_celltypes = adata.obsm['means_cell_abundance_w_sf'].shape[1]
            celltype_names = [f'celltype_{i}' for i in range(n_celltypes)]
        
        # Add prefix to match column names
        celltype_names = [f'meanscell_abundance_w_sf_{name}' for name in celltype_names]
        
        # Get cell abundance data
        cell_abundance = adata.obsm['means_cell_abundance_w_sf']
        if issparse(cell_abundance):
            cell_abundance = cell_abundance.toarray()
        
        # Create DataFrame
        df = pd.DataFrame(cell_abundance, index=adata.obs_names, columns=celltype_names)
        df = df.fillna(0)
        
        # Add sample ID to index
        df.index = pd.MultiIndex.from_tuples(
            [(sample_id, idx) for idx in df.index],
            names=['sample', 'barcode']
        )
        
        # Compute total cells per spot
        df['total_cells'] = df.sum(axis=1)
        
        return df
    except Exception as e:
        return None

# %%
# Test: Load one sample
import time

test_center = 'CHUV'
test_cell2loc_path = DATA_PATH_CHUV_FILES["Cell2location"]
test_base_path = S3Path(str(test_cell2loc_path))
test_sample_dirs = [d.name for d in test_base_path.iterdir() if d.is_dir()]
test_sample_id = test_sample_dirs[0]
test_h5ad_file = test_base_path / test_sample_id / "cohort_level_2.h5ad"

sns.histplot(test_df['total_cells'])

# %%
# Load all centers, all samples per center in parallel
centers = {
    'CHUV': DATA_PATH_CHUV_FILES["Cell2location"],
    # 'UKER': DATA_PATH_UKER_FILES["Cell2location"],
    # 'CR': DATA_PATH_CR_FILES["Cell2location"],
}

all_results = []

for center_name, cell2loc_path in centers.items():
    print(f"\nProcessing {center_name}...")
    base_path = S3Path(str(cell2loc_path))
    sample_dirs = [d.name for d in base_path.iterdir() if d.is_dir()]
    
    # Prepare sample paths
    samples_to_process = [
        (center_name, sample_id, base_path / sample_id / "cohort_level_2.h5ad")
        for sample_id in sample_dirs
    ]
    
    # Process samples in parallel with tqdm for this center
    center_results = Parallel(n_jobs=-1)(
        delayed(load_sample)(center_name, sample_id, h5ad_file)
        for center_name, sample_id, h5ad_file in tqdm(samples_to_process, desc=f"  {center_name}")
    )
    
    # Filter out None results and add center info
    center_results = [df for df in center_results if df is not None]
    
    for df in center_results:
        # Add center to index
        new_index = pd.MultiIndex.from_tuples(
            [(center_name, *idx) for idx in df.index],
            names=['center', 'sample', 'barcode']
        )
        df.index = new_index
        df['center'] = center_name
    
    all_results.extend(center_results)
    print(f"  ✓ Loaded {len(center_results)} samples from {center_name}")

# Combine all results
df_all = pd.concat(all_results, axis=0)
print(f"\n✅ Loaded {len(df_all)} spots from {len(all_results)} samples across all centers")

# %%
# Plot histogram
sns.histplot(data=df_all, x='total_cells', bins=50)
plt.xlabel('Total Cells per Spot')
plt.ylabel('Frequency')
plt.title('Cell Count Distribution per Spot - All Centers')
plt.axvline(df_all['total_cells'].median(), color='red', linestyle='--', 
            label=f'Median: {df_all["total_cells"].median():.1f}')
plt.legend()
plt.show()

print(f"\nSummary:")
print(f"Total spots: {len(df_all)}")
print(f"Mean cells per spot: {df_all['total_cells'].mean():.2f}")
print(f"Median cells per spot: {df_all['total_cells'].median():.2f}")
print(f"\nBy center:")
for center in df_all['center'].unique():
    center_data = df_all[df_all['center'] == center]
    print(f"  {center}: {len(center_data)} spots, median: {center_data['total_cells'].median():.2f}")
