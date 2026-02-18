<div align="center">
<h1>GBM Collab repository</h1>
Data exploration and analyses of the GBM Mosaic data
</div>

# Table of Contents
- [Installation](#installation)
- [Load clean data](#load-clean-data)

# Installation
1. Start by initialising Conda.
```bash
conda init bash && bash
```

2. Create the environment.
```bash
conda create -n gbmcollaboration -y python=3.10
conda activate gbmcollaboration
```

3. Clone the repository and install the dependencies.
```bash
git clone https://github.com/owkin/gbm-collaboration.git
cd gbm-collaboration
make install-poetry
make install-all
pre-commit install
```

# Load clean data
Each modality loader constains arguments to filter the data based on patient's characteristics.
For all tabular data and single-cell data, those 3 main arguments can be used:
- exclude_idh1_mut: bool,
- exclude_patient_ids: bool,
- exclude_relapse: bool
A particularity for clinical data, is the exclusion of additionnal CR samples for which
OS and/or PFS are pending.

As of February 2026, If you want to load the clinical data, you can run the following code,
that allows flexibility in excluding particular samples at Charite.
```python
from gbmcollaboration import CHUVDataCenter, UKERDataCenter, CRDataCenter, SURVIVAL_EXCLUSION

# -- Loop over the three data centers and merge the clinical data
merged_df = []

for center in [CHUVDataCenter, UKERDataCenter, CRDataCenter]:
    center_data_dict = center.load_tabular(
        exclude_idh1_mut=True, # exclude IDH1/2 mutated patients
        exclude_patient_ids=True, # exclude sample completely after clinical checking, maybe reintroduce after initial submission
        exclude_relapse=True # exclude relapse samples
    )
    clin_dict = center_data_dict['clinical']
    df = clin_dict['processed gbm clinical']
    # Exclude CR samples for which survival data is pending, but rest of modalities can be used
    if center == CRDataCenter:
        df = df.drop(index=SURVIVAL_EXCLUSION)
    # To be sure check if we only have Baseline samples
    print(f"---Loading the clinical data from {center.name}---")
    print(df['Sample collection chronology'].value_counts())
    print(f"Found {len(df)} primary cases")
    # append to the merged dataframe
    merged_df.append(df)

# merge the three dataframes
df_concat = pd.concat(merged_df, ignore_index=False)
print("Dataframes from the three centers have been merged.")
```

Same arguments can be found when loading single-cells
```python
from gbmcollaboration import AllDataCenters
adata = AllDataCenters.load_singlecell(
    exclude_idh1_mut=True,
    exclude_patient_ids = True,
    exclude_relapse = True,
    )
```

For Visium, low quality samples can be additionnaly excluded, for example:
```python
from gbmcollaboration import AllDataCenters
adata = AllDataCenters.load_spatial_adata(
    exclude_low_quality_visium=True,
    exclude_idh1_mut=True,
    exclude_patient_ids = True,
    exclude_relapse = True
)
```

# Additional documentation
Additional information on how to run the InferCNVpy analyses can be found in this [README](./analysis/infercnv/README.md)
