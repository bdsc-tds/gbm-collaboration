# %%
import subprocess
from pathlib import Path
import pandas as pd
from gbmcollaboration.utils import s3_ls, get_datasets

def get_optitype_samplesheet_s3(bucket, paths, seq_type='dna'):
    """Creates a DataFrame from a list of S3 paths."""
    df = []
    for p in paths:
        # Create full s3:// path
        s3_full_path = f's3://{bucket}/{p.as_posix()}'
        df.append(
            [
                p.stem,      # 'sample_name' from 'sample_name.bam'
                s3_full_path,
                seq_type
            ]
        )
    df = pd.DataFrame(df, columns=['sample', 'bam', 'seq_type'])
    return df

def sync_bam_files_s3(s3_prefix, local_dir):
    """
    Syncs all .bam files from an S3 prefix to a local directory using 'aws s3 sync'.
    """
    print(f"Syncing all .bam files from {s3_prefix} to {local_dir}...")
    local_dir.mkdir(parents=True, exist_ok=True)

    # The aws s3 sync command with filters
    # 1. Exclude everything
    # 2. Include only files ending in .bam
    command = [
        "aws", "s3", "sync",
        s3_prefix,
        str(local_dir),
        "--exclude", "*",
        "--include", "*.bam",
        "--no-follow-symlinks"
    ]

    try:
        # Using check=True will raise an exception if the command returns a non-zero exit code
        subprocess.run(command, check=True, capture_output=True, text=True)
        print("Sync complete.")
    except subprocess.CalledProcessError as e:
        print("--- AWS CLI Error ---")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        raise RuntimeError(f"Failed to sync files from {s3_prefix}. AWS CLI returned a non-zero exit code.") from e


# --- Configuration ---
DOWNLOAD_BAM_FILES = True
MODALITIES = ['WES']


for modality in MODALITIES:
    # --- Outputs ---
    RESULTS_DIR = Path(f'/home/sagemaker-user/custom-file-systems/efs/fs-09913c1f7db79b6fd/results/{modality}/hlatyping/')
    OUT_FILE_SAMPLESHEET_LOCAL = RESULTS_DIR / 'samplesheet_local.csv'
    OUT_FILE_SAMPLESHEET_S3 = RESULTS_DIR / 'samplesheet_s3.csv'

    # --- Main Logic ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    if modality == 'WES':
        seq_type = 'dna'
    elif modality == 'RNAseq':
        seq_type = 'rna'
    else:
        raise ValueError('modality should be WES or RNAseq')

    # 1. Get WES S3 paths and create initial samplesheet
    print("Gathering S3 paths for all datasets...")
    datasets = get_datasets()
    all_s3_paths = []
    dataset_prefixes = {}
                                                # not present on abstra but keep the loop in case future cohorts are added
    for dataset_name in ['MOSAIC - UKER - GBM', 'MOSAIC - CR - GBM']: #, 'MOSAIC - CHUV - GBM']:
        dataset_info = datasets[dataset_name]
        s3_prefix_str = f"s3://{dataset_info['bucket']}/{dataset_info['prefix']}{modality}/"

        # Store the prefix for the sync command later
        dataset_prefixes[dataset_name] = s3_prefix_str

        # List all bam files to build the initial dataframe
        files = s3_ls(dataset_info['bucket'], dataset_info['prefix'] + modality, recursive=True)
        files_bam = [Path(p.key) for p in files if p.key.endswith('.bam')]

        df_ = get_optitype_samplesheet_s3(dataset_info['bucket'], files_bam, seq_type=seq_type)
        all_s3_paths.append(df_)

    # 2. Save the s3 samplesheet
    df_s3 = pd.concat(all_s3_paths, ignore_index=True)
    df_s3.to_csv(OUT_FILE_SAMPLESHEET_S3, index=False)
    print(f"S3 samplesheet saved to {OUT_FILE_SAMPLESHEET_S3}")


    # 3. Sync all required BAM files from all datasets in one go
    if DOWNLOAD_BAM_FILES:
        print("\nStarting BAM file synchronization...")
        for dataset_name, s3_prefix in dataset_prefixes.items():
            sync_bam_files_s3(s3_prefix, RESULTS_DIR)


        # 4. Update the DataFrame with the new local paths
        print("\nUpdating samplesheet with local file paths...")
        local_bam_paths = []
        for s3_bam_path in df_s3['bam']:
            # The original S3 path is like 's3://bucket/prefix/to/file.bam'
            # We just need the filename part.
            file_name = Path(s3_bam_path).name
            local_path = RESULTS_DIR / f'bam/{file_name}'

            # A quick check to ensure the file was actually synced
            if not local_path.exists():
                raise FileNotFoundError(f"Synced file not found at expected location: {local_path}. Check the sync log.")

            local_bam_paths.append(str(local_path))

        # 5. Save the local samplesheet
        df_local = df_s3.copy()
        df_local['bam'] = local_bam_paths
        df_local.to_csv(OUT_FILE_SAMPLESHEET_LOCAL, index=False)
        print(f"Local samplesheet saved to {OUT_FILE_SAMPLESHEET_LOCAL}")

# %%
