#!/bin/bash -l

# Check if modality is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <modality> (e.g., WES, RNAseq, etc.)"
  exit 1
fi

MODALITY=$1

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR=../../results/${MODALITY}/hlatyping/
OUT_DIR=$RESULTS_DIR/results
WORK_DIR=$RESULTS_DIR/work
INPUT_SAMPLESHEET=$RESULTS_DIR/samplesheet_copy_with_duplicates_G.csv

# Activate conda environment
conda activate nf-core
echo "Activated Nextflow environment"

# Run Nextflow pipeline
nextflow run nf-core-hlatyping_2.1.0/2_1_0/ \
    -profile slurm \
    -resume \
    -work-dir $WORK_DIR \
    --input $INPUT_SAMPLESHEET \
    --outdir $OUT_DIR
