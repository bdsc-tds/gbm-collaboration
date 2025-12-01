#!/bin/bash

conda config --add channels bioconda
conda config --add channels conda-forge
conda env create --name nf-core -f nf-core.yml -y
conda activate nf-core

# install singularity
# sudo apt update
# sudo apt install -y software-properties-common
# sudo add-apt-repository -y ppa:apptainer/ppa
# sudo apt update
# sudo apt install -y apptainer-suid

 # need to do this out of sagemaker then copy over
nf-core pipelines download hlatyping -r 2.1.0 -x none -s singularity -o ./nf-core-hlatyping_2.1.0/ # or -s singularity, -s docker, ...
# nextflow pull nf-core/hlatyping -r 2.1.0 # alternative download command. Works on sagemaker but w/o singularity
# nextflow plugin install nf-amazon@2.15.0  # copy to ~/.nextflow/plugins/nf-amazon-2.15.0. Only needed if using s3 paths in pipeline samplesheet
# nextflow plugin install nf-schema@2.3.0  # copy to ~/.nextflow/plugins/nf-schema-2.3.0. Only needed if using offline cluster

# put this in ~/.bashrc
# export NXF_OFFLINE='true'
