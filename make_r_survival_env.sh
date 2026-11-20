#!/bin/bash

conda activate gbm_main
conda install r-base rpy2 -y
Rscript -e 'install.packages("survRM2", repos="https://cloud.r-project.org")'
