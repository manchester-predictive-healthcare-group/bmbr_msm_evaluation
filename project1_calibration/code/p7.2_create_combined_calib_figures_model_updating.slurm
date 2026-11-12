#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 04:00:00
#SBATCH -n 1

set -euo pipefail
module purge
module load apps/gcc/R/4.5.0
Rscript p7.2_create_combined_calib_figures_model_updating.R > p7.2_create_combined_calib_figures_model_updating.out 2>&1
