#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 04:00:00
#SBATCH -n 1

set -euo pipefail
module purge
module load apps/gcc/R/4.5.0
Rscript p6.3_create_combined_calib_figures.R > p6.3_create_combined_calib_figures.out 2>&1
