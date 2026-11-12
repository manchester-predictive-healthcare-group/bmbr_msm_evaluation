#!/bin/bash --login
#SBATCH -p multicore_small
#SBATCH -t 1-00:00:00
#SBATCH -n 4

set -euo pipefail
module purge
module load apps/gcc/R/4.5.0
Rscript p6.1_assess_calibration_csh.R > p6.1_assess_calibration_csh.out 2>&1
