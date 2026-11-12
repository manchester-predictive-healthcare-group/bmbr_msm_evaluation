#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 04:00:00
#SBATCH -n 1

set -euo pipefail
module purge
module load apps/gcc/R/4.5.0
Rscript p5.2_pred_tp_combine.R > p5.2_pred_tp_combine.out 2>&1
