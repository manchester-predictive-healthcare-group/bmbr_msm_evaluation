#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 1-00:00:00
#SBATCH -a 1-4
#SBATCH -n 1

set -euo pipefail
module purge
module load apps/gcc/R/4.5.0
Rscript "p4_fit_model_dgm${SLURM_ARRAY_TASK_ID}.R" > "p4_fit_model_dgm${SLURM_ARRAY_TASK_ID}.out" 2>&1
