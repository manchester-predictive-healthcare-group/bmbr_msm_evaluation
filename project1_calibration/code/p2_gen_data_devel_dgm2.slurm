#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 1-0
#SBATCH -a 1-2
#SBATCH -n 1

module purge
module load apps/gcc/R/4.5.0

set -euo pipefail

Rscript p2_gen_data_devel.R 2 ${SLURM_ARRAY_TASK_ID} > p2_gen_data_devel_dgm2.${SLURM_ARRAY_TASK_ID}.out