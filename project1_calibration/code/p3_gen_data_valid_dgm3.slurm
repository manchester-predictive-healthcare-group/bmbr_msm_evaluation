#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 1-0
#SBATCH -a 1-8
#SBATCH -n 1

module purge
module load apps/gcc/R/4.5.0

Rscript p3_gen_data_valid.R 3 ${SLURM_ARRAY_TASK_ID} > p3_gen_data_valid_dgm3.${SLURM_ARRAY_TASK_ID}.out