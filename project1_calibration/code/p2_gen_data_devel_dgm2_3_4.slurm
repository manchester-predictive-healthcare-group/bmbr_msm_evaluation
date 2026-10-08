#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 1-0
#SBATCH -a 1-2
#SBATCH -n 1

module purge
module load apps/gcc/R/4.5.0

Rscript p2_gen_data_devel.R 2 ${SLURM_ARRAY_TASK_ID} > p2_gen_data_devel_dgm2.${SLURM_ARRAY_TASK_ID}.out
Rscript p2_gen_data_devel.R 3 ${SLURM_ARRAY_TASK_ID} > p2_gen_data_devel_dgm3.${SLURM_ARRAY_TASK_ID}.out
Rscript p2_gen_data_devel.R 4 ${SLURM_ARRAY_TASK_ID} > p2_gen_data_devel_dgm4.${SLURM_ARRAY_TASK_ID}.out