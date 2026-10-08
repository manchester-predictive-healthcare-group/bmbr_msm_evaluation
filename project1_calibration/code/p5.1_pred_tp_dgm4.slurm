#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 01:00:00
#SBATCH -a 1-200
#SBATCH -n 1

module purge
module load apps/gcc/R/4.5.0

devel_num=$(( (SLURM_ARRAY_TASK_ID - 1) / 100 + 1 )) # gives 1-2
set=$(( (SLURM_ARRAY_TASK_ID - 1) % 100 + 1 )) # gives 1-100

Rscript p5.1_pred_tp.R 4 ${devel_num} 5000 ${set} > p5.1_pred_tp_dgm4.${devel_num}.${set}.out