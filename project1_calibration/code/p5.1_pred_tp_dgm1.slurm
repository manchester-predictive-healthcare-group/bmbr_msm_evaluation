#!/bin/bash --login
#SBATCH -p serial
#SBATCH -t 01:00:00
#SBATCH -a 1-800
#SBATCH -n 1

module purge
module load apps/gcc/R/4.5.0

devel_num=$(( (SLURM_ARRAY_TASK_ID - 1) / 100 + 1 )) # gives 1-8 
set=$(( (SLURM_ARRAY_TASK_ID - 1) % 100 + 1 )) # gives 1-100

Rscript p5.1_pred_tp.R 1 ${devel_num} 5000 ${set} > p5.1_pred_tp_dgm1.${devel_num}.${set}.out