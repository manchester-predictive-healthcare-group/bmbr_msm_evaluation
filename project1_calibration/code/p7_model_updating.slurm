#!/bin/bash --login
#SBATCH -p multicore
#SBATCH -t 24:00:00
#SBATCH -a 1-64
#SBATCH -n 3

module purge
module load apps/gcc/R/4.5.0

# --- Define how many scenarios for each DGM ---
# -a range should be equal to the sum of scenarios_per_dgm
scenarios_per_dgm=(40 6 8 10)

# --- Map linear task ID to dgm, scenario
task=${SLURM_ARRAY_TASK_ID}

dgm=0
scenario=0
cumulative=0

# --- Create input for each Rscript
for n_scenarios in "${scenarios_per_dgm[@]}"; do
    # Increment devel_num each iteration (started at 0, so becomes 1, 2, ...)
    dgm=$(( dgm + 1 ))
    # Keep a running total of tasks accounted for so far
    cumulative=$(( cumulative + n_scenarios ))
    # Check if our task ID falls within this devel_num's block
    if [ "$task" -le "$cumulative" ]; then
        # How far into this dgm's block are we?
        offset=$(( cumulative - n_scenarios ))
        scenario=$(( task - offset ))
        # Stop looping — we found our (dgm, scenario) pair
        break
    fi
done

Rscript p7.1_model_updating.R ${dgm} ${scenario} > p7.1_model_updating.dgm${dgm}.${scenario}.out