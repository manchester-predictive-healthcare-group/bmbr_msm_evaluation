#!/bin/bash --login

set -euo pipefail

### Move to the workspace containing all R and Slurm files
cd "PATH/TO/BATCH/WORKSPACE"

### Define simulation inputs
job1=$(sbatch --parsable p1_define_inputs.slurm)
echo "Submitted input-definition job: ${job1}"

### Simulate development datasets in parallel after inputs are defined
job2a=$(sbatch --parsable --dependency="afterok:${job1}" p2_gen_data_devel_dgm1.slurm)
job2b=$(sbatch --parsable --dependency="afterok:${job1}" p2_gen_data_devel_dgm2.slurm)
job2c=$(sbatch --parsable --dependency="afterok:${job1}" p2_gen_data_devel_dgm3.slurm)
job2d=$(sbatch --parsable --dependency="afterok:${job1}" p2_gen_data_devel_dgm4.slurm)
echo "Submitted development arrays: ${job2a}, ${job2b}, ${job2c}, ${job2d}"

### Simulate validation datasets in parallel after inputs are defined
job3a=$(sbatch --parsable --dependency="afterok:${job1}" p3_gen_data_valid_dgm1.slurm)
job3b=$(sbatch --parsable --dependency="afterok:${job1}" p3_gen_data_valid_dgm2.slurm)
job3c=$(sbatch --parsable --dependency="afterok:${job1}" p3_gen_data_valid_dgm3.slurm)
job3d=$(sbatch --parsable --dependency="afterok:${job1}" p3_gen_data_valid_dgm4.slurm)
echo "Submitted validation arrays: ${job3a}, ${job3b}, ${job3c}, ${job3d}"

### Fit all four DGMs after every development-data array succeeds
job4=$(sbatch --parsable --dependency="afterok:${job2a}:${job2b}:${job2c}:${job2d}" p4_fit_models.slurm)
echo "Submitted model-fitting array: ${job4}"

### Predict transition probabilities after models and validation data are available
prediction_dependency="afterok:${job4}:${job3a}:${job3b}:${job3c}:${job3d}"
job5a=$(sbatch --parsable --dependency="${prediction_dependency}" p5.1_pred_tp_dgm1.slurm)
job5b=$(sbatch --parsable --dependency="${prediction_dependency}" p5.1_pred_tp_dgm2.slurm)
job5c=$(sbatch --parsable --dependency="${prediction_dependency}" p5.1_pred_tp_dgm3.slurm)
job5d=$(sbatch --parsable --dependency="${prediction_dependency}" p5.1_pred_tp_dgm4.slurm)
echo "Submitted prediction arrays: ${job5a}, ${job5b}, ${job5c}, ${job5d}"

### Combine prediction batches only after every prediction array succeeds
job6=$(sbatch --parsable --dependency="afterok:${job5a}:${job5b}:${job5c}:${job5d}" p5.2_pred_tp_combine.slurm)
echo "Submitted prediction-combination job: ${job6}"

### Assess cause-specific calibration after models and validation data are available
job7=$(sbatch --parsable --dependency="afterok:${job4}:${job3a}:${job3b}:${job3c}:${job3d}" p6.1_assess_calibration_csh.slurm)
echo "Submitted cause-specific calibration job: ${job7}"

### Assess transition-probability calibration after combined predictions are available
job8=$(sbatch --parsable --dependency="afterok:${job6}" p6.2_assess_calibration_tp.slurm)
echo "Submitted transition-probability calibration array: ${job8}"

### Create original-model combined figures after both calibration stages succeed
job9=$(sbatch --parsable --dependency="afterok:${job7}:${job8}" p6.3_create_combined_calib_figures.slurm)
echo "Submitted original-model figure job: ${job9}"

### Run model updating after cause-specific calibration is available
job10=$(sbatch --parsable --dependency="afterok:${job7}" p7.1_model_updating.slurm)
echo "Submitted model-updating array: ${job10}"

### Create updated-model figures after updating and the original figure stage succeed
job11=$(sbatch --parsable --dependency="afterok:${job9}:${job10}" p7.2_create_combined_calib_figures_model_updating.slurm)
echo "Submitted updated-model figure job: ${job11}"

echo "Pipeline submitted"
echo "Inputs: ${job1}"
echo "Development: ${job2a} ${job2b} ${job2c} ${job2d}"
echo "Validation: ${job3a} ${job3b} ${job3c} ${job3d}"
echo "Models: ${job4}"
echo "Predictions: ${job5a} ${job5b} ${job5c} ${job5d}"
echo "Combined predictions: ${job6}"
echo "CSH calibration: ${job7}"
echo "TP calibration: ${job8}"
echo "Original figures: ${job9}"
echo "Model updating: ${job10}"
echo "Updated figures: ${job11}"
