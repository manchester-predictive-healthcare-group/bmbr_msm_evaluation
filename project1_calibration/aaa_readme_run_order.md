BMBR multistate calibration simulation: run order

Before running, replace PATH/TO/ROOT in the R and R Markdown files with the project directory containing the data, figures, and R subdirectories. If using the batch submission script, also replace PATH/TO/BATCH/WORKSPACE with the directory containing the R and Slurm files.

The complete Slurm workflow can be submitted using:

    ./aa_submit_calibration_pipeline.sh

The individual stages are listed below in execution order.

1. Define simulation inputs

   R file:
     p1_define_inputs.R

   Slurm file:
     p1_define_inputs.slurm

   Creates the transition matrices and simulation-input files for DGM1-DGM4.

2. Simulate development datasets

   R file:
     p2_gen_data_devel.R

   Slurm files:
     p2_gen_data_devel_dgm1.slurm   DGM1 development scenarios 1-8
     p2_gen_data_devel_dgm2.slurm   DGM2 development scenarios 1-2
     p2_gen_data_devel_dgm3.slurm   DGM3 development scenarios 1-2
     p2_gen_data_devel_dgm4.slurm   DGM4 development scenarios 1-2

3. Simulate validation datasets

   R file:
     p3_gen_data_valid.R

   Slurm files:
     p3_gen_data_valid_dgm1.slurm   DGM1 validation scenarios 1-40
     p3_gen_data_valid_dgm2.slurm   DGM2 validation scenarios 1-6
     p3_gen_data_valid_dgm3.slurm   DGM3 validation scenarios 1-8
     p3_gen_data_valid_dgm4.slurm   DGM4 validation scenarios 1-10

   The development and validation simulations can run in parallel after step 1.

4. Fit the transition-specific models

   R files:
     p4_fit_model_dgm1.R
     p4_fit_model_dgm2.R
     p4_fit_model_dgm3.R
     p4_fit_model_dgm4.R

   Slurm file:
     p4_fit_models.slurm

   Run after all development datasets have been generated.

5. Estimate transition probabilities

   R file:
     p5.1_pred_tp.R

   Slurm files:
     p5.1_pred_tp_dgm1.slurm   8 development models x 100 batches
     p5.1_pred_tp_dgm2.slurm   2 development models x 100 batches
     p5.1_pred_tp_dgm3.slurm   2 development models x 100 batches
     p5.1_pred_tp_dgm4.slurm   2 development models x 100 batches

   Run after model fitting and validation-data generation have completed.

6. Combine transition-probability batches

   R file:
     p5.2_pred_tp_combine.R

   Slurm file:
     p5.2_pred_tp_combine.slurm

   Run after all transition-probability arrays have completed.

7. Assess cause-specific calibration

   R file:
     p6.1_assess_calibration_csh.R

   Slurm file:
     p6.1_assess_calibration_csh.slurm

   Runs all 64 DGM/scenario combinations internally. It can run once model fitting and validation-data generation have completed.

8. Assess transition-probability calibration

   R file:
     p6.2_assess_calibration_tp_csf.R

   Slurm file:
     p6.2_assess_calibration_tp.slurm

   The 64-task array covers DGM1 scenarios 1-40, DGM2 scenarios 1-6, DGM3 scenarios 1-8, and DGM4 scenarios 1-10. Run after the transition-probability batches have been combined.

9. Create combined calibration figures for the original models

   R file:
     p6.3_create_combined_calib_figures.R

   Slurm file:
     p6.3_create_combined_calib_figures.slurm

   Run after both cause-specific and transition-probability calibration have completed.

10. Update selected transition models and reassess calibration

   R file:
     p7.1_model_updating.R

   Slurm file:
     p7.1_model_updating.slurm

   The 64-task array covers DGM1 scenarios 1-40, DGM2 scenarios 1-6, DGM3 scenarios 1-8, and DGM4 scenarios 1-10. Run after cause-specific calibration has completed.

11. Create combined calibration figures for the updated models

   R file:
     p7.2_create_combined_calib_figures_model_updating.R

   Slurm file:
     p7.2_create_combined_calib_figures_model_updating.slurm

   Run after model updating and the original-model figure stage have completed.

12. Create the supplementary results document

   R Markdown file:
     p8_supp_data_file_1.Rmd

   Render after the original and updated combined figures have been created.

13. Optional EBMT example

   R file:
     p9_apply_functions_ebmt.R
