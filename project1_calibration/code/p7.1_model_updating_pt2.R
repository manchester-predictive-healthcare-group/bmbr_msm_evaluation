### Hit an OOM issue on CSF. Re-running p7.1_model_updating.R file, from after where transition probabilities
### have been estimated.

### Setwd and clear workspace
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Extract arguments
args <- commandArgs(trailingOnly = T)
dgm <- as.numeric(args[1])
scenario <- as.numeric(args[2])
print(paste("dgm = ", dgm))
print(paste("scenario = ", scenario))

### Read in sim_inputs for the DGM
sim_inputs <- readRDS(paste("data/sim_inputs_dgm", dgm, ".rds", sep = ""))

### Set seed (set the seed using the scenario?)
set.seed(101)

### Read in devel_num and valid_num
devel_num <- as.numeric(sim_inputs[scenario, "devel_num"])
print(paste("devel num = ", devel_num))
valid_num <- as.numeric(sim_inputs[scenario, "valid_num"])
print(paste("valid num = ", valid_num))

### Read in validation data
df_valid_mstate <- readRDS(paste("data/dgm", dgm, "_df_valid_mstate", scenario, ".rds", sep = ""))
df_valid_raw <- readRDS(paste("data/dgm", dgm, "_df_valid_raw", scenario, ".rds", sep = ""))

### Save
tp <- readRDS(paste("data/dgm", dgm, "_pt_model_update_scenario", scenario, "_all.rds", sep = ""))

###
### Estimate calibration of transition probabilities
###

### Write a function to create the calibration object and plots and save them
estimate_calibration_tp <- function(t_in){
  
  ### A very small number of observations have negative predictions very close to zero
  ### These are when xa is simulated to be very large
  ### Set these to close to zero
  tp <- dplyr::mutate(
    tp,
    dplyr::across(
      .cols = dplyr::starts_with("pstate"),
      .fns  = ~ pmax(.x, 1e-6)
    )
  )
  
  ### Extract predictions at time of interest
  tp <- tp |> 
    dplyr::filter(t == t_in) |>
    dplyr::select(dplyr::starts_with("pstate"))
  
  ### Estimate calibration
  calib_object <- calibmsm::calib_msm(data_ms = df_valid_mstate,
                                      data_raw = df_valid_raw,
                                      j = 1,
                                      s = 0,
                                      t = t_in,
                                      tp_pred = tp,
                                      w_covs = c("xa", "xb"),
                                      calib_type = "blr",
                                      curve_type = "rcs",
                                      rcs_nk = 5,
                                      CI = FALSE,
                                      assess_moderate = TRUE,
                                      assess_mean = TRUE)
  
  ### Save object
  saveRDS(calib_object, paste0("data/calib_object_tp_model_update_dgm", dgm, "_s", scenario, "_t", t_in, ".rds"))
  
  ### Save plot
  ## First need to get number of states that can be reached
  ## Going to assign this manually
  
  # Get number of transitions
  if (dgm == 1){
    states_n <- 3
  } else {
    states_n <- 4
  }
  
  # Get number of rows (always assume three columns for consistency)
  plot_nrows <- ceiling(states_n/3)
  # Every row should have an extra inch, so can just use plot_nrows as the height argument in ragg::agg_png,
  # and the nrows argument for plot.calib_msm
  
  # Save plot
  ragg::agg_png(paste0("figures/calib_plot_tp_model_update_dgm", dgm, "_s", scenario, "_t", t_in, ".png"), 
                res = 600, 
                width = 3, 
                height = plot_nrows, 
                scaling = 1/3, 
                unit = "in")
  grid::grid.draw(plot(calib_object, nrow = plot_nrows, ncol = 3))
  dev.off()
  
  print(paste("FINISHED", scenario, "t = ", t_in, Sys.time()))
  
}
