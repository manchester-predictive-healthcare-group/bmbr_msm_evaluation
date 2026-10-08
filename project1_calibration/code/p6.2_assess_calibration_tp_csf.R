###
### Assess calibration of transition probabilities
### Arguments are read in at command line so it can be parallelised on CSF
###

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

### Write a function to create the calibration object and plots and save them
estimate_calibration_tp <- function(dgm, scenario, t_in){
  
  print(paste("scenario = ", scenario))
  
  ### Read in sim_inputs for the DGM
  sim_inputs <- readRDS(paste("data/sim_inputs_dgm", dgm, ".rds", sep = ""))
  
  ### For each validation dataset, there is a corresponding development dataset, stored in the variable devel_num
  ### This defines which model is read in and used to make predictions
  devel_num <- as.numeric(sim_inputs[scenario, "devel_num"])
  print(paste("devel num = ", devel_num))
  
  ### Also extract the validation dataset type number for printing purposes only
  valid_num <- as.numeric(sim_inputs[scenario, "valid_num"])
  print(paste("valid num = ", valid_num))
  
  ### Load transition probabilities
  tp <- readRDS(paste("data/dgm", dgm, "_pt_scenario", devel_num, "_all.rds", sep = ""))
  
  ### A very small number of observations have negative predictions very close to zero (e.g. -7*10^312)
  ### These are when xa is simulated to be very large
  ### Set these to close to zero
  tp <- dplyr::mutate(
    tp,
    dplyr::across(
      .cols = dplyr::starts_with("pstate"),
      .fns  = ~ pmax(.x, 1e-6)
    )
  )
  
  ### Load the validation datasets
  df_valid_mstate <- readRDS(paste("data/dgm", dgm, "_df_valid_mstate", scenario, ".rds", sep = ""))
  df_valid_raw <- readRDS(paste("data/dgm", dgm, "_df_valid_raw", scenario, ".rds", sep = ""))
  
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
  saveRDS(calib_object, paste0("data/calib_object_tp_dgm", dgm, "_s", scenario, "_t", t_in, ".rds"))
  
  ### Save plot
  ## First need to get number of states that can be reached
  ## Going to assign this manually
  
  # Get number of states than can be transitioned into
  if (dgm == 1){
    states_n <- 3
  } else {
    states_n <- 4
  }

  # Get number of columns
  plot_ncols <- 3
  # Get number of rows 
  plot_nrows <- ceiling(states_n/plot_ncols)
  # Every row should have an extra inch, so can just use plot_nrows as the height argument
  
  # Every row should have an extra inch, so can just use plot_nrows as the height argument in ragg::agg_png,
  # and the nrows argument for plot.calib_msm
  
  # Save plot
  ragg::agg_png(paste0("figures/calib_plot_tp_dgm", dgm, "_s", scenario, "_t", t_in, ".png"), 
                res = 600, 
                width = plot_ncols, 
                height = plot_nrows, 
                scaling = 1/plot_ncols, 
                unit = "in")
  grid::grid.draw(plot(calib_object, nrow = plot_nrows, ncol = plot_ncols))
  dev.off()
  
  print(paste("FINISHED", scenario, "t = ", t_in, Sys.time()))
  
}

### Run function
lapply(c(5), function(x){estimate_calibration_tp(dgm = dgm, scenario = scenario, t_in = x)})
