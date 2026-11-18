###
### Assess calibration of transition probabilities
###

### Setwd and clear workspace
rm(list = ls())
setwd("PATH/TO/ROOT")

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

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
  
  ### Enforce a small positive lower bound as a numerical safeguard
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
  ragg::agg_png(paste0("figures/calib_plot_tp_dgm", dgm, "_s", scenario, "_t", t_in, ".png"), 
                res = 600, 
                width = 3, 
                height = plot_nrows, 
                scaling = 1/3, 
                unit = "in")
  grid::grid.draw(plot(calib_object, nrow = plot_nrows, ncol = 3))
  dev.off()
  
  print(paste("FINISHED", scenario, "t = ", t_in, Sys.time()))
  
}

###
### Run function for DGM 2 (6 scenarios)
###
lapply(1:6, function(x){
  lapply(c(5), function(y){
    estimate_calibration_tp(dgm = 2, x, y)
  })
})

###
### Run function for DGM 3 (8 scenarios)
###
lapply(1:8, function(x){
  lapply(c(5), function(y){
    estimate_calibration_tp(dgm = 3, x, y)
  })
})

###
### Run function for DGM 4 (10 scenarios)
###
lapply(1:10, function(x){
  lapply(c(5), function(y){
    estimate_calibration_tp(dgm = 4, x, y)
  })
})

###
### Run function for DGM 1 (40 scenarios)
###
# c(sapply(seq(1,40, by = 8), function(s) s:(s+4)))
lapply(1:40, function(x){
  lapply(c(5), function(y){
    estimate_calibration_tp(dgm = 1, x, y)
  })
})

