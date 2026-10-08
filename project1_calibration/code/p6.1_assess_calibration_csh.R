###
### Assess calibration of cause-specific hazards
###

### Setwd and clear workspace
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Write a function to create the calibration object and plots and save them
estimate_calibration_csh <- function(dgm, scenario, t_in){
  
  ### Read in sim_inputs for the DGM
  sim_inputs <- readRDS(paste("data/sim_inputs_dgm", dgm, ".rds", sep = ""))
  
  ### Paste scenario
  print(paste("scenario = ", scenario))
  
  ### For each validation dataset, there is a corresponding development dataset, stored in the variable devel_num
  ### This defines which model is read in and used to make predictions
  devel_num <- as.numeric(sim_inputs[scenario, "devel_num"])
  print(paste("devel num = ", devel_num))
  
  ### Also extract the validation dataset type number for printing purposes only
  valid_num <- as.numeric(sim_inputs[scenario, "valid_num"])
  print(paste("valid num = ", valid_num))
  
  ### Read in model fit
  msm_flexsurv_fit <- readRDS(paste("data/dgm", dgm, "_msm_flexsurv_fit", devel_num, ".rds", sep = ""))
  
  ### Load the validation dataset
  df_valid_mstate <- readRDS(paste("data/dgm", dgm, "_df_valid_mstate", scenario, ".rds", sep = ""))
  
  ### Assess calibration
  calib_object <- calib_csh_flexsurv(data_mstate = df_valid_mstate, fit_msm = msm_flexsurv_fit, t = t_in, nk = 5)
  
  ### Save plot
  ## First need to get number of transitions, so I know number of rows in plot, so I can save correct
  ## height (always same number of columns)
  
  # Get number of transitions
  transitions_n <- max(df_valid_mstate$trans)
  # Get number of columns
  plot_ncols <- ifelse(transitions_n <= 3, 3, 4)
  # Get number of rows 
  plot_nrows <- ceiling(transitions_n/plot_ncols)
  # Every row should have an extra inch, so can just use plot_nrows as the height argument
  
  # Save plot
  ragg::agg_png(paste0("figures/calib_plot_csh_dgm", dgm, "_s", scenario, "_t", t_in, ".png"), 
                res = 600, 
                width = plot_ncols, 
                height = plot_nrows, 
                scaling = 1/plot_ncols, 
                unit = "in")
  plot(calib_object[["plots_comb"]])
  dev.off()
  
  ### NB: may need to change the following few lines in the future, depending if we make any changes to calib_csh, or est_calib_csh_ph
  ### Remove the combined plot before saving the output object
  calib_object<- calib_object[["calib_list"]]
  
  ### Reduce the plot list component to just plot data, so plots can be reproduced, but don't want to save all the ggplot stuff which 
  ### takes up space on disk
  calib_object <- lapply(calib_object, function(x){
    x[["plot"]] <- x[["plot"]][["data"]]
    names(x)[names(x) == "plot"] <- "plotdata"
    x[["calib_data"]] <- NULL
    return(x)
  })
  
  ### Save object
  saveRDS(calib_object, paste0("data/calib_object_csh_dgm", dgm, "_s", scenario, "_t", t_in, ".rds"))
  
  print(paste("FINISHED", scenario, "t = ", t_in, Sys.time()))
  warnings()
}

###
### Run function for DGM 1 (40 scenarios)
###
# c(sapply(seq(1,40, by = 8), function(s) s:(s+4)))
lapply(1:40, function(x){
  lapply(c(5), function(y){
    estimate_calibration_csh(dgm = 1, x, y)
  })
})

###
### Run function for DGM 2 (6 scenarios)
###
lapply(1:6, function(x){
  lapply(c(5), function(y){
    estimate_calibration_csh(dgm = 2, x, y)
  })
})

###
### Run function for DGM 3 (8 scenarios)
###
lapply(1:8, function(x){
  lapply(c(5), function(y){
    estimate_calibration_csh(dgm = 3, x, y)
  })
})

###
### Run function for DGM 4 (10 scenarios)
###
lapply(1:10, function(x){
  lapply(c(5), function(y){
    estimate_calibration_csh(dgm = 4, x, y)
  })
})
