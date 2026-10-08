### Let's write a program will that do model updating

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

### Set t_in == 5
### We are identifying which model transitions need to be updated based on the cause-specific hazard plots
### at this time point.
t_in <- 5

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

### Read in the model
msm_flexsurv_fit <- readRDS(paste("data/dgm", dgm, "_msm_flexsurv_fit", devel_num, ".rds", sep = ""))

### Read in CSH calibration data to identify which transitions need to be updated
calib_object_csh <- readRDS(paste0("data/calib_object_csh_dgm",dgm,"_s",scenario,"_t",t_in,".rds"))

### Get vector for transitions where ICI > 0.1
### These numbers correspond to the elements of the model fit object
### Note an ICI > 0.1 is not a good choice for whether to update or not. This is just a threshold above/below
### which we know the ICI of the cause-specific hazard models lie depending on whether there was a shift or not.
### This choice is threfore simply to implement the analyes in this simulation.
transitions_update <- which(sapply(calib_object_csh, function(x) {
  x$ICI > 0.1
}))

### Update relevant models (note, this could update multiple transitions, but for most simulation scenarios,
### it's just updating one).
### We only update scale and shape (or just rate for exponential), leaving coefficients untouched, 
### but these could also be updated by altering the fixedpars argument.
if (length(transitions_update) > 0){
  for (transition in transitions_update){
    
    ### Create initial vector with terms on correct scale
    inits_vec <- coef(msm_flexsurv_fit[[transition]])
    
    ### DGM-1 development models 1-4 were fitted as exponential
    if (dgm == 1 && devel_num %in% 1:4){
      
      ### Put the exponential rate onto its natural scale
      inits_vec[1] <- exp(inits_vec[1])
      
      ### Re-estimate the rate while fixing xa and xb
      msm_flexsurv_fit[[transition]] <- flexsurv::flexsurvreg(
        survival::Surv(time, status) ~ xa + xb,
        subset = (trans == transition),
        data = df_valid_mstate,
        dist = "exp",
        inits = inits_vec,
        fixedpars = c(2,3)
      )
      ### All other development models fitted as Weibull
    } else {
      ### Put shape and scale onto natural scales
      inits_vec[c(1,2)] <- exp(inits_vec[c(1,2)])
      
      ## Update model
      msm_flexsurv_fit[[transition]] <- flexsurv::flexsurvreg(
        survival::Surv(time, status) ~ xa + xb,
        subset = (trans == transition),
        data = df_valid_mstate,
        dist = "weibullPH",
        inits = inits_vec,   # start from original fit
        fixedpars = c(3, 4)                     # fix xa, xb;
      )}
    
  }
}

### Save updated model
saveRDS(msm_flexsurv_fit, paste("data/dgm", dgm, "_msm_updated_flexsurv_fit", scenario, ".rds", sep = ""))
print(paste("Model updated and saved", Sys.time()))

###
### Estimate calibration CSH
###

### Write a function to create the calibration object and plots and save them
estimate_calibration_csh <- function(t_in){
  
  print(paste("Calibrate CSH", Sys.time()))
  
  ### Assess calibration
  calib_object <- calib_csh_flexsurv(data_mstate = df_valid_mstate, fit_msm = msm_flexsurv_fit, t = t_in, nk = 5)
  
  ### Save plot
  ## First need to get number of transitions, so I know number of rows in plot, so I can save correct
  ## height (always same number of columns)
  
  # Get number of transitions
  transitions_n <- max(df_valid_mstate$trans)
  # Get number of columns
  plot_ncols <- 3
  # Get number of rows 
  plot_nrows <- ceiling(transitions_n/plot_ncols)
  # Every row should have an extra inch, so can just use plot_nrows as the height argument
  
  # Save plot
  ragg::agg_png(paste0("figures/calib_plot_csh_model_update_dgm", dgm, "_s", scenario, "_t", t_in, ".png"), 
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
  saveRDS(calib_object, paste0("data/calib_object_csh_model_update_dgm", dgm, "_s", scenario, "_t", t_in, ".rds"))
  
  print(paste("FINISHED csh calibration ", scenario, "t = ", t_in, Sys.time()))
  warnings()
}

### Run this function
estimate_calibration_csh(5)
print(paste("Calibration of cause-specific hazard assessed", Sys.time()))

###
### Estimate transition probabilities in the validation dataset using the updated model
###

### Define tmat
tmat <- attributes(df_valid_mstate)$trans

### Calculate transition probabilities
tp <- lapply(FUN = est_risk_flexsurv, 
             1:nrow(df_valid_raw), 
             data_raw = df_valid_raw, 
             tmat = tmat, 
             msm_fit = msm_flexsurv_fit,
             t_vec = c(2.5, 5, 7.5),
             m_iter = 10000)

## Turn into data.frame
tp <- 
  do.call(rbind, lapply(tp, function(x){
    do.call(rbind, x)
  }))

### If any values are equal to 0 and 1, change
tp <- dplyr::mutate(tp, dplyr::across(
  .cols = dplyr::starts_with("pstate"),
  .fns = ~ pmin(pmax(.x, 0.000001), 0.999999)
))


### Save
saveRDS(tp, paste("data/dgm", dgm, "_pt_model_update_scenario", scenario, "_all.rds", sep = ""))
print(paste("Transition probabilities estimated", Sys.time()))

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

### Estimate calibration of transition probabilities
estimate_calibration_tp(5)
print(paste("Calibration of transition probabilities assessed", Sys.time()))