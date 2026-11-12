###
### Fit multistate models
###

### Setwd and clear workspace
rm(list = ls())
setwd("PATH/TO/PROJECT")

### Define dgm
dgm <- 2

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Loop through development datasets
for (devel_num in 1:2){
  
  print(paste("devel_num = ", devel_num, Sys.time()))
  
  ### Read in development data
  df_devel_mstate <- readRDS(paste("data/dgm", dgm, "_df_devel_mstate", devel_num, ".rds", sep = ""))
  df_devel_raw <- readRDS(paste("data/dgm", dgm, "_df_devel_raw", devel_num, ".rds", sep = ""))

  ### Define tmat
  tmat <- attributes(df_devel_mstate)$trans
  
  ### Select the model distribution for this development scenario
  model_dist <- "weibullPH"
  
  ### Identify the from-state and to-state for each numbered transition
  transition_positions <- which(tmat > 0, arr.ind = TRUE)
  
  ### Order the state pairs by their transition numbers
  transition_positions <- transition_positions[order(tmat[transition_positions]),,drop = FALSE]
  
  ### Construct the corresponding unique predictor names
  ### For DGM1, this produces x12, x13 and x23
  transition_predictors <- paste0(
    "x",
    transition_positions[, "from"],
    transition_positions[, "to"]
  )
  
  ### Confirm that the matrix positions are aligned with transition numbers
  stopifnot(
    identical(
      as.integer(tmat[transition_positions]),
      seq_along(transition_predictors)
    )
  )
  
  ### Create an empty list for the fitted transition-specific models
  msm_flexsurv_fit <- vector(
    "list",
    length(transition_predictors)
  )
  
  ### Fit each transition using its corresponding unique predictor
  for (transition in seq_along(transition_predictors)){
    
    ### Select the predictor associated with this transition
    unique_predictor <- transition_predictors[transition]
    
    ### Print the transition, predictor and start time
    message("Fitting development transition ",transition," using ",unique_predictor,"; ",Sys.time())
    
    ### Construct the transition-specific model formula
    model_formula <- stats::reformulate(
      c("xa", "xb", unique_predictor),
      response = "survival::Surv(time, status)"
    )
    
    ### Fit the transition-specific model
    msm_flexsurv_fit[[transition]] <- flexsurv::flexsurvreg(
      model_formula,
      subset = trans == transition,
      data = df_devel_mstate,
      dist = model_dist
    )
  }
  
  ### Save
  saveRDS(msm_flexsurv_fit, paste("data/dgm", dgm, "_msm_flexsurv_fit", devel_num, ".rds", sep = ""))
  print("FINISHED")
  
}



