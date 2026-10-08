
###
### DGM1: Simulate transition data
###
dgm1 <- function(n, #number of patients to simulate
                 max_follow, #maximum follow up
                 shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                 shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                 shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                 beta12_xa, beta12_xb, beta12_xasq, beta12_xbsq, #covariate effects for transition 12
                 beta13_xa, beta13_xb, beta13_xasq, beta13_xbsq, #covariate effects for transition 13
                 beta23_xa, beta23_xb, beta23_xasq, beta23_xbsq, #covariate effects for transition 23
                 x_in, #baseline predictors, dataframe with two columns (xa continuous, xb binary)
                 numsteps, ...) #number of sampler steps in gems data generation process
{
  
  #   n.cohort <- 100
  #   x.baseline <- data.frame("xa" = rnorm(n.cohort, 0, 1), "xb" = rnorm(n.cohort, 0, 1))
  #   n <- n.cohort
  #   max_follow <- ceiling(365.25*7)
  #
  #   ### Baseline hazards
  #   shape12 <- 1
  #   scale12 <- 1588.598
  #
  #   shape13 <- 1
  #   scale13 <- 0.5*1588.598
  #
  #   shape23 <- 1
  #   scale23 <- 5*1588.598
  #
  #   #qweibull(0.8, 1, 1588.598)
  #
  #   ## Covariate effects
  #   beta12_xa <- 1
  #   beta12_xb <- 1
  #   beta13_xa <- 0.5
  #   beta13_xb <- 0.5
  #   beta23_xa <- 1
  #   beta23_xb <- 0.5
  #
  #   x_in <- x.baseline
  #   numsteps <- max_follow
  
  ## Generate a baseline covariate data frame
  bl <- x_in
  
  ## Generate an empty hazard matrix
  hf <- gems::generateHazardMatrix(3)
  #hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[1, 3]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[2, 3]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  print(hf)
  
  
  ## We are using a clock reset approach to generate data
  ## Replace t with (t + sum(history)) to implement a clock forward approach
  
  ## Generate an empty parameter matrix
  par <- gems::generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12,
                      beta_xa = beta12_xa, beta_xb = beta12_xb, beta_xasq = beta12_xasq, beta_xbsq = beta12_xbsq)
  par[[1, 3]] <- list(shape = shape13, scale = scale13,
                      beta_xa = beta13_xa, beta_xb = beta13_xb, beta_xasq = beta13_xasq, beta_xbsq = beta13_xbsq)
  par[[2, 3]] <- list(shape = shape23, scale = scale23,
                      beta_xa = beta23_xa, beta_xb = beta23_xb, beta_xasq = beta23_xasq, beta_xbsq = beta23_xbsq)
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- gems::simulateCohort(transitionFunctions = hf, parameters = par,
                                 cohortSize = n, baseline = bl, to = max_follow, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in
  
  ## Rename the column names
  colnames(cohort@time.to.state) <- paste0("state", 1:ncol(cohort@time.to.state))
  
  ## Get data into the common data model
  cohort_out <- data.frame(cohort@time.to.state, as.data.frame(cohort@baseline), patid = 1:nrow(cohort@time.to.state))
  
  return(cohort_out)
  
}


###
### DGM2: Simulate transition data
###
dgm2 <- function(n, #number of patients to simulate
                 max_follow, #maximum follow up
                 shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                 shape14, scale14, #shape and scale for weibull baseline hazard for transition 1 -> 4
                 shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                 beta12_xa, beta12_xb, beta12_xasq, beta12_xbsq, #covariate effects for transition 1 -> 2
                 beta14_xa, beta14_xb, beta14_xasq, beta14_xbsq, #covariate effects for transition 1 -> 4
                 beta23_xa, beta23_xb, beta23_xasq, beta23_xbsq, #covariate effects for transition 2 -> 3
                 x_in, #baseline predictors, dataframe with two columns (xa continuous, xb binary)
                 numsteps, ...) #number of sampler steps in gems data generation process
{
  
  #   n.cohort <- 100
  #   x.baseline <- data.frame("xa" = rnorm(n.cohort, 0, 1), "xb" = rnorm(n.cohort, 0, 1))
  #   n <- n.cohort
  #   max_follow <- ceiling(365.25*7)
  #
  #   ### Baseline hazards
  #   shape12 <- 1
  #   scale12 <- 1588.598
  #
  #   shape14 <- 1
  #   scale14 <- 0.5*1588.598
  #
  #   shape23 <- 1
  #   scale23 <- 5*1588.598
  #
  #   #qweibull(0.8, 1, 1588.598)
  #
  #   ## Covariate effects
  #   beta12_xa <- 1
  #   beta12_xb <- 1
  #   beta14_xa <- 0.5
  #   beta14_xb <- 0.5
  #   beta23_xa <- 1
  #   beta23_xb <- 0.5
  #
  #   x_in <- x.baseline
  #   numsteps <- max_follow
  
  ## Generate a baseline covariate data frame
  bl <- x_in
  
  ## Generate an empty hazard matrix
  hf <- gems::generateHazardMatrix(4)
  #hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[1, 4]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[2, 3]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  print(hf)
  
  
  ## We are using a clock reset approach to generate data
  ## Replace t with (t + sum(history)) to implement a clock forward approach
  
  ## Generate an empty parameter matrix
  par <- gems::generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12,
                      beta_xa = beta12_xa, beta_xb = beta12_xb, beta_xasq = beta12_xasq, beta_xbsq = beta12_xbsq)
  par[[1, 4]] <- list(shape = shape14, scale = scale14,
                      beta_xa = beta14_xa, beta_xb = beta14_xb, beta_xasq = beta14_xasq, beta_xbsq = beta14_xbsq)
  par[[2, 3]] <- list(shape = shape23, scale = scale23,
                      beta_xa = beta23_xa, beta_xb = beta23_xb, beta_xasq = beta23_xasq, beta_xbsq = beta23_xbsq)
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- gems::simulateCohort(transitionFunctions = hf, parameters = par,
                                 cohortSize = n, baseline = bl, to = max_follow, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in
  
  ## Rename the column names
  colnames(cohort@time.to.state) <- paste0("state", 1:ncol(cohort@time.to.state))
  
  ## Get data into the common data model
  cohort_out <- data.frame(cohort@time.to.state, as.data.frame(cohort@baseline), patid = 1:nrow(cohort@time.to.state))
  
  return(cohort_out)
  
}

###
### DGM3: Simulate transition data
###
dgm3 <- function(n, #number of patients to simulate
                 max_follow, #maximum follow up
                 shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                 shape14, scale14, #shape and scale for weibull baseline hazard for transition 1 -> 4
                 shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                 shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                 beta12_xa, beta12_xb, beta12_xasq, beta12_xbsq, #covariate effects for transition 1 -> 2
                 beta14_xa, beta14_xb, beta14_xasq, beta14_xbsq, #covariate effects for transition 1 -> 4
                 beta23_xa, beta23_xb, beta23_xasq, beta23_xbsq, #covariate effects for transition 2 -> 3
                 beta24_xa, beta24_xb, beta24_xasq, beta24_xbsq, #covariate effects for transition 2 -> 4
                 x_in, #baseline predictors, dataframe with two columns (xa continuous, xb binary)
                 numsteps, ...) #number of sampler steps in gems data generation process
{
  
  ## Generate a baseline covariate data frame
  bl <- x_in
  
  ## Generate an empty hazard matrix
  hf <- gems::generateHazardMatrix(4)
  #hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[1, 4]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[2, 3]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[2, 4]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  print(hf)
  
  
  ## We are using a clock reset approach to generate data
  ## Replace t with (t + sum(history)) to implement a clock forward approach
  
  ## Generate an empty parameter matrix
  par <- gems::generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12,
                      beta_xa = beta12_xa, beta_xb = beta12_xb, beta_xasq = beta12_xasq, beta_xbsq = beta12_xbsq)
  par[[1, 4]] <- list(shape = shape14, scale = scale14,
                      beta_xa = beta14_xa, beta_xb = beta14_xb, beta_xasq = beta14_xasq, beta_xbsq = beta14_xbsq)
  par[[2, 3]] <- list(shape = shape23, scale = scale23,
                      beta_xa = beta23_xa, beta_xb = beta23_xb, beta_xasq = beta23_xasq, beta_xbsq = beta23_xbsq)
  par[[2, 4]] <- list(shape = shape24, scale = scale24,
                      beta_xa = beta24_xa, beta_xb = beta24_xb, beta_xasq = beta24_xasq, beta_xbsq = beta24_xbsq)
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- gems::simulateCohort(transitionFunctions = hf, parameters = par,
                                 cohortSize = n, baseline = bl, to = max_follow, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in
  
  ## Rename the column names
  colnames(cohort@time.to.state) <- paste0("state", 1:ncol(cohort@time.to.state))
  
  ## Get data into the common data model
  cohort_out <- data.frame(cohort@time.to.state, as.data.frame(cohort@baseline), patid = 1:nrow(cohort@time.to.state))
  
  return(cohort_out)
  
}

###
### DGM4: Simulate transition data
###
dgm4 <- function(n, #number of patients to simulate
                 max_follow, #maximum follow up
                 shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                 shape14, scale14, #shape and scale for weibull baseline hazard for transition 1 -> 4
                 shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                 shape24, scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                 shape34, scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                 beta12_xa, beta12_xb, beta12_xasq, beta12_xbsq, #covariate effects for transition 1 -> 2
                 beta14_xa, beta14_xb, beta14_xasq, beta14_xbsq, #covariate effects for transition 1 -> 4
                 beta23_xa, beta23_xb, beta23_xasq, beta23_xbsq, #covariate effects for transition 2 -> 3
                 beta24_xa, beta24_xb, beta24_xasq, beta24_xbsq, #covariate effects for transition 2 -> 4
                 beta34_xa, beta34_xb, beta34_xasq, beta34_xbsq, #covariate effects for transition 3 -> 4
                 x_in, #baseline predictors, dataframe with two columns (xa continuous, xb binary)
                 numsteps, ...) #number of sampler steps in gems data generation process
{
  
  #   n.cohort <- 100
  #   x.baseline <- data.frame("xa" = rnorm(n.cohort, 0, 1), "xb" = rnorm(n.cohort, 0, 1))
  #   n <- n.cohort
  #   max_follow <- ceiling(365.25*7)
  #
  #   ### Baseline hazards
  #   shape12 <- 1
  #   scale12 <- 1588.598
  #
  #   shape14 <- 1
  #   scale14 <- 0.5*1588.598
  #
  #   shape23 <- 1
  #   scale23 <- 5*1588.598
  #
  #   #qweibull(0.8, 1, 1588.598)
  #
  #   ## Covariate effects
  #   beta12_xa <- 1
  #   beta12_xb <- 1
  #   beta14_xa <- 0.5
  #   beta14_xb <- 0.5
  #   beta23_xa <- 1
  #   beta23_xb <- 0.5
  #
  #   x_in <- x.baseline
  #   numsteps <- max_follow
  
  ## Generate a baseline covariate data frame
  bl <- x_in
  
  ## Generate an empty hazard matrix
  hf <- gems::generateHazardMatrix(4)
  #hf
  
  ## Change the entries of the transitions we want to allow
  ## Define the transitions as weibull
  hf[[1, 2]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[1, 4]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[2, 3]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[2, 4]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  hf[[3, 4]] <- function(t, shape, scale, beta_xa, beta_xb, beta_xasq, beta_xbsq) {
    exp(bl["xa"]*beta_xa + bl["xasq"]*beta_xasq + bl["xb"]*beta_xb + bl["xbsq"]*beta_xbsq)*(shape/scale)*(t/scale)^(shape - 1)}
  
  print(hf)
  
  
  ## We are using a clock reset approach to generate data
  ## Replace t with (t + sum(history)) to implement a clock forward approach
  
  ## Generate an empty parameter matrix
  par <- gems::generateParameterMatrix(hf)
  
  ## Use the vector of scales in each transition hazard
  par[[1, 2]] <- list(shape = shape12, scale = scale12,
                      beta_xa = beta12_xa, beta_xb = beta12_xb, beta_xasq = beta12_xasq, beta_xbsq = beta12_xbsq)
  par[[1, 4]] <- list(shape = shape14, scale = scale14,
                      beta_xa = beta14_xa, beta_xb = beta14_xb, beta_xasq = beta14_xasq, beta_xbsq = beta14_xbsq)
  par[[2, 3]] <- list(shape = shape23, scale = scale23,
                      beta_xa = beta23_xa, beta_xb = beta23_xb, beta_xasq = beta23_xasq, beta_xbsq = beta23_xbsq)
  par[[2, 4]] <- list(shape = shape24, scale = scale24,
                      beta_xa = beta24_xa, beta_xb = beta24_xb, beta_xasq = beta24_xasq, beta_xbsq = beta24_xbsq)
  par[[3, 4]] <- list(shape = shape34, scale = scale34,
                      beta_xa = beta34_xa, beta_xb = beta34_xb, beta_xasq = beta34_xasq, beta_xbsq = beta34_xbsq)
  
  ## Generate the cohort
  time.in <- Sys.time()
  cohort <- gems::simulateCohort(transitionFunctions = hf, parameters = par,
                                 cohortSize = n, baseline = bl, to = max_follow, sampler.steps = numsteps)
  time.out <- Sys.time()
  time.diff <- time.out - time.in
  
  ## Rename the column names
  colnames(cohort@time.to.state) <- paste0("state", 1:ncol(cohort@time.to.state))
  
  ## Get data into the common data model
  cohort_out <- data.frame(cohort@time.to.state, as.data.frame(cohort@baseline), patid = 1:nrow(cohort@time.to.state))
  
  return(cohort_out)
  
}

###
### Convert data to mstate format and applying censoring, with the option to round event times to help with computation time
### when fitting models with mstate. When fitting models with flexsurv, rounding not neccesary.
###
convert_mstate_cens <- function(cohort_in,
                                max_follow,
                                cens_shape,
                                cens_scale,
                                cens_beta_xa,
                                cens_beta_xb,
                                dgm_in,
                                rounding = FALSE,
                                dp = 3){
  
  #   cohort_in <- cohort[["cohort"]]
  #   max_follow <- cohort[["max_follow"]]
  #   cens_shape <- 1
  #   cens_scale <- 4000
  #   cens_beta_xa <- 0
  #   cens_beta_xb <- 0
  
  # cohort_in = df_devel
  # max_follow = max_follow
  # cens_shape = cens_shape
  # cens_scale = cens_scale
  # cens_beta_xa = cens_beta_xa
  # cens_beta_xb = cens_beta_xb
  # dgm_in = 2
  
  ### NBNBNB READ
  ### Function to round event times to nearest dp and create a gap of at least (10^dp) between all events
  ### NB: If applying rounding, want to test this function, but I don't foresee it will be used, 
  ### given this set is only to help with computation time, and planning to use flexsurv which has minimal computational issues.
  ### The function to create gaps may not be a perfect solution, but this happens extremely rarely. I think issues will arise 
  ### with this code in models with multiple absorbing states
  adjust_event_times <- function(dat) {
    
    # Get column names and increment value
    state_cols <- sort(names(dat)[grepl("^state\\d+$", names(dat))])
    increment  <- 1 / (10^dp)
    
    # Step 1: Round all state times
    dat <- dat %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(state_cols), ~ round(., dp)))
    
    # Step 2: For each state (from 2nd onwards), add increment for each prior
    #         state it is equal to, accumulating if multiple ties exist
    for (i in seq_along(state_cols)[-1]) {
      # Pick column
      current_col <- state_cols[i]
      # Pick all prior columns
      prior_cols  <- state_cols[seq_len(i - 1)]
      # Get value from picked column
      current_vals <- dat[[current_col]]
      
      # For every column in prior columns
      for (prior_col in prior_cols) {
        # Get prior column value
        prior_vals   <- dat[[prior_col]]
        # Add increment to current value if equal to any prior value,
        # because sequence goes in order, these previous columns will not be equal to eachother
        current_vals <- dplyr::if_else(
          !is.na(current_vals) & !is.na(prior_vals) & current_vals == prior_vals,
          current_vals + increment,
          current_vals
        )
      }
      
      dat[[current_col]] <- current_vals
    }
    
    return(dat)
  }
  
  ### Run function is rounding == TRUE
  if (rounding == TRUE){
    cohort_in <- adjust_event_times(cohort_in)
  }
  
  ###
  ### Generate censoring times
  ###
  
  ## Generate the censoring times for each individual
  cens_times <- simsurv::simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, 
                                 x = dplyr::select(cohort_in, c("xa", "xb")),
                                 betas = c("xa" = cens_beta_xa, "xb" = cens_beta_xb))
  
  ## Make maximum censoring time the time of max follow
  if (rounding == FALSE){
    cens_times$eventtime <- pmin(cens_times$eventtime, rep(max_follow, nrow(cens_times)))
  } else if (rounding == TRUE){
    cens_times$eventtime <- pmin(round(cens_times$eventtime, dp), rep(max_follow, nrow(cens_times)))
    # if applying rounding, ensure nothing rounded to zero
    cens_times$eventtime <- pmax(round(cens_times$eventtime, dp), 1/(10^dp))
  }
  
  ###
  ### Put data in mstate format
  ###
  
  ### Turn event times into a data.frame and select event times only
  if (dgm_in == 1){
    dat_mstate_temp <- dplyr::select(cohort_in, paste0("state", 1:3))
  } else if (dgm_in == 2){
    dat_mstate_temp <- dplyr::select(cohort_in, paste0("state", 1:4))
  } else if (dgm_in == 3){
    dat_mstate_temp <- dplyr::select(cohort_in, paste0("state", 1:4))
  } else if (dgm_in == 4){
    dat_mstate_temp <- dplyr::select(cohort_in, paste0("state", 1:4))
  }
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat_mstate_temp_noNA <- dat_mstate_temp
  dat_mstate_temp_noNA <- data.frame(t(apply(dat_mstate_temp_noNA, 1, function(x) {ifelse(is.na(x), max_follow, x)})))
  #head(dat_mstate_temp_noNA)
  
  ## Add censoring variables
  dat_mstate_temp_noNA <- cbind(dat_mstate_temp_noNA, matrix(0, ncol = ncol(dat_mstate_temp), nrow = nrow(dat_mstate_temp)))
  colnames(dat_mstate_temp_noNA)[(ncol(dat_mstate_temp)+1):(ncol(dat_mstate_temp)*2)] <-
    paste0("state", 1:ncol(dat_mstate_temp), "_s")
  
  ## If it is not an NA value (from original dataset), set the event/censoring indicator to 1 (event happened)
  for (i in (2:ncol(dat_mstate_temp))){
    dat_mstate_temp_noNA[!is.na(dat_mstate_temp[,i]),(i+ncol(dat_mstate_temp))] <- 1
  }
  #head(dat_mstate_temp_noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat_mstate_temp <- dat_mstate_temp_noNA
  
  ## Add the censoring times
  dat_mstate_temp$cens_time <- cens_times$eventtime
  
  ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
  ## Write function to do this. Writing inside the overarching "convert_mstate_cens", because not planning to use anywhere else
  apply_censoring <- function(dat, cens_time_col = "cens_time") {
    
    # Get column names for event times
    state_time_cols <- names(dat)[grepl("^state\\d+$", names(dat))]
    # Get column names for event indicators
    state_ind_cols  <- names(dat)[grepl("^state\\d+_s$", names(dat))]
    
    # Step 1: Cap event times at censoring time
    dat <- dat %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(state_time_cols),
        ~ dplyr::case_when(
          . < .data[[cens_time_col]]  ~ .,
          . >= .data[[cens_time_col]] ~ .data[[cens_time_col]]
        )
      ))
    
    # Step 2: Zero out indicators where time was capped
    for (ind_col in state_ind_cols) {
      time_col <- sub("_s$", "", ind_col)  # e.g. "state2_s" -> "state2"
      dat <- dat %>%
        dplyr::mutate(!!ind_col := dplyr::case_when(
          .data[[time_col]] < .data[[cens_time_col]]  ~ .data[[ind_col]],
          .data[[time_col]] == .data[[cens_time_col]] ~ 0
        ))
    }
    
    return(dat)
  }
  
  ## Apply the censoring
  dat_mstate_temp <- apply_censoring(dat_mstate_temp)
  
  ## Now need to add baseline data
  dat_mstate_temp$xa <- cohort_in$xa
  dat_mstate_temp$xb <- cohort_in$xb
  dat_mstate_temp$xasq <- cohort_in$xasq
  dat_mstate_temp$xbsq <- cohort_in$xbsq
  dat_mstate_temp$patid <- cohort_in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- readRDS("data/tmat_list.rds")[[dgm_in]]
  
  ## Next detect number of states from column names
  N_states <- ncol(tmat)
  
  ## Put data into wide format using msprep from mstate
  dat_mstate_temp_wide <- mstate::msprep(
    dat_mstate_temp,
    trans  = tmat,
    time   = c(NA, paste0("state", 2:N_states)),
    status = c(NA, paste0("state", 2:N_states, "_s")),
    keep   = c("xa","xb","xasq","xbsq","cens_time", "patid")
  )
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("xa", "xb", "xasq", "xbsq")
  dat_mstate_temp_wide <- mstate::expand.covs(dat_mstate_temp_wide, covs, longnames = FALSE)
  head(dat_mstate_temp_wide)
  
  ###
  ### Create raw data object (not mstate format)
  ###
  
  ### Extract the data and add censoring time
  data_raw <- data.frame(cohort_in, "cens_time" = cens_times$eventtime)
  
  ### Create the dtcens and dtcens_s variables, which will be used by calibmsm for validation
  ### The dtcens variable is the minimum of either entering an absorbing state, or being censored
  ### This will therefore differ depending on the DGM, so will hard code this
  if (dgm_in == 1){
    data_raw <- dplyr::mutate(data_raw,
                              dtcens = pmin(cens_time, state3, na.rm = T),
                              dtcens_s = dplyr::case_when(dtcens == state3 & !is.na(state3) ~ 0,
                                                          TRUE ~ 1))
    # Consider to work when rounding applied XXXX
    # data_raw <- dplyr::mutate(data_raw,
    #                           dtcens = pmin(cens_time, state3, na.rm = T),
    #                           dtcens_s = dplyr::case_when(state3 < cens_time & !is.na(state3) ~ 0,
    #                                                       TRUE ~ 1))
  } else if (dgm_in == 2){
    data_raw <- dplyr::mutate(data_raw,
                              dtcens = pmin(cens_time, state3, state4, na.rm = T),
                              dtcens_s = dplyr::case_when(dtcens == state3 & !is.na(state3) ~ 0,
                                                          dtcens == state4 & !is.na(state4) ~ 0,
                                                          TRUE ~ 1))
    # # Consider to work when rounding applied XXXX
    # data_raw <- dplyr::mutate(data_raw,
    #                           dtcens = pmin(cens_time, state3, state4, na.rm = T),
    #                           dtcens_s = dplyr::case_when(state3 < cens_time & !is.na(state3) ~ 0,
    #                                                       state4 < cens_time & !is.na(state4) ~ 0,
    #                                                       TRUE ~ 1))
  } else if (dgm_in == 3){
    data_raw <- dplyr::mutate(data_raw,
                              dtcens = pmin(cens_time, state3, state4, na.rm = T),
                              dtcens_s = dplyr::case_when(dtcens == state3 & !is.na(state3) ~ 0,
                                                          dtcens == state4 & !is.na(state4) ~ 0,
                                                          TRUE ~ 1))
    # # Consider to work when rounding applied XXXX
    # data_raw <- dplyr::mutate(data_raw,
    #                           dtcens = pmin(cens_time, state3, state4, na.rm = T),
    #                           dtcens_s = dplyr::case_when(state3 < cens_time & !is.na(state3) ~ 0,
    #                                                       state4 < cens_time & !is.na(state4) ~ 0,
    #                                                       TRUE ~ 1))
  } else if (dgm_in == 4){
    data_raw <- dplyr::mutate(data_raw,
                              dtcens = pmin(cens_time, state4, na.rm = T),
                              dtcens_s = dplyr::case_when(dtcens == state4 & !is.na(state4) ~ 0,
                                                          TRUE ~ 1))
    # Consider to work when rounding applied XXXX
    # data_raw <- dplyr::mutate(data_raw,
    #                           dtcens = pmin(cens_time, state4, na.rm = T),
    #                           dtcens_s = dplyr::case_when(state4 < cens_time & !is.na(state4) ~ 0,
    #                                                       TRUE ~ 1))
  }
  
  ### Output object
  return(list("data_mstate" = dat_mstate_temp_wide,
              "tmat" = tmat,
              "data_raw" = data_raw))
  
}


###
### DGM1: Convert data to mstate format, applying censoring, round data to help with computation time
### Non generalised version of convert_mstate_cens that always applies rounding.
### Old, not used in simulation, used for tests only.
###
convert_mstate_cens_dgm1 <- function(cohort_in,
                                     max_follow,
                                     cens_shape,
                                     cens_scale,
                                     cens_beta_xa,
                                     cens_beta_xb,
                                     dp = 3){
  
  #   cohort_in <- cohort[["cohort"]]
  #   max_follow <- cohort[["max_follow"]]
  #   cens_shape <- 1
  #   cens_scale <- 4000
  #   cens_beta_xa <- 0
  #   cens_beta_xb <- 0
  
  ## Turn all event times to integers and create a gap of at least one between all events
  cohort_in <- cohort_in %>%
    dplyr::mutate(state1 = round(state1, dp),
                  state2 = round(state2, dp),
                  state3 = round(state3, dp)) %>%
    dplyr::mutate(state2 = dplyr::case_when(!is.na(state2) & !is.na(state1) & state1 == state2 ~ state2 + 1/(10^dp),
                                            TRUE ~ state2),
                  state3 = dplyr::case_when(!is.na(state3) & !is.na(state1) & state3 == state1 ~ state3 + 1/(10^dp),
                                            !is.na(state3) & !is.na(state2) & state3 == state2 ~ state3 + 1/(10^dp),
                                            TRUE ~ state3)
    )
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat_mstate_temp <- dplyr::select(cohort_in, paste0("state", 1:3))
  
  ## Generate the censoring times for each individual
  cens_times <- simsurv::simsurv("weibull", lambdas = 1/cens_scale, gammas = cens_shape, 
                                 x = dplyr::select(cohort_in, c("xa", "xb")),
                                 betas = c("xa" = cens_beta_xa, "xb" = cens_beta_xb))
  
  # ## If censoring time < 1, set to 1
  # ## Why was I doing this?
  # cens_times$eventtime[cens_times$eventtime < 1] <- 1
  
  ## Make maximum censoring time the time of max follow
  cens_times$eventtime <- pmin(round(cens_times$eventtime, dp), rep(max_follow, nrow(cens_times)))
  
  ## Make minimum censoring time the smallest possible event time above zero
  cens_times$eventtime <- pmax(round(cens_times$eventtime, dp), 1/(10^dp))
  
  ## Create raw dataset for output (the data that was inputted, plus cens_time and patid) and ignore this for rest of function
  data_raw <- data.frame(cohort_in, "cens_times" = cens_times$eventtime)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat_mstate_temp_noNA <- dat_mstate_temp
  dat_mstate_temp_noNA <- data.frame(t(apply(dat_mstate_temp_noNA, 1, function(x) {ifelse(is.na(x), max_follow, x)})))
  #head(dat_mstate_temp_noNA)
  
  ## Add censoring variables
  dat_mstate_temp_noNA <- cbind(dat_mstate_temp_noNA, matrix(0, ncol = ncol(dat_mstate_temp), nrow = nrow(dat_mstate_temp)))
  colnames(dat_mstate_temp_noNA)[(ncol(dat_mstate_temp)+1):(ncol(dat_mstate_temp)*2)] <-
    paste0("state", 1:ncol(dat_mstate_temp), "_s")
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
  for (i in (2:ncol(dat_mstate_temp))){
    dat_mstate_temp_noNA[!is.na(dat_mstate_temp[,i]),(i+ncol(dat_mstate_temp))] <- 1
  }
  #head(dat_mstate_temp_noNA)
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat_mstate_temp <- dat_mstate_temp_noNA
  
  ## Add the censoring times
  dat_mstate_temp$cens_time <- cens_times$eventtime
  
  ## If any events happen after censoring, reduce event time to the censoring time, and set the event indicator to 0
  dat_mstate_temp <- dat_mstate_temp %>%
    dplyr::mutate(state2 = dplyr::case_when(state2 < cens_time ~ state2,
                                            state2 >= cens_time ~ cens_time),
                  state3 = dplyr::case_when(state3 < cens_time ~ state3,
                                            state3 >= cens_time ~ cens_time)) %>%
    dplyr::mutate(state2_s = dplyr::case_when(state2 < cens_time ~ state2_s,
                                              state2 == cens_time ~ 0),
                  state3_s = dplyr::case_when(state3 < cens_time ~ state3_s,
                                              state3 == cens_time ~ 0))
  
  ## Now need to add baseline data
  dat_mstate_temp$xa <- cohort_in$xa
  dat_mstate_temp$xb <- cohort_in$xb
  dat_mstate_temp$xasq <- cohort_in$xasq
  dat_mstate_temp$xbsq <- cohort_in$xbsq
  dat_mstate_temp$patid <- cohort_in$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- mstate::transMat(x = list(c(2,3), c(3), c()),
                           names = paste0("state", 1:3))
  
  
  ##############################################################################
  # Found a bug possibly??
  # When there are not categorical covariates, the msprep "keep =" functionality
  # does not appear to work.
  
  # For example if I run the following code, which does not keep xb (the categorical
  # variable), we get all the xa values jumbled up and in the wrong order
  # cohort.dat.wide.clockreset.test <- msprep(cohort.dat2, trans = tmat, time = c(NA, "state2", "state3"),
  #                                      status = c(NA, "state2_s", "state3_s"), keep = c("xa"))
  # head(cohort.dat.wide.clockreset.test)
  
  # However when we include xb, as in the example below, it is fine
  ##############################################################################
  #
  # THIS BUG NO LONGER APPEARS TO BE HAPPENING - FOUND IN PREVIOUS CODE WHICH IS WHERE ABOVE TEXT IS FROM
  # LEAVING FOR FUTURE REFERENCE TO CHECK FOR THIS ISSUE
  #
  ##############################################################################
  
  ## Now can prepare the data into wide format
  dat_mstate_temp_wide <- mstate::msprep(dat_mstate_temp, trans = tmat,
                                         time = c(NA, paste0("state", 2:3)),
                                         status = c(NA, paste0("state", 2:3, "_s")),
                                         keep = c("xa","xb","xasq","xbsq","cens_time","patid"))
  
  ## Want to expand the covariates to allow different covariate effects per transition
  covs <- c("xa", "xb", "xasq", "xbsq")
  dat_mstate_temp_wide <- mstate::expand.covs(dat_mstate_temp_wide, covs, longnames = FALSE)
  head(dat_mstate_temp_wide)
  
  return(list("data_mstate" = dat_mstate_temp_wide,
              "tmat" = tmat,
              "data_raw" = data_raw))
  
}


###
### Write a function to get predicted transition probabilities for an individual
###

### NB: I don't like to use of "msm_cox_fit" as an argument
### Would rather keep it more general
### Especially given this will also be an input for the function to assess calibration of cause-specific hazards
### 
### Note this wasn't included for assessing calibration of tp's, because users were required to estimate transition probabilities themselves
### We also want to extend the models to not just be cox fit's, which calib_msm already is (because users calculate probabilities themselves),
### but for the cause specific hazards, if it's done internally, shouldn't be called "cox" fit.
### Its a fine name for the model fits in this simulation, as they are cox, but shouldn't be the function argument I don't think.

### Can't change yet as I have preliminary simulations running using est_risk
### Also want to change name to est_risk_probtrans
est_risk_probtrans <- function(id, data_mstate, tmat, covs, msm_cox_fit, t_vec){
  
  ### Paste progress
  print(paste("id =", id, Sys.time(), sep = " "))
  
  ### Define output function
  probtrans_out <- vector("list", length(t_vec))
  
  ### Create temporary dataset to make predictors with
  ## Get location of individual
  person_id_loc <- which(data_mstate$id == id)
  ## Extract dataset with number of rows = number of transition
  temp_data <- data_mstate[rep(person_id_loc[1], max(tmat[!is.na(tmat)])), c("id", covs)]
  ## Add trans
  temp_data$trans <- 1:max(tmat[!is.na(tmat)])
  ## Attribute transition matrix
  attr(temp_data, "trans") <- tmat
  ## Expand the covariates
  temp_data <- mstate::expand.covs(temp_data, covs, longnames = FALSE)
  ## Add strata
  temp_data$strata <- temp_data$trans
  
  ### Calculate transition specific hazards
  msm_fit <- msfit(object = msm_cox_fit, trans = tmat,
                   newdata = temp_data)
  
  ### Calculate probtrans
  probtrans <- probtrans(msm_fit, predt = 0, variance = FALSE)
  
  ### Extract probtrans and se for each yearly interval and put into appropriate datasets
  for (i in 1:length(probtrans_out)){
    
    ### Extract prob trans
    pt_se <- probtrans[[1]] %>%
      ## Extract row corresponding to correct time point
      dplyr::slice(min(which(probtrans[[1]]$time > t_vec[i]))) %>%
      ## Extract all the correct columns
      dplyr::select(paste("pstate",seq_len(ncol(tmat)), sep = ""))
    
    ### Create the row to add by adding patid to pt_se
    pt_se <- data.frame("id" = id, pt_se, "t" = t_vec[i])
    
    ### Bind into output dataset
    probtrans_out[[i]] <- pt_se
  }
  
  return(probtrans_out)
  
}


### Same function but estimating risks using mssample instead of probtrans
### Include variable for M
### Currently taking a data frame and splitting it into list to match output from est_risk, change this in the future
est_risk_mssample <- function(id, data_mstate, tmat, covs, msm_cox_fit, t_vec, m_iter = 1000){
  
  ### Paste progress
  print(paste("id =", id, Sys.time(), sep = " "))
  
  ### Define output function
  probtrans_out <- vector("list", length(t_vec))
  
  ### Create temporary dataset to make predictors with
  ## Get location of individual
  person_id_loc <- which(data_mstate$id == id)
  ## Extract dataset with 27 rows of this individual
  temp_data <- data_mstate[rep(person_id_loc[1], max(tmat[!is.na(tmat)])), c("id", covs)]
  ## Add trans
  temp_data$trans <- 1:max(tmat[!is.na(tmat)])
  ## Attribute transition matrix
  attr(temp_data, "trans") <- tmat
  ## Expand the covariates
  temp_data <- mstate::expand.covs(temp_data, covs, longnames = FALSE)
  ## Add strata
  temp_data$strata <- temp_data$trans
  
  ### Calculate transition specific hazards
  msm_fit <- msfit(object = msm_cox_fit, trans = tmat,
                   newdata = temp_data)
  
  ### Calculate probtrans
  probtrans <- mstate::mssample(Haz = msm_fit$Haz[msm_fit$Haz$time <= (max(t_vec)+0.1), ], 
                                trans = tmat, 
                                clock = "reset", 
                                tvec = t_vec, 
                                M = m_iter)
  
  ### Add id and rename time to t, and re-order
  probtrans <- dplyr::mutate(probtrans, id = id) |>
    dplyr::rename(t = time) |>
    dplyr::relocate(id, paste("pstate",1:3, sep = ""), t)
  
  ### Put into same format as output from est_risk
  probtrans_out <- unname(split(probtrans, seq(nrow(probtrans))))
  
  return(probtrans_out)
  
}

###
### Write a function to get predicted transition probabilities for an individual, when model was developed using flexible parametric
###

### NB: calibmsm::calib_msm is agnostic to way model is developed, but calibration of cause-specific hazards will not be?
### Consider this!
est_risk_flexsurv <- function(id, data_raw, tmat, msm_fit, t_vec, m_iter = 10000){
  
  ### Paste progress
  print(paste("id =", id, Sys.time(), sep = " "))
  
  ## Select individual for predictions (only works when id = 1:nrow())
  temp_data <- data_raw[id, ]
  ## Safer but slower code to do the reduction
  # temp_data <- data_raw[fastmatch::fmatch(id, data_raw$id), ]
  # if (anyNA(temp_data$id)) stop("id not found in data raw")
  
  ### Calculate probtrans
  flexsurv_pt <- flexsurv::pmatrix.simfs(msm_fit, trans = tmat, t = t_vec, M = m_iter, newdata = temp_data)
  
  ### A single prediction time returns a matrix rather than 3D array
  ### Add the missing time dimension so the extract below works for when length(t_vec)==1
  if (length(dim(flexsurv_pt)) == 2){
    flexsurv_pt <- array(
      flexsurv_pt, 
      dim = c(dim(flexsurv_pt), 1),
      dimnames = c(dimnames(flexsurv_pt), list(NULL))
    )
  }
  
  ### Reduce to predictions out of state 1, and add id and t variables
  ### Get into format that is compatible with the rest of code
  probtrans_out <- 
    lapply(1:dim(flexsurv_pt)[3], function(x){
      
      ### Reduce to first rows
      out <- flexsurv_pt[1,,x]
      
      ### Create data.frame
      out <- as.data.frame(t(out))
      
      ### Name probability columns by state number
      names(out) <- paste0("pstate", seq_len(ncol(out)))
      
      ### Add extra variables
      out <- dplyr::mutate(out, id = id, t = t_vec[x]) |>
        dplyr::relocate(id) 
      
      return(out)})
  
  return(probtrans_out)
  
}


est_risk_flexsurv_OLD <- function(id, data_raw, tmat, msm_fit, t_vec, m_iter = 10000){
  
  ### Paste progress
  print(paste("id =", id, Sys.time(), sep = " "))
  
  ## Select individual for predictions (only works when id = 1:nrow())
  temp_data <- data_raw[id, ]
  ## Safer but slower code to do the reduction
  # temp_data <- data_raw[fastmatch::fmatch(id, data_raw$id), ]
  # if (anyNA(temp_data$id)) stop("id not found in data raw")
  
  ### Calculate probtrans
  flexsurv_pt <- flexsurv::pmatrix.simfs(msm_fit, trans = tmat, t = t_vec, M = m_iter, newdata = temp_data)
  
  ### Reduce to predictions out of state 1, and add id and t variables
  ### Get into format that is compatible with the rest of code
  probtrans_out <- 
    lapply(1:dim(flexsurv_pt)[3], function(x){
      ### Reduce to first rows
      out <- flexsurv_pt[1,,x]
      ### Create data.frame
      out <- as.data.frame(t(out))
      ### Add extra variables
      out <- dplyr::mutate(out, id = id, t = t_vec[x]) |>
        dplyr::relocate(id) |>
        ### Rename state to pstate so compatible with rest of code
        dplyr::rename_with(~ paste0("p", .), dplyr::starts_with("state"))
      return(out)})
  
  return(probtrans_out)
  
}

###
### Write a function to get CSH hazard calibration plots
###
create_plot_csh_old <- function(data_mstate, t){
  
  ### Create plots list object
  plots_list <- vector("list", 3)
  
  for (transition in 1:3){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Fit the cox model, using lp1 as the only predictors
    calib.model <- coxph(Surv(time, status) ~ lp, data = data_csh)
    
    ### Now generate predicted-observed values
    ### Extract the baseline hazard at t = 0.25
    basehaz <- basehaz(calib.model, centered = FALSE)
    basehaz <- basehaz$hazard[min(which(basehaz$time >= t))]
    # if (transition %in% c(1,2)){
    #   basehaz <- basehaz$hazard[min(which(basehaz$time >= 0.25))]
    # } else if (transition == 3){
    #   basehaz <- basehaz$hazard[min(which(basehaz$time >= 0.25))]
    # }
    
    ### Calculate predicted-observed risks
    data_csh$pred.obs <- 1-exp(-basehaz*exp(data_csh$lp))
    
    ### Plot predo.obs1, vs tp12
    plots_list[[transition]] <- ggplot(aes(x = netrisk, y = pred.obs), data = data_csh) +
      geom_line(color = "red") +
      geom_abline(slope = 1, intercept = 0, lty = "dotted") +
      ggtitle(paste("Transition ", c("1->2", "1->3", "2->3")[transition], sep = "")) +
      xlab("Net risk") + ylab("Observed net risk")
    
  }
  
  ### Combine plots into single ggplot
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = 1, ncol = 3, common.legend = TRUE)
  
  return(plots_comb)
  # ggsave(paste("plots/calibplot_CSH_scen", scenario, ".png", sep = ""), plots_comb, dpi = 300)
  
}

create_plot_csh_old2 <- function(data_mstate, t){
  
  ### Create plots list object
  calib_list <- vector("list", 3)
  
  for (transition in 1:3){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Create calibration object
    calib_object <- est_calib_csh_ph(data = data_csh, t = t, surv = data_csh$predsurv, nk = 4)
    
    ### Add title
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggtitle(paste("Transition ", c("1->2", "1->3", "2->3")[transition], sep = ""))
    
    ### Add to output
    calib_list[[transition]] <- calib_object
    
  }
  
  ### Create list of plots
  plots_list <- lapply(calib_list, function(x){x[["plot"]]})
  
  ### Combine plots into single ggplot
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = 1, ncol = 3, common.legend = TRUE)
  
  ### Create output object
  output_object <- list("calib_list" = calib_list, "plots_comb" = plots_comb)
  return(output_object)
  
}

###
### Function to estimate calibration of the cause-specific hazards functions.
### It will generate predictions from fit_msm for individuals in the validation cohort.
### then apply est_calib_csh_ph to calculate calibration curves, ICI, E50 and E90 for each transition.
##
### This will need generalising to:
### - interval censored data
### - any number of transitions
### - defining which transitions you want to look at
### - using pseudo-values instead of proportional hazards
### - when the cox fit hasn't been fitted with "model = TRUE"
###
### - Need to think about whether this works with variables that change at intermediate transitions (makes particular sense for clock reset)
### NB: I have generalised calib_csh_flexsurv to any number of transitions, but it was slightly easier,
### given each model is stored in it's own list element. Probably want to align these two into the same function eventually?
###
### NB: Should change name to calib_csh_cox, but need to check not used anywhere in code before doing this
### (I think its now obselete since move to flexible parametric)
calib_csh <- function(data_mstate, fit_msm, t, nk = 4){
  
  ### Create basehaz table (uncentered)
  basehaz_table <- basehaz(fit_msm, centered = FALSE)
  
  ### Get baseline hazards at t
  bhaz1 <- basehaz_table$hazard[max(which(basehaz_table$time <= t & basehaz_table$strata == "trans=1"))]
  bhaz2 <- basehaz_table$hazard[max(which(basehaz_table$time <= t & basehaz_table$strata == "trans=2"))]
  bhaz3 <- basehaz_table$hazard[max(which(basehaz_table$time <= t & basehaz_table$strata == "trans=3"))]
  
  ### Get linear predictor for each transition
  data_mstate$lp <- predict(fit_msm, newdata = data_mstate, type = "lp", reference = "zero")
  
  ### Turn into risk scores for the specific transition (net risks)
  data_mstate <- dplyr::mutate(data_mstate, 
                               predsurv = exp(-exp(lp)*c(bhaz1, bhaz2, bhaz3)[trans]),
                               predrisk = 1 - exp(-exp(lp)*c(bhaz1, bhaz2, bhaz3)[trans])
  )
  
  ### Create plots list object
  calib_list <- vector("list", 3)
  names(calib_list) <- c("1->2", "1->3", "2->3")
  
  for (transition in 1:3){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Create calibration object
    calib_object <- est_calib_csh_ph(data = data_csh, t = t, surv = data_csh$predsurv, nk = nk)
    
    ### Add title
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggtitle(paste("Transition ", c("1->2", "1->3", "2->3")[transition], sep = ""))
    
    ### Add to output
    calib_list[[transition]] <- calib_object
    
  }
  
  ### Create list of plots
  plots_list <- lapply(calib_list, function(x){x[["plot"]]})
  
  ### Combine plots into single ggplot
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = 1, ncol = 3, common.legend = TRUE)
  
  ### Create output object
  output_object <- list("calib_list" = calib_list, "plots_comb" = plots_comb)
  return(output_object)
  
}


calib_csh_attempt2_DEL <- function(data_valid_mstate, data_devel_mstate, fit_msm, t){
  
  data_valid_mstate <- df_valid_mstate
  data_devel_mstate <- df_devel_mstate
  fit_msm <- msm_cox_fit
  t <- 5
  
  ### Function to extract baseline hazard at time t for a given transition
  get_bhaz_csh <- function(transition, time){
    
    ### Subset development data
    data_devel_subset <- subset(data_devel_mstate, trans == transition) 
    
    print("test1")
    ### Need to fit the same cox model in stratified data to get the appropriate cumulative baseline hazard
    ### Remove strata term from formula as we have subsetted
    fit_transition <- coxph(update(msm_cox_fit$formula, ~ . - strata(trans)), data = data_devel_subset, model = TRUE) # NB: Coefficients for other transition variables will be NA
    print("test2")
    ### Extract cumulative basehaz
    basehaz_obj_transition <- basehaz(fit_transition, centered = FALSE)
    basehaz_t_transition <- basehaz_obj_transition$hazard[max(which(basehaz_obj_transition$time <= t))]
    
    return(basehaz_t_transition)
    
  }
  
  basehaz_t_vec <- lapply(1:3, function(x){get_bhaz_csh(x, time = t)})
  basehaz_t_vec
  
  ### Get baseline hazards at t
  bhaz1 <- basehaz_table$hazard[max(which(basehaz_table$time <= t & basehaz_table$strata == "trans=1"))]
  bhaz2 <- basehaz_table$hazard[max(which(basehaz_table$time <= t & basehaz_table$strata == "trans=2"))]
  bhaz3 <- basehaz_table$hazard[max(which(basehaz_table$time <= t & basehaz_table$strata == "trans=3"))]
  
  ### Get linear predictor for each transition
  data_mstate$lp <- predict(fit_msm, newdata = data_mstate, type = "lp", reference = "zero")
  
  ### Turn into risk scores for the specific transition (net risks)
  data_mstate <- dplyr::mutate(data_mstate, 
                               predsurv = exp(-exp(lp)*c(bhaz1, bhaz2, bhaz3)[trans]),
                               predrisk = 1 - exp(-exp(lp)*c(bhaz1, bhaz2, bhaz3)[trans])
  )
  
  ### Create plots list object
  calib_list <- vector("list", 3)
  names(calib_list) <- c("1->2", "1->3", "2->3")
  
  for (transition in 1:3){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Create calibration object
    calib_object <- est_calib_csh_ph(data = data_csh, t = t, surv = data_csh$predsurv, nk = 4)
    
    ### Add title
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggtitle(paste("Transition ", c("1->2", "1->3", "2->3")[transition], sep = ""))
    
    ### Add to output
    calib_list[[transition]] <- calib_object
    
  }
  
  ### Create list of plots
  plots_list <- lapply(calib_list, function(x){x[["plot"]]})
  
  ### Combine plots into single ggplot
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = 1, ncol = 3, common.legend = TRUE)
  
  ### Create output object
  output_object <- list("calib_list" = calib_list, "plots_comb" = plots_comb)
  return(output_object)
  
}


###
### Function to estimate calibration of the cause-specific hazards functions.
### This is specific to a multistate model developed using "flexsurv"
### It will generate predictions from fit_msm for individuals in the validation cohort.
### then apply est_calib_csh_ph to calculate calibration curves, ICI, E50 and E90 for each transition.
### Can be applied to a model with any number of transitions
##
### This will need generalising to:
### - interval censored data
### - defining which transitions you want to look at
### - using pseudo-values instead of proportional hazards
### - when the cox fit hasn't been fitted with "model = TRUE"
###
### - Need to think about whether this works with variables that change at intermediate transitions (makes particular sense for clock reset)
calib_csh_flexsurv <- function(data_mstate, fit_msm, t, nk = 4){
  
  ### Get number of transitions
  transitions_n <- max(data_mstate$trans)
  
  ### Create plots list object
  calib_list <- vector("list", transitions_n)
  
  ###
  ### Get names of transitions
  ###
  
  ### NB: ADD A TEST TO ENSURE THE TRANSITION MATRIX ATTRIBUTE IS THERE
  
  ### Define tmat
  tmat <- attributes(data_mstate)$trans
  
  # Pre-allocate the result vector
  transition_names <- character(transitions_n)
  
  # Loop over each integer value 1..max_val and find which cell contains it
  for (transition in seq_len(transitions_n)) {
    
    # which(..., arr.ind = TRUE) returns the row/col index of the matching cell
    idx <- which(tmat == transition, arr.ind = TRUE)
    
    # Extract the row and column *names* using the index
    row_name <- rownames(tmat)[idx[, 1]]
    col_name <- colnames(tmat)[idx[, 2]]
    
    # Build the "from -> to" string
    transition_names[transition] <- paste0(row_name, " -> ", col_name)
    
  }
  
  ###
  ### Loop through for each transition and assess calibration
  ###
  for (transition in 1:transitions_n){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Get predicted risks
    ### Note, the transitions must match the order of the fitted models in fit_msm
    data_csh$predsurv <- predict(fit_msm[[transition]], newdata = data_csh, type = "survival", times = t)$.pred_survival
    
    ### Create calibration object
    calib_object <- est_calib_csh_ph(data = data_csh, t = t, surv = data_csh$predsurv, nk = nk)
    
    ### Add title
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggplot2::ggtitle(paste("Transition: ", transition_names[transition], sep = "")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(1)))
    
    ### Add to output
    calib_list[[transition]] <- calib_object
    
  }
  
  ### Create list of plots
  plots_list <- lapply(calib_list, function(x){x[["plot"]]})
  
  ### Combine plots into single ggplot
  # Get number of rows (always assume three columns for consistency)
  plot_nrows <- ceiling(transitions_n/3)
  # Create plots
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = plot_nrows, ncol = 3, common.legend = TRUE)
  
  ### Create output object
  output_object <- list("calib_list" = calib_list, "plots_comb" = plots_comb)
  return(output_object)
  
}


calib_csh_flexsurv_test <- function(data_mstate, fit_msm, t){
  
  ### Get number of transitions
  transitions_n <- max(data_mstate$trans)
  
  ### Create plots list object
  calib_list <- vector("list", transitions_n)
  
  ###
  ### Get names of transitions
  ###
  
  ### NB: ADD A TEST TO ENSURE THE TRANSITION MATRIX ATTRIBUTE IS THERE
  
  ### Define tmat
  tmat <- attributes(data_mstate)$trans
  
  # Pre-allocate the result vector
  transition_names <- character(transitions_n)
  
  # Loop over each integer value 1..max_val and find which cell contains it
  for (transition in seq_len(transitions_n)) {
    
    # which(..., arr.ind = TRUE) returns the row/col index of the matching cell
    idx <- which(tmat == transition, arr.ind = TRUE)
    
    # Extract the row and column *names* using the index
    row_name <- rownames(tmat)[idx[, 1]]
    col_name <- colnames(tmat)[idx[, 2]]
    
    # Build the "from -> to" string
    transition_names[transition] <- paste0(row_name, " -> ", col_name)
    
  }
  
  ###
  ### Loop through for each transition and assess calibration
  ###
  for (transition in 1:transitions_n){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Get predicted risks
    ### Note, the transitions must match the order of the fitted models in fit_msm
    data_csh$predsurv <- predict(fit_msm[[transition]], newdata = data_csh, type = "survival", times = t)$.pred_survival
    
    ### Create calibration object
    calib_object <- est_calib_csh_ph(data = data_csh, t = t, surv = data_csh$predsurv, nk = 4)
    
    ### Add title
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggplot2::ggtitle(paste("Transition: ", transition_names[transition], sep = "")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(1)))
    
    ## Add invisible scatter for marginal
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs),
                          col = grDevices::rgb(0, 0, 0, alpha = 0)) +
      ## Remove legend
      ggplot2::theme(legend.position = "none")
    
    ### Add marginal
    calib_object[["plot"]] <- ggExtra::ggMarginal(calib_object[["plot"]], 
                                                  margins = "x", 
                                                  size = 5, 
                                                  type = "density", 
                                                  colour = "red")
    
    ### Add to output
    calib_list[[transition]] <- calib_object
    
  }
  
  ### Create list of plots
  plots_list <- lapply(calib_list, function(x){x[["plot"]]})
  
  ### Combine plots into single ggplot
  # Get number of rows (always assume three columns for consistency)
  plot_nrows <- ceiling(transitions_n/3)
  # Create plots
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = plot_nrows, ncol = 3, common.legend = TRUE)
  
  ### Create output object
  output_object <- list("calib_list" = calib_list, "plots_comb" = plots_comb)
  return(output_object)
  
}


calib_csh_flexsurv_OLD_not_general <- function(data_mstate, fit_msm, t){
  
  ### Get number of transitions
  ### Need to generalise this and the names for package
  transitions_n <- 3
  
  ### Create plots list object
  calib_list <- vector("list", 3)
  names(calib_list) <- c("1->2", "1->3", "2->3")
  
  for (transition in 1:transitions_n){
    
    ### Create a dataset upon fit to fit the calibration model, where we will model time until outcome
    print(transition)
    data_csh <- subset(data_mstate, trans == transition)
    
    ### Get predicted risks
    data_csh$predsurv <- predict(fit_msm[[transition]], newdata = data_csh, type = "survival", times = t)$.pred_survival
    
    ### Create calibration object
    calib_object <- est_calib_csh_ph(data = data_csh, t = t, surv = data_csh$predsurv, nk = 4)
    
    ### Add title
    calib_object[["plot"]] <- calib_object[["plot"]] +
      ggtitle(paste("Transition ", c("1->2", "1->3", "2->3")[transition], sep = ""))
    
    ### Add to output
    calib_list[[transition]] <- calib_object
    
  }
  
  ### Create list of plots
  plots_list <- lapply(calib_list, function(x){x[["plot"]]})
  
  ### Combine plots into single ggplot
  plots_comb <- ggpubr::ggarrange(plotlist = plots_list, nrow = 1, ncol = 3, common.legend = TRUE)
  
  ### Create output object
  output_object <- list("calib_list" = calib_list, "plots_comb" = plots_comb)
  return(output_object)
  
}

###
### Assessing calibration using proportional hazards regression approach (graphical calibration curves, Austin et al)
###
est_calib_csh_ph <- function(data, fit, bhaz, t, surv = NULL, nk = 4, pred.plot.range = NULL, plot = TRUE){
  
  ### Get the survival probabilities
  if (is.null(surv)){
    data$surv <- as.numeric(est_surv(newdata = data, fit = fit, bhaz = bhaz, t = t))
  } else {
    data$surv <- as.numeric(surv)
  }
  
  ### Add complementary log-log of predicted survival probabilities to data
  data$cloglog <- log(-log(data$surv))
  
  ### Fit calibration model
  fit.calib <- survival::coxph(survival::Surv(time, status) ~ rms::rcs(cloglog, nk), data = data)
  bhaz.calib <- survival::basehaz(fit.calib, centered = TRUE)
  
  ###
  ### Generate predicted observed values,
  ###
  
  ### Create a dataframe with the calibration data
  pred.obs <- est_surv(newdata = data, fit = fit.calib, bhaz = bhaz.calib, t = t)
  output.data <- data.frame("id" = data$id, "pred.obs" = 1 - pred.obs, "pred" = 1 - data$surv)
  
  ### Calculate ICI, E50 and E90
  ICI <- mean(abs(output.data$pred.obs - output.data$pred))
  E50 <- median(abs(output.data$pred.obs - output.data$pred))
  E90 <- as.numeric(quantile(abs(output.data$pred.obs - output.data$pred), probs = .9, na.rm = TRUE))
  
  ### Create dataframe for plotting
  ### Do so over the range of values pred.plot.range, if specified
  if (!is.null(pred.plot.range)){
    ### Create temporary data frame
    tmp.data <- data.frame("pred" = pred.plot.range, cloglog = log(-log(1 - pred.plot.range)))
    
    ### Calculate predicted observed values over pred.plot.range
    pred.obs <- 1 - est_surv(newdata = tmp.data, fit = fit.calib, bhaz = bhaz.calib, t = t)
    
    ### Create plot data
    plot.data <- data.frame("pred.obs" = pred.obs, "pred" = tmp.data$pred)
    
  } else {
    
    ### Create plot data
    plot.data <- output.data
    
    ### Now set output.data to NA, 
    ### this is so we don't save it twice (it will be an element in the plot object, which will be outputted)
    output.data <- NA
  }
  
  ### Create plot
  if (plot == TRUE){
    plot <- ggplot2::ggplot(data = plot.data) +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), color = "red") +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") + 
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk")
  } else {
    plot <- NULL
  }
  
  ### Create output.object
  output.object <- list("plot" = plot,
                        "ICI" = ICI,
                        "E50" = E50,
                        "E90" = E90,
                        "calib_data" = output.data)
  
  return(output.object)
  
}

###
### Write a function to estimate survival probabilities based on a fitted model (fit) and baseline hazard (bhaz).
### Baseline hazard must have been fitted using basehaz(surv.obj, centered = TRUE).
###
est_surv <- function(newdata, fit, bhaz = NULL, t){
  
  ### Get bhaz if not specified
  if (is.null(bhaz)){bhaz <- survival::basehaz(fit, centered = TRUE)}
  
  ### Get the lp
  lp <- predict(fit, newdata = newdata, reference = "sample")
  
  ### Get the linear predictor for ne wdata
  surv <- as.numeric(exp(-exp(lp)*bhaz$hazard[max(which(bhaz$time <= t))]))
  
  return(surv)
  
}

####################################################
### Define functions for creating combined plots ###
####################################################
### XXXX

###
### Helper function to create ggtitle subtitles
###
make_title <- function(title_text){
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = title_text, 
                      fontface = "bold", size = 4) +
    ggplot2::theme_void()
}

###
### Helper function to thin calibration curves
###
thin_curve_data <- function(dat, n = 5000){
  
  # Return unchanged if already small
  if (nrow(dat) <= n){
    return(dat)
  }
  
  # Select evenly spaced rows
  idx <- unique(round(seq(1, nrow(dat), length.out = n)))
  
  dat[idx, , drop = FALSE]
  
}

###
### Helper function to calculate ICI, E50, E90
###
calc_metrics <- function(pred, obs){
  diffs <- abs(pred - obs)
  c(ICI = mean(diffs, na.rm = TRUE), E50 = as.numeric(median(diffs, na.rm = TRUE)), E90 = as.numeric(quantile(diffs, 0.9, na.rm = TRUE)))
}

###
### Create histogram
###
make_histogram <- function(plotdata){
  
  ggplot2::ggplot(
    plotdata,
    ggplot2::aes(x = pred)
  ) +
    ggplot2::geom_density(
      fill = "grey80",
      colour = "grey50"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 0.25,
      plot.margin = ggplot2::margin(t = 1, r = 5, b = 0, l = 5),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::ylab("Density")
  
}

###
### Create TP panel
###
make_tp_panel <- function(plotdata, state_num, thin_n = NULL){
  
  # Calculate metrics
  metrics <- calc_metrics(plotdata$pred, plotdata$obs)
  
  # Thin only the calibration curve
  curve_dat <- plotdata |>
    dplyr::arrange(pred) 
  
  ### Thin data
  if (!is.null(thin_n)){
    curve_dat <- thin_curve_data(curve_dat, n = thin_n)
    }
  
  # Get max value for x-axis
  xmax <- max(curve_dat$pred, na.rm = TRUE)
  
  # Get margin as a fraction of the panels range (for inserting calibration metrics onto plot)
  metric_margin <- 0.025*xmax
  
  # Histogram
  hist_plot <- make_histogram(plotdata) +
    ggplot2::ggtitle(paste("State: ", state_num, sep = "")) +
    ggplot2::xlim(c(0,xmax))
  
  # Calibration curve
  calib_plot <-
    ggplot2::ggplot(curve_dat) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = pred,
        y = obs
      ),
      colour = "red"
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::coord_fixed(xlim = c(0,xmax),
                         ylim = c(0,xmax)) +
    # Add calibration metrics
    ggplot2::annotate("text", x = xmax - metric_margin, y = 0 + metric_margin,
                      label = paste0("ICI: ", round(metrics[["ICI"]], 3), "\n",
                                     "E50: ", round(metrics[["E50"]], 3), "\n",
                                     "E90: ", round(metrics[["E90"]], 3)),
                      hjust = 1, vjust = 0, size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1,
                   plot.margin = ggplot2::margin(t = 0)) +
    ggplot2::xlab("Predicted risk") +
    ggplot2::ylab("Observed risk")
  
  # Combine histogram and calibration plot
  hist_plot /
    calib_plot +
    patchwork::plot_layout(
      heights = c(1, 4)
    )
  
}

###
### Create CSH panel
###
make_csh_panel <- function(plotdata, transition_num){
  
  # Calculate metrics
  metrics <- calc_metrics(plotdata$pred, plotdata$pred.obs)
  
  ### Read in tmat
  tmat <- readRDS("data/tmat_list.rds")[[dgm]]
  
  ### Get number of transitions
  transitions_n <- sum(!is.na(tmat))
  
  ### Pre-allocate the transition_names vector
  transition_names <- character(transitions_n)
  
  ### Create transition names
  # Loop over each integer value 1..max_val and find which cell contains it
  for (transition in seq_len(transitions_n)) {
    
    # which(..., arr.ind = TRUE) returns the row/col index of the matching cell
    idx <- which(tmat == transition, arr.ind = TRUE)
    
    # Extract the row and column *names* using the index
    row_name <- rownames(tmat)[idx[, 1]]
    col_name <- colnames(tmat)[idx[, 2]]
    
    # Build the "from -> to" string
    transition_names[transition] <- paste0(row_name, " -> ", col_name)
    
  }
  
  # Thin only the calibration curve
  curve_dat <- plotdata |>
    dplyr::arrange(pred) |>
    thin_curve_data(n = thin_n)
  
  # Get max value for x-axis
  xmax <- max(curve_dat$pred, na.rm = TRUE)
  
  # Get margin as a fraction of the panels range (for inserting calibration metrics onto plot)
  metric_margin <- 0.025*xmax
  
  # Histogram
  hist_plot <- make_histogram(plotdata) + 
    ggplot2::ggtitle(paste("Transition: ", transition_names[transition_num], sep = ""))
  
  # Calibration curve
  calib_plot <-
    ggplot2::ggplot(curve_dat) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = pred,
        y = pred.obs
      ),
      colour = "red"
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::coord_fixed(xlim = c(0,xmax),
                         ylim = c(0,xmax)) +
    # Add calibration metrics
    ggplot2::annotate("text", x = xmax - metric_margin, y = 0 + metric_margin,
                      label = paste0("ICI: ", round(metrics[["ICI"]], 3), "\n",
                                     "E50: ", round(metrics[["E50"]], 3), "\n",
                                     "E90: ", round(metrics[["E90"]], 3)),
                      hjust = 1, vjust = 0, size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1,
                   plot.margin = ggplot2::margin(t = 0)) +
    ggplot2::xlab("Predicted risk") +
    ggplot2::ylab("Observed risk")
  
  # Combine histogram and calibration plot
  hist_plot /
    calib_plot +
    patchwork::plot_layout(
      heights = c(1, 4)
    )
  
}

###
### Function to combined figure for one scenario
###
create_combined_calibration_plot <- function(dgm,
                                             scenario,
                                             t_in,
                                             model_updated = FALSE,
                                             thin_n = NULL){
  
  ### Read in sim_inputs for the DGM, devel num and valid num, used for plot titles
  sim_inputs <- readRDS(paste("data/sim_inputs_dgm", dgm, ".rds", sep = ""))
  devel_num <- as.numeric(sim_inputs[scenario, "devel_num"])
  valid_num <- as.numeric(sim_inputs[scenario, "valid_num"])
  
  ###
  ### Read in calibration data
  ###
  if (isFALSE(model_updated)){
    
    calib_tp_obj <- readRDS(paste0("data/calib_object_tp_dgm",dgm,"_s",scenario,"_t",t_in,".rds"))
    calib_csh_obj <- readRDS(paste0("data/calib_object_csh_dgm",dgm,"_s",scenario,"_t",t_in,".rds"))
    
  } else if (isTRUE(model_updated)){
    
    calib_tp_obj <- readRDS(paste0("data/calib_object_tp_model_update_dgm",dgm,"_s",scenario,"_t",t_in,".rds"))
    calib_csh_obj <- readRDS(paste0("data/calib_object_csh_model_update_dgm",dgm,"_s",scenario,"_t",t_in,".rds"))
    
  }
  
  
  ###
  ### Create TP panels
  ###
  tp_panels <- lapply(
    seq_along(calib_tp_obj$plotdata),
    function(i){
      make_tp_panel(
        plotdata = calib_tp_obj$plotdata[[i]],
        state_num = i,
        thin_n = thin_n
      )
    }
  )
  
  ###
  ### Create CSH panels
  ###
  csh_panels <- lapply(
    seq_along(calib_csh_obj),
    function(i){
      make_csh_panel(
        plotdata = calib_csh_obj[[i]]$plotdata,
        transition_num = i
      )
    }
  )
  
  ###
  ### Create titles
  ###
  tp_title <- make_title(paste0("Calibration of the Transition Probabilities"))
  csh_title <- make_title(paste0("Calibration of the Cause-Specific Failure Functions"))
  
  ###
  ### Helper to build a row of panels, padded with equal-width spacers on
  ### either side so it visually centers against a row of n_ref panels.
  ### Building each row as its own patchwork object (rather than one shared
  ### design string for the whole figure) keeps its internal width-solving
  ### independent of the other rows, since the panels use coord_fixed() +
  ### aspect.ratio = 1 and mixing aspect-locked panels into one shared
  ### design grid seems to distort column widths across rows.
  ###
  make_centered_row <- function(panels_list, n_ref){
    
    n_panels <- length(panels_list)
    pad_total <- n_ref - n_panels
    
    if (pad_total <= 0){
      # Row already fills the reference width (or exceeds it) - no padding needed
      row <- patchwork::wrap_plots(panels_list, nrow = 1, ncol = n_panels)
      return(row)
    }
    
    pad_each <- pad_total / 2
    
    row <- patchwork::wrap_plots(
      c(list(patchwork::plot_spacer()), panels_list, list(patchwork::plot_spacer())),
      nrow = 1,
      widths = c(pad_each, rep(1, n_panels), pad_each)
    )
    
    return(row)
    
  }
  
  ###
  ### Layout by DGM
  ###
  
  if (dgm == 1){
    grid_cols <- 3
    grid_rows <- 2
    tp_row  <- make_centered_row(tp_panels, n_ref = 3)
    csh_row <- make_centered_row(csh_panels, n_ref = 3)
    final_plot <- (tp_title / tp_row / csh_title / csh_row) +
      patchwork::plot_layout(heights = c(0.05, 1, 0.05, 1))
  }
  
  if (dgm == 2){
    grid_cols <- 4
    grid_rows <- 2
    tp_row  <- make_centered_row(tp_panels, n_ref = 4)
    csh_row <- make_centered_row(csh_panels, n_ref = 4)
    final_plot <- (tp_title / tp_row / csh_title / csh_row) +
      patchwork::plot_layout(heights = c(0.05, 1, 0.05, 1))
  }
  
  if (dgm == 3){
    grid_cols <- 4
    grid_rows <- 2
    tp_row  <- make_centered_row(tp_panels, n_ref = 4)
    csh_row <- make_centered_row(csh_panels, n_ref = 4)
    final_plot <- (tp_title / tp_row / csh_title / csh_row) +
      patchwork::plot_layout(heights = c(0.05, 1, 0.05, 1))
  }
  
  if (dgm == 4){
    # CSH panels split across two lines (3 then 2), each centered against
    # the 4-wide TP row independently
    grid_cols <- 4
    grid_rows <- 3
    tp_row   <- make_centered_row(tp_panels, n_ref = 4)
    csh_row1 <- make_centered_row(csh_panels[1:3], n_ref = 4)
    csh_row2 <- make_centered_row(csh_panels[4:5], n_ref = 4)
    final_plot <- (tp_title / tp_row / csh_title / csh_row1 / csh_row2) +
      patchwork::plot_layout(heights = c(0.05, 1, 0.05, 1, 1))
  }
  
  ###
  ### Save figure
  ###
  
  ### Dynamically choose sizing
  base_panel_width <- 3.5
  height_width_ratio <- 1.25
  title_row_height <- base_panel_width*height_width_ratio*0.05
  width_in <- grid_cols*base_panel_width
  height_in <-
    # Main bulk
    grid_rows*base_panel_width*height_width_ratio + 
    # Space for subtitles
    2*title_row_height + 
    # space for title
    0.5
  
  ### Save
  if (isFALSE(model_updated)){
    png(
      filename = paste0("figures/calib_plot_comb_DGM",dgm,"_s",scenario,"_t",t_in,".png"),
      width = width_in,
      height = height_in,
      units = "in",
      res = 90
    )
  } else if (isTRUE(model_updated)){
    png(
      filename = paste0("figures/calib_plot_comb_model_updated_DGM",dgm,"_s",scenario,"_t",t_in,".png"),
      width = width_in,
      height = height_in,
      units = "in",
      res = 90
    )
  }
  
  print(final_plot)
  grDevices::dev.off()
  
  ### Save high res figures
  if (isFALSE(model_updated)){
    if (dgm == 1 & devel_num == 5 & valid_num == 2 & t_in == 5){
      png(filename = paste0("figures/Figure2.png"),
          width = width_in,
          height = height_in,
          units = "in",
          res = 300)
      print(final_plot)
      grDevices::dev.off()
    } else if (dgm == 2 & devel_num == 1 & valid_num == 3 & t_in == 5){
      png(filename = paste0("figures/Figure3.png"),
          width = width_in,
          height = height_in,
          units = "in",
          res = 300)
      print(final_plot)
      grDevices::dev.off()
    } else if (dgm == 3 & devel_num == 1 & valid_num == 3 & t_in == 5){
      png(filename = paste0("figures/Figure4.png"),
          width = width_in,
          height = height_in,
          units = "in",
          res = 300)
      print(final_plot)
      grDevices::dev.off()
    } else if (dgm == 4 & devel_num == 1 & valid_num == 3 & t_in == 5){
      png(filename = paste0("figures/Figure5.png"),
          width = width_in,
          height = height_in,
          units = "in",
          res = 300)
      print(final_plot)
      grDevices::dev.off()
    }
  }
  
}

