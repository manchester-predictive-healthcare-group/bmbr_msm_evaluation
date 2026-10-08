###
### Simulate validation datasets
###

### Setwd and clear workspace
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

### Extract arguments (dgm and scenario)
args <- commandArgs(trailingOnly = T)
dgm <- as.numeric(args[1])
scenario <- as.numeric(args[2])
print(paste("dgm = ", dgm))
print(paste("scenario = ", scenario))

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Read in sim_inputs for the DGM
sim_inputs <- readRDS(paste("data/sim_inputs_dgm", dgm, ".rds", sep = ""))

### Set seed
### Fix seed deliberately so all validation datasets share baseline covariates
set.seed(102)

### Set general parameters
n_cohort <- 500000
max_follow <- 10.001
numsteps <- 10001

###
### Define all the input parameters
###

### Shape and scale for transitions
shape12 <- sim_inputs[scenario, "valid_shape12"]
shape13 <- sim_inputs[scenario, "valid_shape13"]
shape14 <- sim_inputs[scenario, "valid_shape14"]
shape23 <- sim_inputs[scenario, "valid_shape23"]
shape24 <- sim_inputs[scenario, "valid_shape24"]
shape34 <- sim_inputs[scenario, "valid_shape34"]

scale12 <- sim_inputs[scenario, "valid_scale12"]
scale13 <- sim_inputs[scenario, "valid_scale13"]
scale14 <- sim_inputs[scenario, "valid_scale14"]
scale23 <- sim_inputs[scenario, "valid_scale23"]
scale24 <- sim_inputs[scenario, "valid_scale24"]
scale34 <- sim_inputs[scenario, "valid_scale34"]

### Coefficients for predictors
beta12_xa <- sim_inputs[scenario, "beta12_xa"]
beta13_xa <- sim_inputs[scenario, "beta13_xa"]
beta14_xa <- sim_inputs[scenario, "beta14_xa"]
beta23_xa <- sim_inputs[scenario, "beta23_xa"]
beta24_xa <- sim_inputs[scenario, "beta24_xa"]
beta34_xa <- sim_inputs[scenario, "beta34_xa"]

beta12_xb <- sim_inputs[scenario, "beta12_xb"]
beta13_xb <- sim_inputs[scenario, "beta13_xb"]
beta14_xb <- sim_inputs[scenario, "beta14_xb"]
beta23_xb <- sim_inputs[scenario, "beta23_xb"]
beta24_xb <- sim_inputs[scenario, "beta24_xb"]
beta34_xb <- sim_inputs[scenario, "beta34_xb"]

beta12_xasq <- sim_inputs[scenario, "beta12_xasq"]
beta13_xasq <- sim_inputs[scenario, "beta13_xasq"]
beta14_xasq <- sim_inputs[scenario, "beta14_xasq"]
beta23_xasq <- sim_inputs[scenario, "beta23_xasq"]
beta24_xasq <- sim_inputs[scenario, "beta24_xasq"]
beta34_xasq <- sim_inputs[scenario, "beta34_xasq"]

beta12_xbsq <- sim_inputs[scenario, "beta12_xbsq"]
beta13_xbsq <- sim_inputs[scenario, "beta13_xbsq"]
beta14_xbsq <- sim_inputs[scenario, "beta14_xbsq"]
beta23_xbsq <- sim_inputs[scenario, "beta23_xbsq"]
beta24_xbsq <- sim_inputs[scenario, "beta24_xbsq"]
beta34_xbsq <- sim_inputs[scenario, "beta34_xbsq"]

### Shape and scale for censoring mechanism
cens_shape <- sim_inputs[scenario, "cens_shape"]
cens_scale <- sim_inputs[scenario, "cens_scale"]
cens_beta_xa <- sim_inputs[scenario, "cens_beta_xa"]
cens_beta_xb <- sim_inputs[scenario, "cens_beta_xb"]

###
### Simulate data
###

### Generate baseline data
### cap xa and xb at 4 to stop extreme values resulting in NA predictions
x_baseline <- data.frame("xa" = rnorm(n_cohort, 0, 1), "xb" = rnorm(n_cohort, 0, 1)) |>
  dplyr::mutate(xa = dplyr::case_when(xa > 4 ~ 4,
                                      xa < -4 ~ -4,
                                      TRUE ~ xa),
                xb = dplyr::case_when(xb > 4 ~ 4,
                                      xb < -4 ~ -4,
                                      TRUE ~ xb),
                xasq = xa^2, xbsq = xb^2)


### Assign function for generating transition data
if (dgm == 1){
  dgm_func <- dgm1
} else if (dgm == 2){
  dgm_func <- dgm2
} else if (dgm == 3){
  dgm_func <- dgm3
} else if (dgm == 4){
  dgm_func <- dgm4
} else {
  stop(paste("Unrecognised dgm:", dgm))
}

### Generate transition data
### Note some arguments are NULL and will not be used (e.g. transition 1 -> 4 for DGM1)
### Implementing it this way means the inputs don't have to be changed manually for each DGM, which could be error prone
df_valid <- dgm_func(n = n_cohort, #number of patients to simulate
                     max_follow = max_follow, #maximum follow up
                     shape12 = shape12, scale12 = scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                     shape13 = shape13, scale13 = scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                     shape14 = shape14, scale14 = scale14, #shape and scale for weibull baseline hazard for transition 1 -> 4
                     shape23 = shape23, scale23 = scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                     shape24 = shape24, scale24 = scale24, #shape and scale for weibull baseline hazard for transition 2 -> 4
                     shape34 = shape34, scale34 = scale34, #shape and scale for weibull baseline hazard for transition 3 -> 4
                     beta12_xa = beta12_xa, beta12_xb = beta12_xb, beta12_xasq = beta12_xasq, beta12_xbsq = beta12_xbsq, #covariate effects for transition 1 -> 2
                     beta13_xa = beta13_xa, beta13_xb = beta13_xb, beta13_xasq = beta13_xasq, beta13_xbsq = beta13_xbsq, #covariate effects for transition 1 -> 3
                     beta14_xa = beta14_xa, beta14_xb = beta14_xb, beta14_xasq = beta14_xasq, beta14_xbsq = beta14_xbsq, #covariate effects for transition 1 -> 4
                     beta23_xa = beta23_xa, beta23_xb = beta23_xb, beta23_xasq = beta23_xasq, beta23_xbsq = beta23_xbsq, #covariate effects for transition 2 -> 3
                     beta24_xa = beta24_xa, beta24_xb = beta24_xb, beta24_xasq = beta24_xasq, beta24_xbsq = beta24_xbsq, #covariate effects for transition 2 -> 4
                     beta34_xa = beta34_xa, beta34_xb = beta34_xb, beta34_xasq = beta34_xasq, beta34_xbsq = beta34_xbsq, #covariate effects for transition 3 -> 4
                     x_in = x_baseline, #baseline predictors, data frame with two columns
                     numsteps = numsteps)


### Convert to mstate format and apply censoring
df_valid_object <- convert_mstate_cens(cohort_in = df_valid,
                                       max_follow = max_follow,
                                       cens_shape = cens_shape,
                                       cens_scale = cens_scale,
                                       cens_beta_xa = cens_beta_xa,
                                       cens_beta_xb = cens_beta_xb,
                                       dgm_in = dgm)

### Extract data_mstate and data_raw
df_valid_mstate <- df_valid_object[["data_mstate"]]
df_valid_raw  <- df_valid_object[["data_raw"]]

### Rename patid to id in df_valid_raw
df_valid_raw <- dplyr::rename(df_valid_raw, id = patid)

### Save relevant objects
saveRDS(df_valid_mstate, paste("data/dgm", dgm, "_df_valid_mstate", scenario, ".rds", sep = ""))
saveRDS(df_valid_raw, paste("data/dgm", dgm, "_df_valid_raw", scenario, ".rds", sep = ""))
print("FINISHED")