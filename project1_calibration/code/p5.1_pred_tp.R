###
### Program to estimate transition probabilities
###

### Setwd and clear workspace
rm(list = ls())
setwd("PATH/TO/ROOT")

### Extract arguments
args <- commandArgs(trailingOnly = T)
dgm <- as.numeric(args[1])
devel_num <- as.numeric(args[2])
set_size <- as.numeric(args[3])
set <- as.numeric(args[4])
print(paste("dgm = ", dgm))
print(paste("devel_num = ", devel_num))
print(paste("set = ", set))

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Read in model fit
msm_flexsurv_fit <- readRDS(paste("data/dgm", dgm, "_msm_flexsurv_fit", devel_num, ".rds", sep = ""))

### Read in validation data
### Baseline predictor values are identical across validation scenarios because the same seed is used (predictors stay same, outcomes change)
### Predicted risks therefore depend only on the development model, so they are calculated once per development dataset.
df_valid_mstate <- readRDS(paste("data/dgm", dgm, "_df_valid_mstate", 1, ".rds", sep = ""))
df_valid_raw <- readRDS(paste("data/dgm", dgm, "_df_valid_raw", 1, ".rds", sep = ""))

###
### Generate transition probs for the set/batch of individuals defined by the variable "set"
###

### Define tmat
tmat <- attributes(df_valid_mstate)$trans

### Calculate transition probabilities
probtrans_out <- lapply(FUN = est_risk_flexsurv, ((set-1)*set_size + 1):(set*set_size), 
                        data_raw = df_valid_raw, 
                        tmat = tmat, 
                        msm_fit = msm_flexsurv_fit,
                        t_vec = c(2.5, 5, 7.5),
                        m_iter = 50000)

### Save them
saveRDS(probtrans_out, paste("data/dgm", dgm, "_pt_scenario", devel_num, "_set", set, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))
