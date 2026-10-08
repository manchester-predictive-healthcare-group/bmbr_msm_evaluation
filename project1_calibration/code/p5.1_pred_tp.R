###
### Program to estimate transition probabilities
###

### Setwd and clear workspace
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

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
# msm_cox_fit <- readRDS(paste("data/dgm", dgm, "_msm_cox_fit", devel_num, ".rds", sep = ""))
msm_flexsurv_fit <- readRDS(paste("data/dgm", dgm, "_msm_flexsurv_fit", devel_num, ".rds", sep = ""))

# ### Read in development data (necessary or predict function won't work)
# df_devel_mstate <- readRDS(paste("data/dgm", dgm, "_df_devel_mstate", devel_num, ".rds", sep = ""))
# df_devel_raw <- readRDS(paste("data/dgm", dgm, "_df_devel_raw", devel_num, ".rds", sep = ""))

### Read in validation data
### Note the values of xa and xb are the same for every validation dataset, meaning so are the predicted risks.
### We can therefore use the same validation dataset each time, and only need to estimate risks once per development dataset (and developed model)
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
                        m_iter = 10000)


### Save them
saveRDS(probtrans_out, paste("data/dgm", dgm, "_pt_scenario", devel_num, "_set", set, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))
