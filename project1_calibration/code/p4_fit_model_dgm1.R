###
### Fit multistate models
###

### Setwd and clear workspace
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

### Define dgm
dgm <- 1

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

###
### Fit multistate models
###

### Loop through development datasets
for (devel_num in 1:8){
  
  print(paste("devel_num = ", devel_num, Sys.time()))
  
  ### Read in development data
  df_devel_mstate <- readRDS(paste("data/dgm", dgm, "_df_devel_mstate", devel_num, ".rds", sep = ""))
  df_devel_raw <- readRDS(paste("data/dgm", dgm, "_df_devel_raw", devel_num, ".rds", sep = ""))
  
  ### Create list to store
  flexsurv_model_list <- vector("list", 3)
  
  ### For development numbers 1 - 4, which are exponential
  if (devel_num %in% 1:4){
    for (i in 1:3){
      print(paste("transition = ", i, sep = ""))
      flexsurv_model_list[[i]] <- flexsurv::flexsurvreg(Surv(time, status) ~ xa + xb, 
                                                        subset = (trans == i), 
                                                        data = df_devel_mstate, 
                                                        dist = "exp")
    }
  } else if (devel_num %in% 5:8){
    ### For development numbers 5 - 8, which are weibull
    for (i in 1:3){
      print(paste("transition = ", i, sep = ""))
      flexsurv_model_list[[i]] <- flexsurv::flexsurvreg(Surv(time, status) ~ xa + xb, 
                                                        subset = (trans == i), 
                                                        data = df_devel_mstate, 
                                                        dist = "weibullPH")
    }
  }
  
  ### Save
  saveRDS(flexsurv_model_list, paste("data/dgm", dgm, "_msm_flexsurv_fit", devel_num, ".rds", sep = ""))
  print("FINISHED")
  
}
