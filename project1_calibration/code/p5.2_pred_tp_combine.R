###
### Combine transition probabilities estimated in section 5.1
###

### Setwd and clear workspace
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

### Write a function to extract and save for a specific development dataset
combine_tp <- function(devel_num){
  
  ## print devel_nums
  print(paste(devel_num, Sys.time()))
  
  ## Read in
  tp <- lapply(1:100, function(x){readRDS(paste("data/dgm", dgm, "_pt_scenario", devel_num, "_set", x, ".rds", sep = ""))})
  
  ## Combine
  tp <- 
    do.call(rbind, lapply(tp, function(y){
      do.call(rbind, lapply(y, function(x){
        do.call(rbind, x)
      }))
    })
    )
  
  ### If any values are equal to 0 and 1, change
  tp <- dplyr::mutate(tp, dplyr::across(
    .cols = dplyr::starts_with("pstate"),
    .fns = ~ pmin(pmax(.x, 0.000001), 0.999999)
  ))
  
  ## Save
  saveRDS(tp, paste("data/dgm", dgm, "_pt_scenario", devel_num, "_all.rds", sep = ""))
}

### Run this function for relevant scenario in each DGM 
### (8 development dataset for DGM1, 2 development datasets for the others)
for (dgm in 1:4){
  print(paste("dgm = ", dgm, Sys.time()))
  if (dgm == 1){
    lapply(1:8, combine_tp)
  } else {
    lapply(1:2, combine_tp)
  }
}
