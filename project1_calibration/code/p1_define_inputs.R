### Set working directory
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

###
### Define tmat_list
### List of transition matricies for each DGM
###

### Clear workspace
rm(list=ls())

### Create transition matricies
tmat1 <- mstate::transMat(x = list(c(2, 3), c(3), c()),
                          names = paste0("state", 1:3))

tmat2 <- mstate::transMat(x = list(c(2, 4), c(3), c(), c()),
                          names = paste0("state", 1:4))

tmat3 <- mstate::transMat(x = list(c(2, 4), c(3, 4), c(), c()),
                          names = paste0("state", 1:4))

tmat4 <- mstate::transMat(x = list(c(2, 4), c(3, 4), c(4), c()),
                          names = paste0("state", 1:4))


### Create list
tmat_list <- list(tmat1, tmat2, tmat3, tmat4)

### Save list
saveRDS(tmat_list, "data/tmat_list.rds")


###
### Define simulation inputs for DGM1
###

### Clear workspace
rm(list=ls())

### Start by creating the development input parameters for different scenarios
sim_inputs_dgm1 <- data.frame("devel_shape12" = c(rep(1, 4), rep(1.5, 4)), 
                              "devel_shape13" = c(rep(1, 4), rep(1.5, 4)), 
                              "devel_shape23" = c(rep(1, 4), rep(1.5, 4)), 
                              "devel_scale12" = rep(c(10, 5, 10, 10), 2), 
                              "devel_scale13" = rep(c(10, 10, 5, 10), 2), 
                              "devel_scale23" = rep(c(10, 10, 10, 5), 2),
                              "devel_num" = 1:8)

### Add the other input parameters for these scenarios
sim_inputs_dgm1 <- dplyr::mutate(sim_inputs_dgm1,
                                 # Validation datasets shape and scale to be the same by default
                                 valid_shape12 = devel_shape12, 
                                 valid_shape13 = devel_shape13, 
                                 valid_shape23 = devel_shape23, 
                                 valid_scale12 = devel_scale12, 
                                 valid_scale13 = devel_scale13, 
                                 valid_scale23 = devel_scale23,
                                 # Log-hazard ratios are fixed
                                 beta12_xa = 0.5, beta13_xa = 0.5, beta23_xa = 0.5,
                                 beta12_xb = -0.5, beta13_xb = -0.5, beta23_xb = -0.5,
                                 beta12_xasq = 0, beta13_xasq = 0, beta23_xasq = 0,
                                 beta12_xbsq = 0, beta13_xbsq = 0, beta23_xbsq = 0,
                                 # Input parameters for censoring mechanism are fixed
                                 cens_shape = 1,
                                 cens_scale = 50,
                                 cens_beta_xa = 0,
                                 cens_beta_xb = 0)

### Create rows for scenarios where validation dataset scale differs from development dataset
sim_inputs_dgm1 <- dplyr::bind_rows(
  # Original rows
  dplyr::mutate(sim_inputs_dgm1, valid_num = 1),
  # valid_scale12 halved
  dplyr::mutate(sim_inputs_dgm1, valid_scale12 = valid_scale12 / 2, valid_num = 2),
  # valid_scale13 halved
  dplyr::mutate(sim_inputs_dgm1, valid_scale13 = valid_scale13 / 2, valid_num = 3),
  # valid_scale23 halved
  dplyr::mutate(sim_inputs_dgm1, valid_scale23 = valid_scale23 / 2, valid_num = 4),
  # All three halved simultaneously
  dplyr::mutate(sim_inputs_dgm1, valid_scale12 = valid_scale12 / 2,
                valid_scale13 = valid_scale13 / 2,
                valid_scale23 = valid_scale23 / 2,
                valid_num = 5)
) 

### Relocate devel and valid scenario numbers
sim_inputs_dgm1 <- dplyr::relocate(sim_inputs_dgm1, devel_num, valid_num)

### Save the inputs
saveRDS(sim_inputs_dgm1, "data/sim_inputs_dgm1.rds")

###
### Define simulation inputs for DGM2
###

### Clear workspace
rm(list=ls())

### Start by creating the development input parameters for different scenarios
sim_inputs_dgm2 <- data.frame("devel_shape12" = rep(1.5, 2), 
                              "devel_shape14" = rep(1.5, 2), 
                              "devel_shape23" = rep(1.5, 2), 
                              "devel_scale12" = c(10, 20), 
                              "devel_scale14" = c(10, 5), 
                              "devel_scale23" = c(10, 10),
                              "devel_num" = 1:2)

### Add the other input parameters for these scenarios
sim_inputs_dgm2 <- dplyr::mutate(sim_inputs_dgm2,
                                 # Validation datasets shape and scale to be the same by default
                                 valid_shape12 = devel_shape12, 
                                 valid_shape14 = devel_shape14, 
                                 valid_shape23 = devel_shape23, 
                                 valid_scale12 = devel_scale12, 
                                 valid_scale14 = devel_scale14, 
                                 valid_scale23 = devel_scale23,
                                 # Log-hazard ratios are fixed
                                 beta12_xa = 0.5, beta14_xa = 0.5, beta23_xa = 0.5,
                                 beta12_xb = -0.5, beta14_xb = -0.5, beta23_xb = -0.5,
                                 beta12_xasq = 0, beta14_xasq = 0, beta23_xasq = 0,
                                 beta12_xbsq = 0, beta14_xbsq = 0, beta23_xbsq = 0,
                                 # Input parameters for censoring mechanism are fixed
                                 cens_shape = 1,
                                 cens_scale = 50,
                                 cens_beta_xa = 0,
                                 cens_beta_xb = 0)

### Create rows for scenarios where validation dataset scale differs from development dataset
### Do not include rows for where validation dataset is same as devel
sim_inputs_dgm2 <- dplyr::bind_rows(
  # valid_scale12 halved
  dplyr::mutate(sim_inputs_dgm2, valid_scale12 = valid_scale12 / 2, valid_num = 1),
  # valid_scale14 halved
  dplyr::mutate(sim_inputs_dgm2, valid_scale14 = valid_scale14 / 2, valid_num = 2),
  # valid_scale23 halved
  dplyr::mutate(sim_inputs_dgm2, valid_scale23 = valid_scale23 / 2, valid_num = 3)
) 

### Relocate devel and valid scenario numbers
sim_inputs_dgm2 <- dplyr::relocate(sim_inputs_dgm2, devel_num, valid_num)

### Save the inputs
saveRDS(sim_inputs_dgm2, "data/sim_inputs_dgm2.rds")

###
### Define simulation inputs for DGM3
###

### Clear workspace
rm(list=ls())

### Start by creating the development input parameters for different scenarios
sim_inputs_dgm3 <- data.frame("devel_shape12" = rep(1.5, 2), 
                              "devel_shape14" = rep(1.5, 2), 
                              "devel_shape23" = rep(1.5, 2), 
                              "devel_shape24" = rep(1.5, 2), 
                              "devel_scale12" = c(10, 20), 
                              "devel_scale14" = c(10, 5), 
                              "devel_scale23" = c(10, 10),
                              "devel_scale24" = c(10, 10),
                              "devel_num" = 1:2)

### Add the other input parameters for these scenarios
sim_inputs_dgm3 <- dplyr::mutate(sim_inputs_dgm3,
                                 # Validation datasets shape and scale to be the same by default
                                 valid_shape12 = devel_shape12, 
                                 valid_shape14 = devel_shape14, 
                                 valid_shape23 = devel_shape23, 
                                 valid_shape24 = devel_shape24, 
                                 valid_scale12 = devel_scale12, 
                                 valid_scale14 = devel_scale14, 
                                 valid_scale23 = devel_scale23,
                                 valid_scale24 = devel_scale24,
                                 # Log-hazard ratios are fixed
                                 beta12_xa = 0.5, beta14_xa = 0.5, beta23_xa = 0.5, beta24_xa = 0.5,
                                 beta12_xb = -0.5, beta14_xb = -0.5, beta23_xb = -0.5, beta24_xb = -0.5,
                                 beta12_xasq = 0, beta14_xasq = 0, beta23_xasq = 0, beta24_xasq = 0,
                                 beta12_xbsq = 0, beta14_xbsq = 0, beta23_xbsq = 0, beta24_xbsq = 0,
                                 # Input parameters for censoring mechanism are fixed
                                 cens_shape = 1,
                                 cens_scale = 50,
                                 cens_beta_xa = 0,
                                 cens_beta_xb = 0)

### Create rows for scenarios where validation dataset scale differs from development dataset
### Do not include rows for where validation dataset is same as devel
sim_inputs_dgm3 <- dplyr::bind_rows(
  # valid_scale12 halved
  dplyr::mutate(sim_inputs_dgm3, valid_scale12 = valid_scale12 / 2, valid_num = 1),
  # valid_scale14 halved
  dplyr::mutate(sim_inputs_dgm3, valid_scale14 = valid_scale14 / 2, valid_num = 2),
  # valid_scale23 halved
  dplyr::mutate(sim_inputs_dgm3, valid_scale23 = valid_scale23 / 2, valid_num = 3),
  # valid_scale24 halved
  dplyr::mutate(sim_inputs_dgm3, valid_scale24 = valid_scale24 / 2, valid_num = 4)
) 

### Relocate devel and valid scenario numbers
sim_inputs_dgm3 <- dplyr::relocate(sim_inputs_dgm3, devel_num, valid_num)

### Save inputs
saveRDS(sim_inputs_dgm3, "data/sim_inputs_dgm3.rds")

###
### Define simulation inputs for DGM4
###

### Clear workspace
rm(list=ls())

### Start by creating the development input parameters for different scenarios
sim_inputs_dgm4 <- data.frame("devel_shape12" = rep(1.5, 2), 
                              "devel_shape14" = rep(1.5, 2), 
                              "devel_shape23" = rep(1.5, 2), 
                              "devel_shape24" = rep(1.5, 2), 
                              "devel_shape34" = rep(1.5, 2), 
                              "devel_scale12" = c(10, 20), 
                              "devel_scale14" = c(10, 5), 
                              "devel_scale23" = c(10, 10),
                              "devel_scale24" = c(10, 10),
                              "devel_scale34" = c(10, 10),
                              "devel_num" = 1:2)

### Add the other input parameters for these scenarios
sim_inputs_dgm4 <- dplyr::mutate(sim_inputs_dgm4,
                                 # Validation datasets shape and scale to be the same by default
                                 valid_shape12 = devel_shape12, 
                                 valid_shape14 = devel_shape14, 
                                 valid_shape23 = devel_shape23, 
                                 valid_shape24 = devel_shape24, 
                                 valid_shape34 = devel_shape34,
                                 valid_scale12 = devel_scale12, 
                                 valid_scale14 = devel_scale14, 
                                 valid_scale23 = devel_scale23,
                                 valid_scale24 = devel_scale24,
                                 valid_scale34 = devel_scale34,
                                 # Log-hazard ratios are fixed
                                 beta12_xa = 0.5, beta14_xa = 0.5, beta23_xa = 0.5, beta24_xa = 0.5, beta34_xa = 0.5,
                                 beta12_xb = -0.5, beta14_xb = -0.5, beta23_xb = -0.5, beta24_xb = -0.5, beta34_xb = -0.5,
                                 beta12_xasq = 0, beta14_xasq = 0, beta23_xasq = 0, beta24_xasq = 0, beta34_xasq = 0,
                                 beta12_xbsq = 0, beta14_xbsq = 0, beta23_xbsq = 0, beta24_xbsq = 0, beta34_xbsq = 0,
                                 # Input parameters for censoring mechanism are fixed
                                 cens_shape = 1,
                                 cens_scale = 50,
                                 cens_beta_xa = 0,
                                 cens_beta_xb = 0)

### Create rows for scenarios where validation dataset scale differs from development dataset
### Do not include rows for where validation dataset is same as devel
sim_inputs_dgm4 <- dplyr::bind_rows(
  # valid_scale12 halved
  dplyr::mutate(sim_inputs_dgm4, valid_scale12 = valid_scale12 / 2, valid_num = 1),
  # valid_scale14 halved
  dplyr::mutate(sim_inputs_dgm4, valid_scale14 = valid_scale14 / 2, valid_num = 2),
  # valid_scale23 halved
  dplyr::mutate(sim_inputs_dgm4, valid_scale23 = valid_scale23 / 2, valid_num = 3),
  # valid_scale24 halved
  dplyr::mutate(sim_inputs_dgm4, valid_scale24 = valid_scale24 / 2, valid_num = 4),
  # valid_scale34 halved
  dplyr::mutate(sim_inputs_dgm4, valid_scale34 = valid_scale34 / 2, valid_num = 5)
) 

### Relocate devel and valid scenario numbers
sim_inputs_dgm4 <- dplyr::relocate(sim_inputs_dgm4, devel_num, valid_num)

### Save inputs
saveRDS(sim_inputs_dgm4, "data/sim_inputs_dgm4.rds")