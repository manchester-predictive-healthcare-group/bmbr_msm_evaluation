###
### Create combined calibration figures
### Transition probabilities + cause-specific hazards
###

### Clear workspace and setwd
rm(list = ls())
setwd("/mnt/bmh01-rds/Pate_bmbr/project1")

### Source functions and prelim
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

###
### User settings
###

### Time point to plot
t_in <- 5

### Number of points to thin to in plots
thin_n = 5000

###
### Number of scenarios for each DGM
###
scenario_max <- c(
  "1" = 40,
  "2" = 6,
  "3" = 8,
  "4" = 10
)

###
### Run all scenarios
###
for (dgm in 1:4){

  for (scenario in seq_len(scenario_max[as.character(dgm)])){

    print(paste("DGM",dgm,"scenario",scenario))

    create_combined_calibration_plot(
      dgm = dgm,
      scenario = scenario,
      t_in = t_in,
      model_updated = FALSE,
      thin_n = thin_n
    )

  }

}

### Testers
# scenario <- 1
# for (dgm in 2:2){
# 
#     print(paste("DGM",dgm,"scenario",scenario))
# 
#     create_combined_calibration_plot(
#       dgm = dgm,
#       scenario = scenario,
#       t_in = t_in,
#       model_updated = FALSE
#     )
# 
#   }