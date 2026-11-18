###
### Apply multistate calibration functions to the EBMT data
###
### This script fits a multistate model using parametric survival models for each transition. 
### Calibration is then assessed in the same cohort, evaluating calibration of the transition probabilities,
### and the cause-specific failure functions.
###

###
### Prelim
###

### Clear workspace
rm(list = ls())

### Set working directory here
### Or use packages such as 'here' package
setwd("PATH/TO/ROOT")

### Load packages
library(survival)
library(dplyr)
library(mstate)
library(ggplot2)
library(calibmsm)
library(patchwork)

### Directory in which to save results
output_dir <- "ebmt_example"

### Time point at which to assess calibration (five years, in days)
t_in <- 5 * 365.25

### Number of simulations used to estimate transition probabilities
m_iter <- 10000

### Set seed so the transition-probability simulations are reproducible
set.seed(101)

### Create output directory
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

###
### Load functions
### These functions are available on the following GitHub page: 
### https://github.com/manchester-predictive-healthcare-group/bmbr_msm_evaluation/tree/main/project1_calibration 
### File: R/z_functions.R
### However they are pasted here so the script remains standalone.
###

###
### Write a function to get predicted transition probabilities for an individual, when model was developed using flexible parametric
###
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

###
### Assessing calibration using proportional hazards regression approach 
### (graphical calibration curves, Austin et al, DOI: 10.1002/sim.8570)
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
### Function to estimate calibration of the cause-specific hazards functions.
### This is specific to a multistate model developed using "flexsurv"
### It will generate predictions from fit_msm for individuals in the validation cohort.
### then apply est_calib_csh_ph to calculate calibration curves, ICI, E50 and E90 for each transition.
### Can be applied to a model with any number of transitions
###
calib_csh_flexsurv <- function(data_mstate, fit_msm, t, nk = 4){
  
  ### Get number of transitions
  transitions_n <- max(data_mstate$trans)
  
  ### Create plots list object
  calib_list <- vector("list", transitions_n)
  
  ###
  ### Get names of transitions
  ###
  
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
    
    ### Restrict the multistate data to the current transition for calibration-model fitting
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
### Note that we don't re-use make_csh_panel from the simuation (like we did with make_tp_panel). Instead this is written directly below 
### because the study function make_csh_panel reads the simulation transition matrix, whereas this example uses the 
### EBMT matrix. tmat is used to create the ggplot titles for each panel.
###

###
### Load and prepare the EBMT data
###

### Load the one-row-per-person and mstate-format datasets
utils::data("ebmtcal", package = "calibmsm")
utils::data("msebmtcal", package = "calibmsm")

### Copy the datasets so the package objects are not modified
data_raw <- ebmtcal
data_mstate <- msebmtcal

### Define the transition matrix used for the EBMT model
### States: Tx, Rec, AE, Rec+AE, Rel, Death
tmat <- attr(data_mstate, "trans")

### Attach the transition matrix because the study functions use this attribute
attr(data_mstate, "trans") <- tmat

### Add the baseline predictors to the mstate-format data
data_mstate <- dplyr::left_join(
  data_mstate,
  dplyr::select(data_raw, id, year, agecl, proph, match),
  by = "id"
)

### Remove individuals with missing values in predictors used by the model
complete_ids <- data_raw$id[stats::complete.cases(data_raw[, c("year", "agecl", "proph", "match")])]
data_raw <- data_raw[data_raw$id %in% complete_ids, ]
data_mstate <- data_mstate[data_mstate$id %in% complete_ids, ]

### Use complete case dataset for development and validation
devel_ids <- complete_ids
valid_ids <- complete_ids

data_devel_raw <- data_raw[data_raw$id %in% devel_ids, ]
data_valid_raw <- data_raw[data_raw$id %in% valid_ids, ]
data_devel_mstate <- data_mstate[data_mstate$id %in% devel_ids, ]
data_valid_mstate <- data_mstate[data_mstate$id %in% valid_ids, ]

### Re-number validation IDs because est_risk_flexsurv selects rows by row number
data_valid_raw$id_original <- data_valid_raw$id
data_valid_raw$id <- seq_len(nrow(data_valid_raw))

### Apply the same new IDs to the mstate-format validation data
id_key <- dplyr::select(data_valid_raw, id_original, id)
data_valid_mstate <- dplyr::left_join(
  data_valid_mstate,
  id_key,
  by = c("id" = "id_original"),
  suffix = c("_original", "")
)

### Restore the transition matrix after subsetting and joining
attr(data_devel_mstate, "trans") <- tmat
attr(data_valid_mstate, "trans") <- tmat

###
### Fit the transition-specific models
###

### Create list to store one fitted model for each transition
transitions_n <- max(data_devel_mstate$trans)
msm_flexsurv_fit <- vector("list", transitions_n)

### Fit a Weibull proportional hazards model for each transition
for (transition in seq_len(transitions_n)){
  print(paste("Fitting transition", transition, "of", transitions_n))

  msm_flexsurv_fit[[transition]] <- flexsurv::flexsurvreg(
    survival::Surv(time, status) ~ year + agecl + proph + match,
    subset = (trans == transition),
    data = data_devel_mstate,
    dist = "weibullPH"
  )
}

### Save fitted multistate model
saveRDS(msm_flexsurv_fit, file.path(output_dir, "ebmt_flexsurv_fit.rds"))

###
### Estimate transition probabilities in the validation cohort
###

### Estimate probabilities of occupying each state at five years
### The outer lapply returns one list element per individual, with the
### prediction function returning one nested element per prediction time
tp <- lapply(
  X = seq_len(nrow(data_valid_raw)),
  FUN = est_risk_flexsurv,
  data_raw = data_valid_raw,
  tmat = tmat,
  msm_fit = msm_flexsurv_fit,
  t_vec = t_in,
  m_iter = m_iter
)

### Combine the individual predictions into one data frame
tp <- do.call(rbind, lapply(tp, function(x){
  do.call(rbind, x)
}))

### Keep probabilities away from exactly zero and one
tp <- dplyr::mutate(
  tp,
  dplyr::across(
    .cols = dplyr::starts_with("pstate"),
    .fns = ~ pmin(pmax(.x, 0.000001), 0.999999)
  )
)

### Extract the state-occupation probabilities required by calib_msm
tp_pred <- dplyr::select(tp, dplyr::starts_with("pstate"))

### Save predicted transition probabilities
saveRDS(tp, file.path(output_dir, "ebmt_transition_probabilities.rds"))

###
### Assess calibration of transition probabilities
###

### Estimate flexible calibration curves for each state
calib_tp_obj <- calibmsm::calib_msm(
  data_ms = data_valid_mstate,
  data_raw = data_valid_raw,
  j = 1,
  s = 0,
  t = t_in,
  tp_pred = tp_pred,
  w_covs = c("year", "agecl", "proph", "match"),
  w_max_follow = t_in,
  calib_type = "blr",
  curve_type = "rcs",
  rcs_nk = 3,
  CI = FALSE,
  assess_moderate = TRUE,
  assess_mean = TRUE
)

### Save transition-probability calibration results
saveRDS(calib_tp_obj, file.path(output_dir, "ebmt_calibration_tp.rds"))

###
### Assess calibration of cause-specific failure functions
###

### Apply the function implemented in the study
calib_csh_full <- calib_csh_flexsurv(
  data_mstate = data_valid_mstate,
  fit_msm = msm_flexsurv_fit,
  t = t_in,
  nk = 3
)

### Retain the transition-specific calibration objects
calib_csh_obj <- calib_csh_full[["calib_list"]]

### Save cause-specific calibration results
saveRDS(calib_csh_obj, file.path(output_dir, "ebmt_calibration_csh.rds"))

###
### Create combined calibration figure in the study format
###

### Create transition-probability panels
### plotdata is generated by calib_msm and used by the study plotting functions
tp_panels <- lapply(
  seq_along(calib_tp_obj$plotdata),
  function(i){
    make_tp_panel(
      plotdata = calib_tp_obj$plotdata[[i]],
      state_num = colnames(tmat)[i]
    )
  }
)

### Create cause-specific failure-function panels
### This is written directly here because the study function make_csh_panel reads the simulation
### transition matrix, whereas this example uses the EBMT matrix. tmat is used to create the ggplot titles
### for each panel.
csh_panels <- vector("list", length(calib_csh_obj))

for (transition in seq_along(calib_csh_obj)){
  ### Extract the calibration data stored in the ggplot object
  plotdata <- calib_csh_obj[[transition]][["plot"]][["data"]]

  ### Calculate calibration metrics
  metrics <- calc_metrics(plotdata$pred, plotdata$pred.obs)

  ### Thin only the calibration curve
  curve_dat <- plotdata |>
    dplyr::arrange(pred)

  ### Get the transition name
  transition_location <- which(tmat == transition, arr.ind = TRUE)
  transition_name <- paste0(
    rownames(tmat)[transition_location[, 1]],
    " -> ",
    colnames(tmat)[transition_location[, 2]]
  )

  ### Get the upper plotting limit and margin for the calibration metrics
  xmax <- max(curve_dat$pred, na.rm = TRUE)
  metric_margin <- 0.025 * xmax

  ### Create density plot using the same format as the study figures
  hist_plot <- make_histogram(plotdata) +
    ggplot2::ggtitle(paste("Transition:", transition_name)) +
    ggplot2::xlim(c(0, xmax))

  ### Create calibration plot using the same format as the study figures
  calib_plot <- ggplot2::ggplot(curve_dat) +
    ggplot2::geom_line(
      ggplot2::aes(x = pred, y = pred.obs),
      colour = "red"
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::coord_fixed(
      xlim = c(0, xmax),
      ylim = c(0, xmax)
    ) +
    ggplot2::annotate(
      "text",
      x = xmax - metric_margin,
      y = metric_margin,
      label = paste0(
        "ICI: ", round(metrics[["ICI"]], 3), "\n",
        "E50: ", round(metrics[["E50"]], 3), "\n",
        "E90: ", round(metrics[["E90"]], 3)
      ),
      hjust = 1,
      vjust = 0,
      size = 3
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(t = 0)
    ) +
    ggplot2::xlab("Predicted risk") +
    ggplot2::ylab("Observed risk")

  ### Combine density and calibration plots into one panel
  csh_panels[[transition]] <- hist_plot /
    calib_plot +
    patchwork::plot_layout(heights = c(1, 4))
}

### Create row titles used in the main study figures
tp_title <- make_title("Calibration of the Transition Probabilities")
csh_title <- make_title("Calibration of the Cause-Specific Failure Functions")

### Arrange the state panels over two rows
### Six states fit naturally into two rows of three panels
tp_plot <- patchwork::wrap_plots(tp_panels, ncol = 3)

### Arrange the eleven transition panels over four rows
csh_plot <- patchwork::wrap_plots(csh_panels, ncol = 3)

### Combine titles and panels using the same relative title height as the study
final_plot <- (tp_title / tp_plot / csh_title / csh_plot) +
  patchwork::plot_layout(heights = c(0.05, 2, 0.05, 4))

### Save the combined figure
ragg::agg_png(
  filename = file.path(output_dir, "ebmt_combined_calibration_plot.png"),
  width = 10.5,
  height = 26.75,
  units = "in",
  res = 300
)
print(final_plot)
grDevices::dev.off()

### Save the complete combined plot object so it can be edited without re-fitting models
saveRDS(final_plot, file.path(output_dir, "ebmt_combined_calibration_plot.rds"))

print("FINISHED")
