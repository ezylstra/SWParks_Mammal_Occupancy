################################################################################
# Functions to be used in any script in SWParks_Mammal_Occupancy repo

# ER Zylstra
# Updated 2023-02-02
################################################################################

# List of functions included:
  # create_cov_list()
  # vegclass_estimates()
  # marginal_plot_occ()
  # marginal_plot_det()
  # paNA() 
  # propNA() 

#------------------------------------------------------------------------------#
# paNA: Aggregate daily detection data during each multi-day sampling occasion
#------------------------------------------------------------------------------#

# INPUTS
# x: vector with daily detection data during one sampling occasion for one site

# RETURNS
# Single value: 
  # NA if camera wasn't operational throughout entire occasion (all values = NA)
  # 1 if species was detected one or more times (even if there are NAs)
  # 0 if species was never detected

paNA <- function(x) {
  if (sum(is.na(x)) == length(x)) {NA} else 
    if (sum(x, na.rm = TRUE) == 0) {0} else {1} 
}

#------------------------------------------------------------------------------#
# propNA: Calculate the proportion of a sampling occasion camera was operational
#------------------------------------------------------------------------------#

# INPUTS
# x: vector with daily detection data during one sampling occasion for one site

# RETURNS
# Single value in [0, 1]

# Create a function to calculate the proportion of a sampling period a camera
# was operational
propNA <- function(x) {
  (occasions$duration[1] - sum(is.na(x))) / occasions$duration[1]
}

#------------------------------------------------------------------------------#
# create_cov_list: Convert a model formula to a vector with covariate names
#------------------------------------------------------------------------------#

# INPUTS
# object: a character string representing a formula for occurrence or 
  # detection probability in an occupancy model (should begin with a ~)

# RETURNS
# z_list: a character vector that contains all covariates in a model formula 
  # (or returns a "1" if there are no covariates)

create_cov_list <- function(object) {
  if (!methods::is(object, "character")) {
    stop("Object should be a character.")
  }

  z_list <- object %>%
    str_remove(pattern = "~ ") %>% 
    str_remove_all(pattern = "I[(]") %>%
    str_remove_all(pattern = "[)]") %>%
    str_remove_all(pattern = "\\^2") %>%
    str_split_1(pattern = " [+] ")

    return(z_list)
}
  
#------------------------------------------------------------------------------#
# vegclass_estimates: Create a table with occupancy or detection probabilities 
# in each vegetation class
#------------------------------------------------------------------------------#

# INPUTS
# model: output from spOccupancy single-season model
# parameter: character indicating whether to calculate occupancy or detection 
  # probabilities (note: vegclasses must have been included as a covariate
  # in the model formula for that parameter)
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)

# RETURNS
# vegclass_table: a table with mean, SD, and 95% CI for occupancy/detection
  # probabilities in each vegetation class

vegclass_estimates <- function(model, 
                               parameter = c("occ", "det"),
                               lower_ci = 0.025,
                               upper_ci = 0.975) {
  
  parameter <- match.arg(arg = parameter)
  
  # Create table to hold results
  vegclass_table <- data.frame(vegclass = 1:3,
                               mean_prob = NA,
                               sd_prob = NA,
                               ci_lower = NA,
                               ci_upper = NA)
  
  if (parameter == "occ") {
    samples <- model$beta.samples
    submodel <- "occupancy"
  } else {
    samples <- model$alpha.samples
    submodel <- "detection"
  }
  
  if (sum(str_detect(colnames(samples), "vegclass")) == 0) {
    stop("vegclasses must be included in model for ", submodel)
  }
  
  # Probability of occupancy/detection in vegclass1 (reference level)
  vegclass1 <- exp(samples[,"(Intercept)"])/(1 + exp(samples[,"(Intercept)"])) 
  vegclass_table$mean_prob[1] <- mean(vegclass1)
  vegclass_table$sd_prob[1] <- sd(vegclass1)
  vegclass_table$ci_lower[1] <- quantile(vegclass1, lower_ci)
  vegclass_table$ci_upper[1] <- quantile(vegclass1, upper_ci)
  
  # Probability of occupancy/detection in vegclass2
  vegclass2 <- samples[,"(Intercept)"] + samples[,"vegclass2"]
  vegclass2 <- exp(vegclass2)/(1 + exp(vegclass2)) 
  vegclass_table$mean_prob[2] <- mean(vegclass2)
  vegclass_table$sd_prob[2] <- sd(vegclass2)
  vegclass_table$ci_lower[2] <- quantile(vegclass2, lower_ci)
  vegclass_table$ci_upper[2] <- quantile(vegclass2, upper_ci)
  
  # Probability of occupancy/detectin in vegclass3
  vegclass3 <- samples[,"(Intercept)"] + samples[,"vegclass3"]
  vegclass3 <- exp(vegclass3)/(1 + exp(vegclass3)) 
  vegclass_table$mean_prob[3] <- mean(vegclass3)
  vegclass_table$sd_prob[3] <- sd(vegclass3)
  vegclass_table$ci_lower[3] <- quantile(vegclass3, lower_ci)
  vegclass_table$ci_upper[3] <- quantile(vegclass3, upper_ci)
  
  return(vegclass_table)
}

#------------------------------------------------------------------------------#
# marginal_plot_occ: Create figure depicting marginal effect of a covariate on 
# occurrence probability (predicted covariate effect assuming all other 
# covariates held constant)
#------------------------------------------------------------------------------#

# INPUTS
# covariate: character name of covariate
# model: output from spOccupancy single-season model
# data_list: list of data required to run model in spOccupancy 
# covariate_table: table with information and axis labels for covariates
# central_meas: metric for summarizing predictions (mean [default] or median)
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)
# line_color: color for line depicting mean/median predicted probability
# transparency: alpha value specifying transparency of shaded CI

# RETURNS
# marg_plot: ggplot object depicting marginal effect of covariate on orig scale

marginal_plot_occ <- function(covariate, 
                              model, 
                              data_list,
                              covariate_table,
                              central_meas = c(mean, median),
                              lower_ci = 0.025,
                              upper_ci = 0.975,
                              line_color = "forestgreen",
                              transparency = 0.2) {
  
  cols <- str_subset(colnames(best$beta.samples), pattern = covariate)
  beta_samples <- best$beta.samples[,c("(Intercept)", cols)]
  X_cov <- seq(from = min(data_list$occ.covs[, covariate]), 
               to = max(data_list$occ.covs[, covariate]),
               length = 100)
  X_cov <- cbind(1, X_cov)
  
  # If there are quadratic effects, add column in X_cov
  if (ncol(beta_samples) == 3) {
    X_cov <- cbind(X_cov, X_cov[,2]^2)
  } 
  
  preds <-  X_cov %*% t(beta_samples)
  preds <- exp(preds)/(1 + exp(preds))
  preds_cent <- apply(preds, 1, central_meas)
  preds_lcl <- apply(preds, 1, quantile, lower_ci)
  preds_ucl <- apply(preds, 1, quantile, upper_ci)
  
  # Identify covariate values for x-axis (on original scale)
  if (str_detect(covariate, "_z")) {
    cov_mn <- mean(data_list$occ.covs[,str_remove(covariate, "_z")])
    cov_sd <- sd(data_list$occ.covs[,str_remove(covariate, "_z")])
    cov_plot <- X_cov[,2] * cov_sd + cov_mn
  } else {
    cov_plot <- X_cov[,2] 
  }
  
  # Create and save plots for later viewing
  data_plot <- data.frame(x = cov_plot,
                          cent = preds_cent,
                          lcl = preds_lcl,
                          ucl = preds_ucl)
  cov_name <- str_remove(covariate, "_z")
  
  cred_interval <- (upper_ci - lower_ci) * 100
  yaxis_label <- paste0("Predicted occurrence probability (", 
                        cred_interval, 
                        "% CI)")
  
  marg_plot <- ggplot(data = data_plot, aes(x = x)) + 
    geom_line(aes(y = cent), col = line_color) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = transparency) +
    labs(x = covariates$axis_label[covariates$short_name == cov_name],
         y = yaxis_label) +
    theme_classic()
  
  return(marg_plot)
}

#------------------------------------------------------------------------------#
# marginal_plot_det: Create figure depicting marginal effect of a covariate on 
# detection probability (predicted covariate effect assuming all other 
# covariates held constant)
#------------------------------------------------------------------------------#

# INPUTS
# covariate: character name of covariate
# model: output from spOccupancy single-season model
# data_list: list of data required to run model in spOccupancy 
# covariate_table: table with information and axis labels for covariates
# central_meas: metric for summarizing predictions (mean [default] or median)
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)
# line_color: color for line depicting mean/median predicted probability
# transparency: alpha value specifying transparency of shaded CI

# RETURNS
# marg_plot: ggplot object depicting marginal effect of covariate on orig scale

marginal_plot_det <- function(covariate, 
                              model, 
                              data_list,
                              covariate_table,
                              central_meas = c(mean, median),
                              lower_ci = 0.025,
                              upper_ci = 0.975,
                              line_color = "forestgreen",
                              transparency = 0.2) {
  
  cols <- str_subset(colnames(best$alpha.samples), pattern = covariate)
  alpha_samples <- best$alpha.samples[,c("(Intercept)", cols)]
  X_cov <- seq(from = min(data_list$det.covs[[covariate]], na.rm = TRUE), 
               to = max(data_list$det.covs[[covariate]], na.rm = TRUE),
               length = 100)
  X_cov <- cbind(1, X_cov)
  
  # If there are quadratic effects, add column in X_cov
  if (ncol(alpha_samples) == 3) {
    X_cov <- cbind(X_cov, X_cov[,2]^2)
  } 
  
  preds <-  X_cov %*% t(alpha_samples)
  preds <- exp(preds)/(1 + exp(preds))
  preds_cent <- apply(preds, 1, central_meas)
  preds_lcl <- apply(preds, 1, quantile, lower_ci)
  preds_ucl <- apply(preds, 1, quantile, upper_ci)
  
  # Identify covariate values for x-axis (on original scale)
  if (str_detect(covariate, "_z")) {
    if (str_remove(covariate, "_z") %in% colnames(data_list$occ.covs)) {
      cov_mn <- mean(data_list$occ.covs[,str_remove(covariate, "_z")])
      cov_sd <- sd(data_list$occ.covs[,str_remove(covariate, "_z")])
    } else {
      cov_mn <- mean(data_list$det.covs[[str_remove(covariate, "_z")]], na.rm = TRUE)
      cov_sd <- sd(data_list$det.covs[[str_remove(covariate, "_z")]], na.rm = TRUE)
    }
    cov_plot <- X_cov[,2] * cov_sd + cov_mn
  } else {
    cov_plot <- X_cov[,2] 
  }
  
  # Create and save plots for later viewing
  data_plot <- data.frame(x = cov_plot,
                          cent = preds_cent,
                          lcl = preds_lcl,
                          ucl = preds_ucl)
  cov_name <- str_remove(covariate, "_z")
  
  cred_interval <- (upper_ci - lower_ci) * 100
  yaxis_label <- paste0("Predicted detection probability (", 
                        cred_interval, 
                        "% CI)")
  
  marg_plot <- ggplot(data = data_plot, aes(x = x)) + 
    geom_line(aes(y = cent), col = line_color) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = transparency) +
    labs(x = covariates$axis_label[covariates$short_name == cov_name],
         y = yaxis_label) +
    theme_classic()
  
  return(marg_plot)
}
