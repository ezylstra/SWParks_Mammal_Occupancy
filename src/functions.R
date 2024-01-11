################################################################################
# Functions to be used in any script in SWParks_Mammal_Occupancy repo

# ER Zylstra
# Updated 2023-12-06
################################################################################

# List of functions included:
  # create_cov_list()
  # vegclass_estimates()
  # mean_estimate()
  # marginal_plot_occ()
  # marginal_plot_det()
  # paNA() 
  # propNA() 
  # occ_time_plot()

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
    str_remove(pattern = " [+] [(]1 [|] site[)]") %>%
    str_remove(pattern = " [+] [(]1 [|] years[)]") %>%
    str_remove_all(pattern = "I[(]") %>%
    str_remove_all(pattern = "[)]") %>%
    str_remove_all(pattern = "\\^2") %>%
    str_split_1(pattern = " [+] ")

    return(z_list)
}

#------------------------------------------------------------------------------#
# parameter_estimates: Create a table with parameter estimates on logit scale 
# from occurrence or detection part of model
#------------------------------------------------------------------------------#

# INPUTS
# model: output from spOccupancy single-season model
# parameter: character indicating whether to extract estimates for parameters in
  # the occurrence or detection part of the model
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)

# RETURNS
# est_table: a dataframe with mean, median, SD, credible interval, Rhat values,
  # and ESS for each parameter in occurrence or detection part of a model

parameter_estimates <- function(model, 
                                parameter = c("occ", "det"),
                                lower_ci = 0.025,
                                upper_ci = 0.975) {
  
  parameter <- match.arg(arg = parameter)
  
  if (parameter == "occ") {
    samples <- model$beta.samples
    submodel <- "occupancy"
    index <- 1
  } else {
    samples <- model$alpha.samples
    submodel <- "detection"
    index <- 2
  }
  
  # Create table to hold results
  est_table <- data.frame(Parameter = colnames(samples),
                          Mean = NA,
                          SD = NA,
                          Lower = NA,
                          Median = NA,
                          Upper = NA,
                          Rhat = round(model$rhat[[index]], 3),
                          ESS = round(model$ESS[[index]]),
                          f = NA,
                          exclude0 = NA,
                          row.names = NULL)
  
  est_table$Mean <- round(apply(samples, 2, mean), 3)
  est_table$SD <- round(apply(samples, 2, sd), 3)
  est_table$Lower <- round(apply(samples, 2, quantile, lower_ci), 3)
  est_table$Median <- round(apply(samples, 2, median), 3)
  est_table$Upper <- round(apply(samples, 2, quantile, upper_ci), 3)
  est_table$f <- round(apply(samples, 2, 
                             function (x) if_else(mean(x) > 0, 
                                                  sum(x > 0) / length(x),
                                                  sum(x < 0) / length(x))), 3)
  est_table$exclude0 <- ifelse(est_table$Lower * est_table$Upper > 0, "yes", "no")
  ci_name <- (upper_ci - lower_ci) * 100
  names(est_table)[names(est_table) == "Lower"] <- paste0("Lower", ci_name, "%")
  names(est_table)[names(est_table) == "Upper"] <- paste0("Upper", ci_name, "%")
  
  return(est_table)
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
# vegclass_table: a dataframe with mean, SD, and 95% CI for occupancy/detection
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
# mean_estimate: Calculate mean occurrence or detection probability (or  
# probability when all other covariates set to 0 [mean, for standardized covs])
#------------------------------------------------------------------------------#

# INPUTS
# model: output from spOccupancy single-season model
# parameter: character indicating whether to calculate occurrence or detection 
  # probabilities 
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)

# RETURNS
# est_table: a dataframe with mean, SD, and 95% CI for occupancy/detection
  # probability

mean_estimate <- function(model, 
                          parameter = c("occ", "det"),
                          lower_ci = 0.025,
                          upper_ci = 0.975) {
  
  parameter <- match.arg(arg = parameter)
  
  if (parameter == "occ") {
    samples <- model$beta.samples
    submodel <- "occupancy"
  } else {
    samples <- model$alpha.samples
    submodel <- "detection"
  }
  
  # Create table to hold results
  est_table <- data.frame(parameter = submodel,
                          mean_prob = NA,
                          sd_prob = NA,
                          ci_lower = NA,
                          ci_upper = NA)
  
  # Probability of occupancy/detection 
  probs <- exp(samples[,"(Intercept)"])/(1 + exp(samples[,"(Intercept)"])) 
  est_table$mean_prob[1] <- mean(probs)
  est_table$sd_prob[1] <- sd(probs)
  est_table$ci_lower[1] <- quantile(probs, lower_ci)
  est_table$ci_upper[1] <- quantile(probs, upper_ci)
  
  return(est_table)
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
  
  cols <- str_subset(colnames(model$beta.samples), pattern = covariate)
  beta_samples <- model$beta.samples[,c("(Intercept)", cols)]
  X_cov <- seq(from = min(data_list$occ.covs[[covariate]]), 
               to = max(data_list$occ.covs[[covariate]]),
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
  # If it's a spatial covariate:
  if (any(str_detect(colnames(spatial_covs), covariate))) {
    if (str_detect(covariate, "_z")) {
      cov_mn <- mean(spatial_covs[,str_remove(covariate, "_z")])
      cov_sd <- sd(spatial_covs[,str_remove(covariate, "_z")])
      cov_plot <- X_cov[,2] * cov_sd + cov_mn
    } else {
      cov_plot <- X_cov[,2] 
    }
  } else {
  # Otherwise it's an annual covariate (other than years)
    if (str_detect(covariate, "_z")) {
      cov_mn <- mean(data_list$occ.covs[[str_remove(covariate, "_z")]])
      cov_sd <- sd(data_list$occ.covs[[str_remove(covariate, "_z")]])
      cov_plot <- X_cov[,2] * cov_sd + cov_mn
    } else {
      cov_plot <- X_cov[,2] 
    }
  }
  
  # Create and save plots for later viewing
  data_plot <- data.frame(x = cov_plot,
                          cent = preds_cent,
                          lcl = preds_lcl,
                          ucl = preds_ucl)
  cov_name <- ifelse(covariate == "burn_severity_2011_z",
                     "burn",
                     str_remove(covariate, "_z"))

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
  
  cols <- str_subset(colnames(model$alpha.samples), pattern = covariate)
  alpha_samples <- model$alpha.samples[,c("(Intercept)", cols)]
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
  # If it's a spatial covariate:
  if (any(str_detect(colnames(spatial_covs), covariate))) {
    if (str_detect(covariate, "_z")) {
      cov_mn <- mean(spatial_covs[,str_remove(covariate, "_z")], )
      cov_sd <- sd(spatial_covs[,str_remove(covariate, "_z")])
      cov_plot <- X_cov[,2] * cov_sd + cov_mn
    } else {
      cov_plot <- X_cov[,2] 
    }
  } else {
    # Otherwise it's an annual or survey covariate
    if (str_detect(covariate, "_z")) {
      cov_mn <- mean(data_list$det.covs[[str_remove(covariate, "_z")]], 
                     na.rm = TRUE)
      cov_sd <- sd(data_list$det.covs[[str_remove(covariate, "_z")]], 
                   na.rm = TRUE)
      cov_plot <- X_cov[,2] * cov_sd + cov_mn
    } else {
      cov_plot <- X_cov[,2] 
    }
  }  

  # Create and save plots for later viewing
  data_plot <- data.frame(x = cov_plot,
                          cent = preds_cent,
                          lcl = preds_lcl,
                          ucl = preds_ucl)
  cov_name <- ifelse(covariate == "burn_severity_2011_z",
                     "burn",
                     str_remove(covariate, "_z"))

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

#------------------------------------------------------------------------------#
# occ_time_plot: Create figure that depicts mean predicted occurrence
# probability over time (including REs if they're in the model), as well as 
# naive occurrence probabilities and the estimated trend (if "years" is in the 
# model)
#------------------------------------------------------------------------------#

# INPUTS
# model: output from spOccupancy single-season model
# data_list: list of data required to run model in spOccupancy 
# covariate_table: table with information and axis labels for covariates
# central_meas: metric for summarizing predictions (mean [default] or median)
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)
# line_color: color for line depicting mean/median predicted probability
# transparency: alpha value specifying transparency of shaded CI

# RETURNS
# occ_time_plot: ggplot object depicting estimates/raw values on probability scale

occ_time_plot <- function(model, 
                          data_list,
                          covariate_table,
                          raw_occ = TRUE,
                          central_meas = c(mean, median),
                          lower_ci = 0.025,
                          upper_ci = 0.975,
                          line_color = "forestgreen",
                          transparency = 0.2) {
  
  # Identify whether trend and/or yearly REs are in the model
  trend <- ifelse(any(str_detect(colnames(model$beta.samples), "years_z")), 1, 0)
  yrRE <- ifelse(dim(model$X.re)[3] == 2, 1, 0)
  
  cred_interval <- (upper_ci - lower_ci) * 100
  yaxis_label <- paste0("Proportion of area used (", 
                        cred_interval, 
                        "% CI)")
  
  # Calculate the number of sites with at least one detection in a year
  site_dets <- apply(data_list$y, c(1, 2), paNA)
  # Proportion of sites with a detection in each year
  raw_occ_prob <- apply(site_dets, 2, mean, na.rm = TRUE)
  raws <- data.frame(x = YEARS,
                     cent = raw_occ_prob,
                     lcl = NA,
                     ucl = NA,
                     type = "Naive",
                     row.names = NULL)
  
  if (trend == 1) {
    tr_col <- str_subset(colnames(model$beta.samples), pattern = "years_z")
    trend_samples <- model$beta.samples[,c("(Intercept)", tr_col)]
    X_trend <- seq(from = min(data_list$occ.covs[["years_z"]]), 
                   to = max(data_list$occ.covs[["years_z"]]),
                   length = 100)
    X_trend <- cbind(1, X_trend)
    preds_tr <-  X_trend %*% t(trend_samples)
    preds_tr <- exp(preds_tr)/(1 + exp(preds_tr))
    preds_tr_cent <- apply(preds_tr, 1, central_meas)
    preds_tr_lcl <- apply(preds_tr, 1, quantile, lower_ci)
    preds_tr_ucl <- apply(preds_tr, 1, quantile, upper_ci)
    
    yr_mn <- mean(data_list$occ.covs[["years"]])
    yr_sd <- sd(data_list$occ.covs[["years"]])
    yr_plot <- X_trend[,2] * yr_sd + yr_mn
    
    trend_df <- data.frame(x = yr_plot,
                           cent = preds_tr_cent,
                           lcl = preds_tr_lcl,
                           ucl = preds_tr_ucl)
  } 
  
  # Get annual estimates (for any model) 
  # (we're assuming a maximum of one annual covariate in the model)
  ann_covs <- c("years_z", "traffic_z", "visits_z", "monsoon_ppt_z", "ppt10_z")
  ann_cols <- str_subset(colnames(model$beta.samples), 
                         pattern = paste0(ann_covs, collapse = "|"))
  ann_samples <- model$beta.samples[,c("(Intercept)", ann_cols)]
  if (!is.null(dim(ann_samples))) {
    ann_values <- data_list$occ.covs[ann_cols][[1]][1,]
    X_ann <- cbind(1, ann_values)
  } else {
    X_ann <- as.matrix(data.frame(int = rep(1, length(YEARS))))
  }
  preds_ann <- X_ann %*% t(ann_samples)
  if (yrRE == 1) {
    yrREcols <- grepl("years", colnames(model$beta.star.samples))
    yrREs <- t(model$beta.star.samples[,yrREcols])
    preds_ann <- preds_ann + yrREs
  }
  preds_ann <- exp(preds_ann)/(1 + exp(preds_ann))
  preds_ann_cent <- apply(preds_ann, 1, central_meas)
  preds_ann_lcl <- apply(preds_ann, 1, quantile, lower_ci)
  preds_ann_ucl <- apply(preds_ann, 1, quantile, upper_ci)
  
  ann_df <- data.frame(x = YEARS,
                       cent = preds_ann_cent,
                       lcl = preds_ann_lcl,
                       ucl = preds_ann_ucl,
                       type = "Estimate")
  ann_df <- rbind(ann_df, raws)
  
  if (trend == 1 & raw_occ) {
    occ_time_plot <- ggplot(data = trend_df, aes(x = x)) + 
      geom_line(aes(y = cent), col = line_color) +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = transparency) +
      geom_segment(data = filter(ann_df, type == "Estimate"), 
                   aes(x = x, xend = x, y = lcl, yend = ucl)) +
      geom_point(data = ann_df, shape = 21,
                 aes(x = x, y = cent, group = type, fill = type)) +
      scale_fill_manual(values = c("black", "white")) +
      labs(x = "Year", y = yaxis_label) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_classic() +
      theme(legend.position = c(0.02, 0.02),
            legend.justification = c("left", "bottom"),
            legend.title = element_blank())
  }
  if (trend == 1 & !raw_occ) {  
    occ_time_plot <- ggplot(data = trend_df, aes(x = x)) + 
      geom_line(aes(y = cent), col = line_color) +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = transparency) +
      geom_segment(data = filter(ann_df, type == "Estimate"), 
                   aes(x = x, xend = x, y = lcl, yend = ucl)) +
      geom_point(data = filter(ann_df, type == "Estimate"), 
                 aes(x = x, y = cent), shape = 21, fill = "black") +
      labs(x = "Year", y = yaxis_label) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_classic()
  }
  if (trend == 0 & raw_occ) {
    occ_time_plot <- ggplot(ann_df, aes(x = x)) +
      geom_segment(data = filter(ann_df, type == "Estimate"), 
                   aes(xend = x, y = lcl, yend = ucl)) +
      geom_point(data = ann_df, shape = 21,
                 aes(y = cent, group = type, fill = type)) +
      scale_fill_manual(values = c("black", "white")) +
      labs(x = "Year", y = yaxis_label) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_classic() +
      theme(legend.position = c(0.02, 0.02),
            legend.justification = c("left", "bottom"),
            legend.title = element_blank())
  }  
  
  if (trend == 0 & !raw_occ) {
    occ_time_plot <- ggplot(ann_df, aes(x = x)) +
      geom_segment(data = filter(ann_df, type == "Estimate"), 
                   aes(xend = x, y = lcl, yend = ucl)) +
      geom_point(data = filter(ann_df, type == "Estimate"), 
                 aes(y = cent), shape = 21, fill = "black") +
      labs(x = "Year", y = yaxis_label) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_classic()
  }  
  return(occ_time_plot)
}

#------------------------------------------------------------------------------#
# det_cat_estimates: Create a table with detection probabilities for different
# combinations of one or more categorical covariates
#------------------------------------------------------------------------------#

# INPUTS
# model: output from spOccupancy single-season model (note: a categorical 
# covariate, like camera or lens_2023, must have been included as a covariate in 
# the model formula for detection)
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)

# RETURNS
# det_table: a dataframe with mean, SD, and 95% CI for detection probabilities 
# associated with each level of the covariate

det_cat_estimates <- function(model, 
                              lower_ci = 0.025,
                              upper_ci = 0.975) {
  
  samples <- model$alpha.samples
  covs <- colnames(samples)
  
  # Identify covariate(s)
  cols <- covs[covs %in% c("camera", "lens_2023")]
  
  if (length(cols) == 0) {
    stop("A categorical covariate with 2 levels must be included in model for detection.")
  }
  
  name1 <- cols[1]
  name2 <- cols[2]
  
  # Create table to hold results
  if (length(cols) == 1) {
    det_table <- data.frame(det1 = 0:1) %>%
      rename_with(~name1, det1)
  } else {
    det_table <- data.frame(det1 = c(0, 0, 1, 1),
                            det2 = c(0, 1, 0, 1)) %>%
      rename_with(~c(name1, name2), c(det1, det2))
  }
  det_table <- det_table %>%
    mutate(mean_prob = NA,
           sd_prob = NA,
           ci_lower = NA,
           ci_upper = NA)
  
  det_classes <- det_table %>%
    select(-c(mean_prob, sd_prob, ci_lower, ci_upper)) %>%
    as.matrix()
  det_classes <- cbind(1, det_classes)
  predsl <- samples[, c("(Intercept)", cols)] %*% t(det_classes)
  preds <- exp(predsl) / (1 + exp(predsl))
  det_table$mean_prob <- apply(preds, 2, mean)
  det_table$sd_prob <- apply(preds, 2, sd)
  det_table$ci_lower <- apply(preds, 2, quantile, probs = lower_ci)
  det_table$ci_upper <- apply(preds, 2, quantile, probs = upper_ci)
  
  return(det_table)
}


#------------------------------------------------------------------------------#
# DEPRECATED
# trend_occ_plot: Create figure depicting trend in occurrence probability over 
# time (only for multi-season models with implicit dynamics)
#------------------------------------------------------------------------------#

# INPUTS
# model: output from spOccupancy single-season model
# data_list: list of data required to run model in spOccupancy 
# covariate_table: table with information and axis labels for covariates
# central_meas: metric for summarizing predictions (mean [default] or median)
# lower_ci: quantile for lower bound of credible interval (0.025 for 95% CI)
# upper_ci: quantile for upper bound of credible interval (0.975 for 95% CI)
# line_color: color for line depicting mean/median predicted probability
# transparency: alpha value specifying transparency of shaded CI

# RETURNS
# trend_plot: ggplot object depicting the trend on probability scale

trend_plot_occ <- function(model, 
                           data_list,
                           covariate_table,
                           raw_occ = TRUE,
                           central_meas = c(mean, median),
                           lower_ci = 0.025,
                           upper_ci = 0.975,
                           line_color = "forestgreen",
                           transparency = 0.2) {
  
  cols <- str_subset(colnames(model$beta.samples), pattern = "years_z")
  beta_samples <- model$beta.samples[,c("(Intercept)", cols)]
  X_cov <- seq(from = min(data_list$occ.covs[["years_z"]]), 
               to = max(data_list$occ.covs[["years_z"]]),
               length = 100)
  X_cov <- cbind(1, X_cov)
  
  preds <-  X_cov %*% t(beta_samples)
  preds <- exp(preds)/(1 + exp(preds))
  preds_cent <- apply(preds, 1, central_meas)
  preds_lcl <- apply(preds, 1, quantile, lower_ci)
  preds_ucl <- apply(preds, 1, quantile, upper_ci)
  
  # Identify covariate values for x-axis (on original scale)
  cov_mn <- mean(data_list$occ.covs[["years"]])
  cov_sd <- sd(data_list$occ.covs[["years"]])
  cov_plot <- X_cov[,2] * cov_sd + cov_mn
  
  # Create and save plots for later viewing
  data_plot <- data.frame(x = cov_plot,
                          cent = preds_cent,
                          lcl = preds_lcl,
                          ucl = preds_ucl)
  
  cred_interval <- (upper_ci - lower_ci) * 100
  yaxis_label <- paste0("Proportion of area used (", 
                        cred_interval, 
                        "% CI)")
  
  if (raw_occ) {
    
    # Calculate the number of sites with at least one detection in a year
    site_dets <- apply(data_list$y, c(1, 2), paNA)
    # Proportion of sites with a detection in each year
    raw_occ_prob <- apply(site_dets, 2, mean, na.rm = TRUE)
    
    yr_min <- min(data_list$occ.covs$years, na.rm = TRUE)
    yr_max <- max(data_list$occ.covs$years, na.rm = TRUE)
    raws <- data.frame(yr = yr_min:yr_max,
                       raw_occ = raw_occ_prob,
                       row.names = NULL)
    
    trend_plot <- 
      ggplot(data = data_plot, aes(x = x)) + 
      geom_line(aes(y = cent), col = line_color) +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = transparency) +
      geom_point(data = raws[!is.na(raws$raw_occ), ], 
                 aes(x = yr, y = raw_occ), col = "black") +
      labs(x = "Year", y = yaxis_label) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_classic()
    
  } else {
    
    trend_plot <- ggplot(data = data_plot, aes(x = x)) + 
      geom_line(aes(y = cent), col = line_color) +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = transparency) +
      labs(x = "Year", y = yaxis_label) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_classic()
  }
  
  return(trend_plot)
}
