################################################################################
# Template to run and evaluate a suite of single-season occupancy models for 
# a given park, year, and species (using the spOccupancy package)

# Objects that need to be specified to run models (ie, to create 
# src/single-season-model/YEAR/spOccupancy-PARK-SPECIES-YEAR.R)
  # PARK, YEAR, SPECIES (lines 41, 44, 58)
  # Specifications for occurrence models: OCC_NULL, OCC_MODELS1, OCC_MODELS2
  # Specifications for detection models: DET_NULL, DET_MODELS1, DET_MODELS2
  # MCMC parameters: N_SAMPLES, N_BURN, N_THIN, N_CHAINS
  # Statistic used to select a "best" model for inferences: STAT

# ER Zylstra
# Updated 2023-02-02
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(spOccupancy)
library(ggplot2)
library(tidyterra)

#------------------------------------------------------------------------------#
# Load all photo, location, events, species data 
#------------------------------------------------------------------------------#

source("src/photo-data/format-mammal-data.R")

# dat = information about each photo (date, time, species, location)
# events = information about each camera deployment (dates, location, duration)
# event_mat = camera location x day matrix with 1/0 indicating whether camera
#             was deployed or not
# locs = information about each camera location (park, lat/long, name)
# species = table with species observed (species code, common name, # of obs)

#------------------------------------------------------------------------------#
# Specify model parameters
#------------------------------------------------------------------------------#

# Select park of interest ("CHIR", "ORPI", or "SAGW")
PARK <- "SAGW"

# Select year of interest
YEAR <- 2017

# Look at detection data for various species
detects <- read.csv("output/species-detections-byparkyr.csv", header = TRUE)
detects <- detects %>%
  dplyr::filter(Park == PARK & yr == YEAR) %>%
  arrange(desc(propdetect))
# View just those species with a detection rate of 5% 
# (n = camera location * sampling occasion)
detects %>% 
  dplyr::filter(propdetect >= 0.05) %>%
  select(c(spp, propdetect))

# Select species of interest (ideally with a detection rate of at least 5%)
SPECIES <- "LECA"

#------------------------------------------------------------------------------#
# Prepare detection and covariate data to run occupancy models with spOccupancy
#------------------------------------------------------------------------------#

source("src/single-season-models/spOccupancy-data-prep.R")

# This outputs data_list, which contains:
# y: a matrix with detection histories (nrow = no. sites; ncol = no. occasions)
# occ.covs: a dataframe with covariate values for each camera location (columns
  # whose names end with "_z" are standardized [z-scores])
# det.covs: a list of all covariates that could be used as covariates for
  # detection. Spatial covariates are vectors (length = no. sites). Survey
  # covariates are matrices (no. sites * no. occasions). Named objects that end
  # with "_z" are standardized)

#------------------------------------------------------------------------------#
# Identify covariates to be included in the occurrence and detection portions 
# of candidate models
#------------------------------------------------------------------------------#

# First, view those pairs of continuous covariates that are highly correlated
# (to help inform occupancy models below)
cor_df %>%
  arrange(desc(corr)) %>%
  dplyr::filter(abs(corr) > 0.5)

# Logicals indicating whether null models (for occurrence or detection) should
# be included in candidate model sets
OCC_NULL <- TRUE
DET_NULL <- TRUE

# Load dataframe with information about covariates:
covariates <- read.csv("data/covariates/covariates.csv", header = TRUE)

# Identify covariates to include in the occurrence part of candidate models

  # See what covariates are available
  covariates %>%
    dplyr::filter(parameter %in% c("either", "occupancy")) %>%
    dplyr::filter(park %in% c(PARK, "all")) %>%
    select(-c(parameter, park))
  
  # Pick covariates to include in simple candidate models via the short_name 
  # column in the covariates dataframe
  OCC_MODELS1 <- c("aspect", "elev2", "slope", "veg")
  
  # To combine covariates in a single candidate model, provide a vector of 
  # short_names e.g., c("aspect", "boundary") would create the following model
  # for occurrence: psi ~ east + north + boundary
  OCC_MODELS2 <- list(c("veg", "boundary"),
                      c("veg", "wash", "roads"))
  
# Identify covariates to include in the detection part of candidate models
    
  # See what covariates are available
  covariates %>%
    dplyr::filter(parameter %in% c("either", "detection")) %>%
    dplyr::filter(park %in% c(PARK, "all")) %>%
    select(-c(parameter, park))
  
  # Pick covariates to include in simple candidate models via the short_name 
  # column in the covariates dataframe
  DET_MODELS1 <- c("effort")
  # To combine different covariates in a candidate model, provide a vector of 
  # short_names
  DET_MODELS2 <- list(c("day2", "deploy", "effort"))

  # Use OCC and DET objects to create formulas for candidate models:
  source("src/single-season-models/spOccupancy-create-model-formulas.R")
    
  message("Check candidate models:", sep = "\n")
  model_specs

#------------------------------------------------------------------------------#
# Run models and compare fit
#------------------------------------------------------------------------------#

# Set MCMC parameters
N_SAMPLES <- 5000
N_BURN <- 3000
N_THIN <- 1
N_CHAINS <- 3

source("src/single-season-models/spOccupancy-run-candidate-models.R")
  # Note: this will often take several minutes to run

# View summary table, ranked by WAIC
model_stats %>%
  arrange(waic)

# Description of columns in summary table:
  # psi: formula for occurrence part of model
  # det: formula for detection part of model
  # max.rhat: maximum value of R-hat across model parameters (want value < 1.05)
  # min.ESS: minimum value of ESS (effective sample size) across model
    # parameters (want value > 400)
  # ppc.sites: posterior predictive checks when binning the data across sites. 
    # P-values < 0.1 or > 0.9 can indicate that model fails to adequately
    # represent variation in occurrence or detection across space.
  # ppc.reps: posterior predictive checks when binning the data across
    # replicates. P-values < 0.1 or > 0.9 can indicate that model fails to 
    # adequately represent variation in detection over time.
  # waic: WAIC (Widely Applicable Information Criterion) for comparing models
    # (lower is better)
  # k.fold.dev: Deviance from k-fold cross validation for comparing models 
    # (lower is better)

#------------------------------------------------------------------------------#
# Look at results and predictions from "best" model
#------------------------------------------------------------------------------#

# Identify a model to use for inferences.  Can base this on WAIC or deviance 
# from k-fold CV.  Alternatively, can select another model by specifying 
# the "best_index" directly.

# Specify STAT as either: waic, k.fold.dev, or model_no
STAT <- "model_no"   

if (STAT == "model_no") {
  # If STAT == "model_no", specify model of interest by model number in table
  best_index <- 20  
} else {
  min_stat <- min(model_stats[,STAT])
  best_index <- model_stats$model_no[model_stats[,stat] == min_stat] 
}

# Extract output from "best" model
best <- out_list[[best_index]]

# View covariate structure
best_psi_model <- model_specs[best_index, 1]
best_p_model <- model_specs[best_index, 2]
message("psi ", best_psi_model)
message("p ", best_p_model)

# Extract list of occurrence and detection covariates in best model
psi_covs_z <- best_psi_model %>%
  str_remove(pattern = "~ ") %>% 
  str_remove_all(pattern = "I[(]") %>%
  str_remove_all(pattern = "[)]") %>%
  str_remove_all(pattern = "\\^2") %>%
  str_split_1(pattern = " [+] ")
psi_covs <- psi_covs_z %>%
  str_remove_all(pattern = "_z")
p_covs_z <- best_p_model %>%
  str_remove(pattern = "~ ") %>% 
  str_remove_all(pattern = "I[(]") %>%
  str_remove_all(pattern = "[)]") %>%
  str_remove_all(pattern = "\\^2") %>%
  str_split_1(pattern = " [+] ")
p_covs <- p_covs_z %>%
  str_remove_all(pattern = "_z")

# Parameter estimates
summary(best)
# Note: this is good for viewing, but will want to use other means to create
# a table for reports/publications

# Trace plots
plot(best$beta.samples, density = FALSE)
plot(best$alpha.samples, density = FALSE)

# Posterior predictive checks (want Bayesian p-values between 0.1 and 0.9)
ppc.site <- as.numeric(model_stats$ppc.sites[model_stats$model_no == best_index])
ppc.rep <- as.numeric(model_stats$ppc.reps[model_stats$model_no == best_index])
if (ppc.site < 0.1 | ppc.site > 0.9) {
  warning(paste0("PPC indicates that we have not adequately described spatial ",
                 "variation in occupancy and/or detection."))
} else {
  cat(paste0("PPC indicates that we have adequately described spatial ",
             "variation in occupancy and detection."))
} 
if (ppc.rep < 0.1 | ppc.site > 0.9) {
  warning(paste0("PPC indicates that we have not adequately described temporal ",
                 "variation in detection."))
} else {
  cat(paste0("PPC indicates that we have adequately described temporal ",
                 "variation in detection."))
}

# If there's evidence that spatial variation isn't well explained, plot the 
# difference in the discrepancy measure between the replicate and actual data 
# across each of the sites (identify sites that are causing a lack of fit).
if (ppc.site < 0.1 | ppc.site > 0.9) {
  best_ppcs <- ppc.sites[[best_index]]
  diff_fit <- best_ppcs$fit.y.rep.group.quants[3, ] - best_ppcs$fit.y.group.quants[3, ]
  
  # Plot differences
  par(mfrow = c(1,1))
  plot(diff_fit, pch = 19, xlab = 'Site ID', ylab = "Replicate - True Discrepancy") 
  
  # Identify sites on a map
  prob_sites <- which(abs(diff_fit) > 0.4)
  plot(lat~long, data = spatial_covs, las = 1) # all camera locs
  points(lat~long, data = spatial_covs[prob_sites,], pch = 19, col = "blue")
}

#------------------------------------------------------------------------------#
# Calculate predicted probability of occupancy, across park
#------------------------------------------------------------------------------#

source("src/single-season-models/spOccupancy-predictions.R")
  # Note: this can take several minutes to run.

# This script creates:
  # best_pred: a list with predictions for each raster cell, MCMC sample
  # preds_mn: a raster with mean values in each cell (across MCMC samples)
  # preds_sd: a raster with SDs in each cell (across MCMC samples)
  # plot_preds_mn: a ggplot object with predicted mean values across park
  # plot_preds_sd: a ggplot object with predcited sd values across park

# Plot predicted means
plot_preds_mn 

# Plot predicted sds
plot_preds_sd

#------------------------------------------------------------------------------#
# Calculate and create figures depicting marginal effects of covariates on 
# occurrence probability (predicted covariate effects assuming all other 
# covariates held constant)
#------------------------------------------------------------------------------#

# Identify continuous covariates in occurrence part of the best model
psi_continuous <- psi_covs_z[!psi_covs_z %in% c("1", "vegclass2", "vegclass3")]
psi_cont_unique <- unique(psi_continuous)
psi_n_cont <- length(psi_cont_unique)

if (psi_n_cont > 0) {
  # Loop through each covariate
  for (cov in psi_cont_unique) {
    cols <- str_subset(colnames(best$beta.samples), pattern = cov)
    beta_samples <- best$beta.samples[,c("(Intercept)", cols)]
    X_cov <- seq(from = min(data_list$occ.covs[, cov]), 
                 to = max(data_list$occ.covs[, cov]),
                 length = 100)
    X_cov <- cbind(1, X_cov)
        
    # If there are quadratic effects, add column in X_cov
    if (ncol(beta_samples) == 3) {
      X_cov <- cbind(X_cov, X_cov[,2]^2)
    } 
     
    preds <-  X_cov %*% t(beta_samples)
    preds <- exp(preds)/(1 + exp(preds))
    preds_mn <- apply(preds, 1, mean)
    preds_lcl <- apply(preds, 1, quantile, 0.025)
    preds_ucl <- apply(preds, 1, quantile, 0.975)

    # Identify covariate values for x-axis (on original scale)
    if (str_detect(cov, "_z")) {
      cov_mn <- mean(data_list$occ.covs[,str_remove(cov, "_z")])
      cov_sd <- sd(data_list$occ.covs[,str_remove(cov, "_z")])
      cov_plot <- X_cov[,2] * cov_sd + cov_mn
    } else {
      cov_plot <- X_cov[,2] 
    }
    
    # Create and save plots for later viewing
    data_plot <- data.frame(x = cov_plot,
                            mn = preds_mn,
                            lcl = preds_lcl,
                            ucl = preds_ucl)
    cov_name <- str_remove(cov, "_z")
    assign(paste0("marginal_plot_psi_", cov_name),
           ggplot(data = data_plot, aes(x = x)) + 
             geom_line(aes(y = mn), col = "forestgreen") +
             geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
             labs(x = covariates$axis_label[covariates$short_name == cov_name], 
                  y = "Predicted occupancy probability (95% CI)") +
             theme_classic())
  }
}

# Can now view these plots, calling them by name. Available plots listed here:
str_subset(ls(), "marginal_plot_psi_")
# eg, if elevation was in the occurrence part of the model: 
# marginal_plot_psi_elev

# Extract mean occupancy probability for different vegetation classes, if 
# included as covariates in the best model
if (sum(str_detect(psi_covs, "veg")) > 0) {
  # Create table to hold results
  occprobs_veg <- data.frame(vegclass = 1:3,
                             mean_prob = NA,
                             sd_prob = NA,
                             ci_min = NA,
                             ci_max = NA)
  betas <- best$beta.samples  
  # Probability of occupancy in vegclass1 (reference level)
  vegclass1 <- exp(betas[,"(Intercept)"])/(1 + exp(betas[,"(Intercept)"])) 
  occprobs_veg$mean_prob[1] <- mean(vegclass1)
  occprobs_veg$sd_prob[1] <- sd(vegclass1)
  occprobs_veg$ci_min[1] <- quantile(vegclass1, 0.025)
  occprobs_veg$ci_max[1] <- quantile(vegclass1, 0.975)
  # Probability of occupancy in vegclass2
  vegclass2 <- betas[,"(Intercept)"] + betas[,"vegclass2"]
  vegclass2 <- exp(vegclass2)/(1 + exp(vegclass2)) 
  occprobs_veg$mean_prob[2] <- mean(vegclass2)
  occprobs_veg$sd_prob[2] <- sd(vegclass2)
  occprobs_veg$ci_min[2] <- quantile(vegclass2, 0.025)
  occprobs_veg$ci_max[2] <- quantile(vegclass2, 0.975)
  # Probability of occupancy in vegclass3
  vegclass3 <- betas[,"(Intercept)"] + betas[,"vegclass3"]
  vegclass3 <- exp(vegclass3)/(1 + exp(vegclass3)) 
  occprobs_veg$mean_prob[3] <- mean(vegclass3)
  occprobs_veg$sd_prob[3] <- sd(vegclass3)
  occprobs_veg$ci_min[3] <- quantile(vegclass3, 0.025)
  occprobs_veg$ci_max[3] <- quantile(vegclass3, 0.975)
  occprobs_veg
}

#------------------------------------------------------------------------------#
# Create figures depicting marginal effects of covariates on detection 
# probability (predicted covariate effects assuming all other covariates held
# constant)
#------------------------------------------------------------------------------#

p_covs_z <- best_p_model %>%
  str_remove(pattern = "~ ") %>% 
  str_remove_all(pattern = "I[(]") %>%
  str_remove_all(pattern = "[)]") %>%
  str_remove_all(pattern = "\\^2") %>%
  str_split_1(pattern = " [+] ")
p_covs <- p_covs_z %>%
  str_remove_all(pattern = "_z")

# Identify continuous covariates in occurrence part of the best model
p_continuous <- p_covs_z[p_covs_z != "1"]
p_cont_unique <- unique(p_continuous)
p_n_cont <- length(p_cont_unique)

if (p_n_cont > 0) {
  # Loop through each covariate
  for (cov in p_cont_unique) {
    cols <- str_subset(colnames(best$alpha.samples), pattern = cov)
    alpha_samples <- best$alpha.samples[,c("(Intercept)", cols)]
    X_cov <- seq(from = min(data_list$det.covs[[cov]], na.rm = TRUE), 
                 to = max(data_list$det.covs[[cov]], na.rm = TRUE),
                 length = 100)
    X_cov <- cbind(1, X_cov)
    
    # If there are quadratic effects, add column in X_cov
    if (ncol(alpha_samples) == 3) {
      X_cov <- cbind(X_cov, X_cov[,2]^2)
    } 
    
    preds <-  X_cov %*% t(alpha_samples)
    preds <- exp(preds)/(1 + exp(preds))
    preds_mn <- apply(preds, 1, mean)
    preds_lcl <- apply(preds, 1, quantile, 0.025)
    preds_ucl <- apply(preds, 1, quantile, 0.975)
    
    # Identify covariate values for x-axis (on original scale)
    if (str_detect(cov, "_z")) {
      if (str_remove(cov, "_z") %in% colnames(data_list$occ.covs)) {
        cov_mn <- mean(data_list$occ.covs[,str_remove(cov, "_z")])
        cov_sd <- sd(data_list$occ.covs[,str_remove(cov, "_z")])
      } else {
        cov_mn <- mean(data_list$det.covs[[str_remove(cov, "_z")]], na.rm = TRUE)
        cov_sd <- sd(data_list$det.covs[[str_remove(cov, "_z")]], na.rm = TRUE)
      }
      cov_plot <- X_cov[,2] * cov_sd + cov_mn
    } else {
      cov_plot <- X_cov[,2] 
    }
    
    # Create and save plots for later viewing
    data_plot <- data.frame(x = cov_plot,
                            mn = preds_mn,
                            lcl = preds_lcl,
                            ucl = preds_ucl)
    cov_name <- str_remove(cov, "_z")
    assign(paste0("marginal_plot_p_", cov_name),
           ggplot(data = data_plot, aes(x = x)) + 
             geom_line(aes(y = mn), col = "forestgreen") +
             geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
             labs(x = covariates$axis_label[covariates$short_name == cov_name], 
                  y = "Predicted detection probability (95% CI)") +
             theme_classic())
  }
}
