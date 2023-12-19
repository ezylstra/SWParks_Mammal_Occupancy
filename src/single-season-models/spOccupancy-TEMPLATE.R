################################################################################
# Template to run and evaluate a suite of single-season occupancy models for 
# a given park, year, and species (using the spOccupancy package)

# Objects that need to be specified to run models (ie, to create 
# src/single-season-model/YEAR/spOccupancy-PARK-YEAR-SPECIES.R)
  # PARK, YEAR, SPECIES (lines 50, 53, 67)
  # Specifications for occurrence models: OCC_NULL, OCC_MODELS
  # Specifications for detection models: DET_NULL, DET_MODELS
  # MCMC parameters: N_SAMPLES, N_BURN, N_THIN, N_CHAINS
  # Method used to select a "best" model for inferences: STAT

# ER Zylstra
# Updated 2023-10-13
################################################################################

#------------------------------------------------------------------------------#
# Load packages and custom functions
#------------------------------------------------------------------------------#

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

# Select park of interest ("CHIR", "ORPI", or "SAGW")
PARK <- "SAGW"

source("src/photo-data/format-mammal-data.R")

# dat = information about each photo (date, time, species, location)
# events = information about each camera deployment (dates, location, duration)
# event_mat = camera location x day matrix with 1/0 indicating whether camera
#             was operational or not
# locs = information about each camera location (park, lat/long, name)
# species = table with species observed (species code, common name, # of obs)

# Load functions
source("src/functions.R")

#------------------------------------------------------------------------------#
# Specify model parameters
#------------------------------------------------------------------------------#

# Select year of interest
YEAR <- 2023

# Look at detection data for various species
detects <- read.csv(paste0("output/species-detections-byyr-", PARK, ".csv"))
detects <- detects %>%
  dplyr::filter(yr == YEAR) %>%
  arrange(desc(propdetect))
# View just those species with a detection rate of 5% (propdetect = proportion 
# of nobs [camera locations * sampling occasion] with species detection)
detects %>% 
  left_join(species,by=c("spp" = "Species_code")) %>%
  select(c(spp, Species, Common_name, nobs, propdetect))

# Select species of interest (ideally with a detection rate of at least 5%)
SPECIES <- "URAM"

# Save this script as: 
# src/single-season-models/YEAR/PARK/spOccupancy-PARK-YEAR-SPECIES.R

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

# Load dataframe with information about covariates:
covariates <- read.csv("data/covariates/covariates.csv", header = TRUE)

#------------------------------------------------------------------------------#
# Information about how to build candidate model sets
#------------------------------------------------------------------------------#

# To specify which covariates to include in models for occurrence or detection,
# we'll use the short_name column in the covariates dataframe.

# To create candidate model sets for occurrence, we'll define 1 or 2 objects:

  # OCC_NULL: a logical indicating whether a null model should be included the 
  # candidate model set
  
  # OCC_MODELS: a list, where each element is a vector listing one of more 
  # covariates to include in a single candidate model. For instance, 
  # OCC_MODELS <- list(c("years"), c("years", "roads")) indicates that we will 
  # include two models in our candidate set, one in which occurrence 
  # probabilities vary with year, and one in which occurrence probabilities vary
  # with both year and distance to roads (occ prob ~ years + roads).

  # We'll use the same process to create candidate model sets for detection:
  # DET_NULL, DET_MODELS

#------------------------------------------------------------------------------#
# Explore the detection process
#------------------------------------------------------------------------------#

# See what covariates are available for detection part of model
covariates %>%
  dplyr::filter(parameter %in% c("either", "detection")) %>%
  dplyr::filter(park %in% c(PARK, "all")) %>%
  select(-c(parameter, park))

# Logical indicating whether a null model for detection should be included in 
# the candidate model set
DET_NULL <- FALSE

# Create a full model for detection
DET_MODELS <- list(c("burn", "day2", "deploy_exp", "effort"))

# Run a model model with no covariates on occurrence, full model for detection
OCC_NULL <- TRUE

# No site random effects (for now)
SITE_RE_OCC <- "none"
SITE_RE_DET <- "none"

# Use OCC and DET objects to create formulas for candidate models:
source("src/single-season-models/spOccupancy-create-model-formulas.R")
message("Check candidate models:", sep = "\n")
model_specs

# Set MCMC parameters
N_SAMPLES <- 8000
N_BURN <- 4000
N_THIN <- 8
N_CHAINS <- 3

# Run model using spOccupancy package
source("src/single-season-models/spOccupancy-run-candidate-models.R")

# View summary table, ranked by WAIC
model_stats %>% arrange(waic)

# Look at results
best_index <- 1
summary(out_list[[best_index]])
samps <- cbind(out_list[[best_index]]$beta.samples[, -1],
               out_list[[best_index]]$alpha.samples[, -1])
(f <- apply(samps, 2, function(x) ifelse(mean(x) > 0, sum(x > 0) / length(x),
                                         sum(x < 0) / length(x))))

# Pick which detection covariates to use in subsequent candidate models
# (for now, can use covariates with f >= 0.9)
DET_MODELS <- list(c("burn", "effort")) 

#------------------------------------------------------------------------------#
# Specify the occurrence portion of candidate models
#------------------------------------------------------------------------------#

# See what covariates are available for occurrence part of model
covariates %>%
  dplyr::filter(parameter %in% c("either", "occupancy")) %>%
  dplyr::filter(park %in% c(PARK, "all")) %>%
  select(-c(parameter, park))

# View those pairs of continuous covariates that are highly correlated
cor_df %>%
  arrange(desc(corr)) %>%
  dplyr::filter(abs(corr) >= 0.7)

# Logical indicating whether a null model for occurrence should be included in 
# the candidate model set
OCC_NULL <- FALSE

# There are 4 categories of spatial covariates (though each park only has 
# covariates in 2 or 3 of the categories):
  # topographic: aspect, elev, slope
    # (using linear rather than quadratic forms of elev & slope because SAGW 
    # doesn't span that large of a range and we often get nonsensical results 
    # with highest probabilities at extreme values)
  # veg: vegclasses + wash (for now, only available for SAGW)
  # burn: burn severity classes for 2011 fire (only available in CHIR)
  # anthropogenic: roads, boundary, trails, pois, roadbound, trailpois

# For occurrence part of the models, try including item(s) from each category of
# spatial covariates, excluding any covariates that are highly correlated 
# (|r| >= 0.7).

# Pick covariates to include candidate models
OCC_MODELS <- list(c("aspect", "veg", "wash", "burn", "roads"),
                   c("elev", "veg", "wash", "burn", "roads"),
                   c("slope", "veg", "wash", "burn", "roads"),
                   c("aspect", "veg", "wash", "burn", "boundary"),
                   # c("elev", "veg", "wash", "burn", "boundary"),
                   c("slope", "veg", "wash", "burn", "boundary"),
                   c("aspect", "veg", "wash", "burn", "trail"),
                   c("elev", "veg", "wash", "burn", "trail"),
                   c("slope", "veg", "wash", "burn", "trail"),
                   c("aspect", "veg", "wash", "burn", "pois"),
                   c("elev", "veg", "wash", "burn", "pois"),
                   c("slope", "veg", "wash", "burn", "pois"),
                   c("aspect", "veg", "wash", "burn", "roadbound"),
                   c("elev", "veg", "wash", "burn", "roadbound"),
                   c("slope", "veg", "wash", "burn", "roadbound"),
                   c("aspect", "veg", "wash", "burn", "trailpoi"),
                   c("elev", "veg", "wash", "burn", "trailpoi"),
                   c("slope", "veg", "wash", "burn", "trailpoi"))

# Use OCC and DET objects to create formulas for candidate models:
source("src/single-season-models/spOccupancy-create-model-formulas.R")
message("Check candidate models:", sep = "\n")
model_specs

# Run candidate models
source("src/single-season-models/spOccupancy-run-candidate-models.R")
# Note: this will often take several minutes to run

# View summary table, ranked by WAIC
model_stats %>% arrange(waic)

#------------------------------------------------------------------------------#
# Select "best" model
#------------------------------------------------------------------------------#

# Identify a model to use for inferences.  Can base this on WAIC or deviance 
# from k-fold CV.  Alternatively, can select another model by setting STAT to
# "model_no" and specifying the "best_index" directly.

# Specify STAT as either: waic or model_no
STAT <- "model_no"   
if (STAT == "model_no") {
  # If STAT == "model_no", specify model of interest by model number in table
  best_index <- 6
} else {
  min_stat <- min(model_stats[,STAT])
  best_index <- model_stats$model_no[model_stats[,STAT] == min_stat] 
}

# Look at covariate effects (for both occurrence and detection). In general, 
# covariates with f values close to or above 0.9 have more explanatory power.
summary(out_list[[best_index]])
samps <- cbind(out_list[[best_index]]$beta.samples[, -1],
               out_list[[best_index]]$alpha.samples[, -1])
(f <- apply(samps, 2, function(x) ifelse(mean(x) > 0, sum(x > 0) / length(x),
                                         sum(x < 0) / length(x))))

# If all the covariates have sufficient explanatory power, then move on to the 
# next step (looking at results from best model and saving model object, below). 
# If not, run a model that excludes covariates with little to no explanatory 
# power from the model for occurrence.

  # Change occupancy part of model (if needed)
  # OCC_NULL <- FALSE
  # OCC_MODELS <- list(c("burn", "elev"))

  # Change detection part of model (if needed)
  # DET_NULL <- TRUE
  # DET_MODELS <- list(c("day", "effort"))
  # rm(DET_MODELS)
  
  # source("src/single-season-models/spOccupancy-create-model-formulas.R")
  # message("Check candidate models:", sep = "\n")
  # model_specs
  # 
  # # Run model(s)
  # source("src/single-season-models/spOccupancy-run-candidate-models.R")
  # model_stats %>% arrange(waic)
  # 
  # # Specify STAT as either: waic or model_no
  # STAT <- "waic"
  # if (STAT == "model_no") {
  #   # If STAT == "model_no", specify model of interest by model number in table
  #   best_index <- 4
  # } else {
  #   min_stat <- min(model_stats[,STAT])
  #   best_index <- model_stats$model_no[model_stats[,STAT] == min_stat]
  # }
  # # Look at model output and f values
  # summary(out_list[[best_index]])
  # samps <- cbind(out_list[[best_index]]$beta.samples[, -1],
  #                out_list[[best_index]]$alpha.samples[, -1])
  # (f <- apply(samps, 2, function(x) ifelse(mean(x) > 0, sum(x > 0) / length(x),
  #                                          sum(x < 0) / length(x))))

# View trace plots
# plot(best$beta.samples, density = FALSE)
# plot(best$alpha.samples, density = FALSE)

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

# If you received a warning that PPC indicates you have not adequately described
# spatial variation in occupancy and/or detection (ppc.site < 0.1 or > 0.9), 
# try adding random site effects to the occupancy and/or detection portion of 
# the model (code below is not elegant, but it works)

  # Specify which covariates you want (probably the same as what was in best model 
  # above, but specify to be sure)
  # OCC_MODELS <- list(c("elev", "roads", "burn"))
  # DET_MODELS <- list(c("burn", "day"))
  
  # SITE_RE_DET <- "none"
  # SITE_RE_OCC <- "none"
  # source("src/single-season-models/spOccupancy-create-model-formulas.R")
  # model_specs1 <- model_specs
  # occ_formulas1 <- occ_formulas
  # det_formulas1 <- det_formulas
  # 
  # SITE_RE_DET <- "unstructured"
  # source("src/single-season-models/spOccupancy-create-model-formulas.R")
  # model_specs2 <- model_specs
  # occ_formulas2 <- occ_formulas
  # det_formulas2 <- det_formulas
  # 
  # SITE_RE_DET <- "none"
  # SITE_RE_OCC <- "unstructured"
  # source("src/single-season-models/spOccupancy-create-model-formulas.R")
  # model_specs3 <- model_specs
  # occ_formulas3 <- occ_formulas
  # det_formulas3 <- det_formulas
  # 
  # SITE_RE_DET <- "unstructured"
  # source("src/single-season-models/spOccupancy-create-model-formulas.R")
  # model_specs4 <- model_specs
  # occ_formulas4 <- occ_formulas
  # det_formulas4 <- det_formulas
  # 
  # model_specs <- rbind(model_specs1, model_specs2, model_specs3, model_specs4)
  # occ_formulas <- c(occ_formulas1, occ_formulas2, occ_formulas3, occ_formulas4)
  # det_formulas <- c(det_formulas1, det_formulas2, det_formulas3, det_formulas4)
  # 
  # # Run models with all combinations of site random effects
  # source("src/single-season-models/spOccupancy-run-candidate-models.R")
  # model_stats %>% arrange(waic)
  
  # # Select the best model, but BEWARE of models that didn't converge (high rhat
  # # and/or low ESS)
  # STAT <- "model_no"
  # if (STAT == "model_no") {
  #   # If STAT == "model_no", specify model of interest by model number in table
  #   best_index <- 3
  # } else {
  #   min_stat <- min(model_stats[,STAT])
  #   best_index <- model_stats$model_no[model_stats[,STAT] == min_stat]
  # }
  # # Look at model output
  # summary(out_list[[best_index]])

# Identify model that will be used for inferences, extract output and formulas
best <- out_list[[best_index]]
best_psi_model <- model_specs[best_index, 1]
best_p_model <- model_specs[best_index, 2]
occ_estimates <- parameter_estimates(model = best, 
                                     parameter = "occ",
                                     lower_ci = 0.025,
                                     upper_ci = 0.975)
det_estimates <- parameter_estimates(model = best, 
                                     parameter = "det",
                                     lower_ci = 0.025,
                                     upper_ci = 0.975)
occ_estimates <- occ_estimates %>%
  rename(Covariate = Parameter) %>%
  mutate(Parameter = "Occurrence", .before = "Covariate")
det_estimates <- det_estimates %>%
  rename(Covariate = Parameter) %>%
  mutate(Parameter = "Detection", .before = "Covariate")
estimates <- rbind(occ_estimates, det_estimates)
estimates

# Extract names of covariates (with and without "_z" subscripts) from best model
psi_covs_z <- create_cov_list(best_psi_model)
if (length(psi_covs_z) == 1 & any(psi_covs_z == "1")) {
  psi_covs_z <- character(0)
  psi_covs <- character(0)
} else {
  psi_covs <- psi_covs_z %>% str_remove_all(pattern = "_z")
}
p_covs_z <- create_cov_list(best_p_model)
if (length(p_covs_z) == 1 & any(p_covs_z == "1")) {
  p_covs_z <- character(0)
  p_covs <- character(0)
} else {
  p_covs <- p_covs_z %>% str_remove_all(pattern = "_z")
}

# Save model object to file
model_filename <- paste0("output/single-season-models/", PARK, "-", YEAR, 
                         "-", SPECIES, ".rds")
model_list <- list(model = best, 
                   psi_model = best_psi_model,
                   p_model = best_p_model,
                   data = data_list)
saveRDS(model_list, file = model_filename)

#------------------------------------------------------------------------------#
# Calculate and visualize predicted probabilities of occupancy, across park
#------------------------------------------------------------------------------#

# Create spatial predictions IF there are spatial covariates in the occurrence
# part of the model. 
# NOTE: this can take several minutes to run.

if (length(psi_covs) > 0) {
  source("src/single-season-models/spOccupancy-predictions.R")
  # Note: this can take several minutes to run.
  
  # This script creates:
  # best_pred: a list with predictions for each raster cell, MCMC sample
  # preds_mn: a raster with mean values in each cell (across MCMC samples)
  # preds_sd: a raster with SDs in each cell (across MCMC samples)
  # plot_preds_mn: a ggplot object with predicted mean values across park
  # plot_preds_sd: a ggplot object with predicted sd values across park
  
  # Plot predicted means
  print(plot_preds_mn) 
  
  # Plot predicted sds
  print(plot_preds_sd)
  
  # Can save either of the plots to file (example below):
  # ggsave(filename = "C:/.../SPECIES_MeanOccupancy.jpg",
  #        plot = plot_preds_mn,
  #        device = "jpeg",
  #        width = 4,
  #        height = 4,
  #        units = "in",
  #        dpi = 600)  
}

#------------------------------------------------------------------------------#
# Calculate and create figures depicting marginal effects of covariates on 
# occurrence probability (predicted covariate effects assuming all other 
# covariates held constant)
#------------------------------------------------------------------------------#

# Identify continuous covariates in occurrence part of the best model
psi_continuous <- psi_covs_z[!psi_covs_z %in% c("1", "vegclass2", "vegclass3")]
psi_cont_unique <- unique(psi_continuous)
psi_n_cont <- length(psi_cont_unique)

# If there are any continuous covariates, create a figure for each:
if (psi_n_cont > 0) {
  # Loop through each covariate
  for (cov in psi_cont_unique) {
    # Create name of plot:
    plotname <- paste0("marginal_psi_", str_remove(cov, "_z"))
    # Create plot
    assign(plotname, 
           marginal_plot_occ(covariate = cov, 
                             model = best, 
                             data_list = data_list,
                             covariate_table = covariates,
                             central_meas = mean))
  } 
}

# Can view these plots, calling them by name. Available plots listed here:
str_subset(ls(), "marginal_psi_")

# Can print all to plot window in Rstudio:
for (fig in str_subset(ls(), "marginal_psi_")) {
  print(get(fig))
}
# Could also save any of the plots to file using ggsave()

# If vegetation classes were included as covariates in the model, extract
# occurrence probabilities for each class
if (sum(str_detect(psi_covs, "veg")) > 0) {
  occprobs_veg <- vegclass_estimates(model = best, 
                                     parameter = "occ")
  print(occprobs_veg)
}

# If there are no covariates in the model (ie, a null model), print overall 
# occurrence probability
if (psi_n_cont == 0 & length(psi_covs) == 0) {
  overall_occ <- mean_estimate(model = best, 
                               parameter = "occ",
                               lower_ci = 0.025,
                               upper_ci = 0.975)
  print(overall_occ)
}  

#------------------------------------------------------------------------------#
# Create figures depicting marginal effects of covariates on detection 
# probability (predicted covariate effects assuming all other covariates held
# constant)
#------------------------------------------------------------------------------#

# Identify continuous covariates in detection part of the best model
p_continuous <- p_covs_z[p_covs_z != "1"]
p_cont_unique <- unique(p_continuous)
p_n_cont <- length(p_cont_unique)

# If there are any continuous covariates, create a figure for each:
if (p_n_cont > 0) {
  # Loop through each covariate
  for (cov in p_cont_unique) {
    # Create name of plot:
    plotname <- paste0("marginal_p_", str_remove(cov, "_z"))
    # Create plot
    assign(plotname, 
           marginal_plot_det(covariate = cov, 
                             model = best, 
                             data_list = data_list,
                             covariate_table = covariates,
                             central_meas = mean))
  } 
}
# Can view these plots, calling them by name. Available plots listed here:
str_subset(ls(), "marginal_p_")

# Or print all to plot window:
for (fig in str_subset(ls(), "marginal_p_")) {
  print(get(fig))
}
# Could also save any of the plots to file using ggsave()

# If vegetation classes were included as covariates in the model, extract
# detection probabilities for each class
if (sum(str_detect(p_covs, "veg")) > 0) {
  detprobs_veg <- vegclass_estimates(model = best, 
                                     parameter = "det")
  print(detprobs_veg)
}

# If there are no covariates in the model (a null model), print overall 
# detection probability
if (p_n_cont == 0 & length(p_covs) == 0) {
  overall_det <- mean_estimate(model = best, 
                               parameter = "det",
                               lower_ci = 0.025,
                               upper_ci = 0.975)
  print(overall_det)
}  
