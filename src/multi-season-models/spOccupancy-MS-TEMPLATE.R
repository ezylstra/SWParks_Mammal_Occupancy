################################################################################
# Template to run and evaluate a suite of multi-season occupancy models for 
# a given park, set of years, and species (using the spOccupancy package)

# ER Zylstra
# Updated 2023-12-04
################################################################################

#------------------------------------------------------------------------------#
# Load packages and custom functions
#------------------------------------------------------------------------------#

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(abind)
library(terra)
library(spOccupancy)
library(ggplot2)
library(gridExtra)
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

# Select years of interest
YEARS <- 2017:2023

# Look at detection data for various species
detects <- read.csv(paste0("output/species-detections-", PARK, ".csv"))
detects <- arrange(detects, desc(propdetect))
# View species, noting those with a detection rate > 5% (propdetect = proportion 
# of nobs [camera locations * sampling occasion] with species detection)
detects %>% 
  left_join(species, by = c("spp" = "Species_code")) %>%
  select(c(spp, Species, Common_name, nobs, propdetect))

# Select species of interest (ideally with a detection rate of at least 5%)
SPECIES <- "ODHE"

# Save this script as: 
# src/multi-season-models/PARK/spOccupancy-PARK-FIRSTYEAR-LASTYEAR-SPECIES.R

#------------------------------------------------------------------------------#
# Prepare detection and covariate data to run occupancy models with spOccupancy
#------------------------------------------------------------------------------#

source("src/multi-season-models/spOccupancy-MS-data-prep.R")

# This outputs data_list, which contains:
# y: a array with detection histories (sites * years * occasions)
# occ.covs: a list of all covariates that could be used in the occurrence part 
  # of the model. Can vary among sites and/or years.
# det.covs: a list of all covariates that could be used in the detection part of
  # the model. Can vary among sites, years, and/or occasions
# coords: a matrix with UTM coordinates for each camera location (in WGS 84, 
  # Zone 12)

# Load dataframe with information about covariates:
covariates <- read.csv("data/covariates/covariates-MS.csv", header = TRUE)

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

# Random effects: 

  # It's worth considering whether we include unstructured random effects for 
  # both site and time in the occurrence part of the model. Site REs make sense 
  # in these implicit dynamic models because we want to account for 
  # non-independence of data from sites over multiple years. A yearly RE could 
  # also make sense because we don't have many covariates that can explain 
  # annual variation in occurrence (or temporal variation beyond that explained 
  # by a linear trend). However, we could run into problems with convergence if 
  # we include both spatial and temporal random effects given the relatively 
  # short time series. For now, we'll try including both and modify things, if
  # necessary, at a later step.

  # Specify whether we want spatial or temporal random effects in the occurrence
  # model by setting SITE_RE_OCC and TIME_RE_OCC as either "none" or 
  # "unstructured", respectively.  
  
  # Can use the same process to specify random effects in the detection model
  # using SITE_RE_DET AND TIME_RE_DET

#------------------------------------------------------------------------------#
# Specify and run first set of candidate models, where we will evaluate:
  # which, if any, annual covariates should be included in occurrence models
  # which, if any, covariates should be included in detection models
#------------------------------------------------------------------------------#

# For detection, use a "full" model
DET_NULL <- FALSE
DET_MODELS <- list(c("day2", "deploy_exp", "effort", "camera", "lens_2023"))

# For occurrence, try each annual covariate in a separate model
OCC_NULL <- TRUE
OCC_MODELS <- list("years", 
                   "visits", 
                   "traffic", 
                   "monsoon_ppt", 
                   "ppt10")

# Include a random site effect in occurrence model, but no other random effects
SITE_RE_OCC <- "unstructured"
TIME_RE_OCC <- "unstructured" 
SITE_RE_DET <- "none"
TIME_RE_DET <- "none"

# Create candidate model set
source("src/multi-season-models/spOccupancy-MS-create-model-formulas.R")
message("Check candidate models:", sep = "\n")
model_specs

# Set MCMC parameters
N_CHAINS <- 3
N_BURN <- 4000
N_THIN <- 15

# Instead of specifying the total number of samples (like we did for single-
# season models), we'll split the samples into a set of N_BATCH batches, each
# comprised of BATCH_LENGTH samples to improve mixing for the adaptive 
# algorithm. See documentation for the spOccupancy package for more info.
# Note that we're generating a fair number of samples to accommodate random 
# effects (had previously used burn = 2000, thin = 10, n_batch = 280)

N_BATCH <- 460
BATCH_LENGTH <- 25

n_samples <- N_BATCH * BATCH_LENGTH
message("MCMC specifications will result in a total of ",
        (n_samples - N_BURN) * N_CHAINS / N_THIN, 
        " posterior samples for each parameter.")

# Run first set of candidate models using spOccupancy package
source("src/multi-season-models/spOccupancy-MS-run-candidate-models.R")
# Note: this can take several minutes to run

# View summary table, ranked by WAIC
model_stats %>% arrange(waic)

# Select the annual covariate (years, visits, traffic, monsoon_ppt, or ppt10) 
# that should be included in the next set of candidate models. Typically, we'll
# select the covariate included in the model with the lowest WAIC. If none are 
# better than the null model, set BEST_ANNUAL <- NA, as random effects will 
# allow for variation in occurrence probability among years.
BEST_ANNUAL <- NA

# Look at parameter estimates for detection part of highest-ranking model and 
# decide what detection model we'd like to use in the next set of candidate 
# models. For now, suggest including any covariates with an f-value > 0.9
min_waic <- min(model_stats$waic)
summary(out_list[[model_stats$model_no[model_stats$waic == min_waic]]])
samps <- out_list[[model_stats$model_no[model_stats$waic == min_waic]]]$alpha.samples[, -1]
f_dets <- apply(samps, 2, 
                function(x) ifelse(mean(x) > 0, sum(x > 0) / length(x),
                                   sum(x < 0) / length(x)))
f_dets

# To use a null model for detection:
  # DET_NULL <- TRUE
  # rm(DET_MODELS)
# To use a model with a subset of those covariates, like day2 and effort:
  DET_MODELS <- list(c("day2", "deploy_exp", "effort", "lens_2023"))
# To use the same model, we can leave DET_MODELS as is.

#------------------------------------------------------------------------------#
# Specify and run second set of candidate models, where we will evaluate:
  # which, if any, spatial covariates should be included in occurrence models
#------------------------------------------------------------------------------#

# There are 3 categories of spatial covariates:
  # topographic: aspect, elev, slope
    # (using linear rather than quadratic forms of elev & slope because SAGW 
    # doesn't span that large of a range and we often get nonsensical results 
    # with highest probabilities at extreme values)
  # veg: vegclasses + wash
  # anthropogenic: roads, boundary, trails, pois, roadbound, trailpois

# For occurrence part of the models, try including item(s) from each category of
# spatial covariates + the annual covariate selected above, after removing
# any combinations that are highly correlated

cor_df %>%
  arrange(desc(corr)) %>%
  filter(abs(corr) >= 0.7)
  
scov_combos <- list(c("aspect", "veg", "wash", "burn", "roads"),
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
                    # c("elev", "veg", "wash", "burn", "roadbound"),
                    c("slope", "veg", "wash", "burn", "roadbound"),
                    c("aspect", "veg", "wash", "burn", "trailpoi"),
                    c("elev", "veg", "wash", "burn", "trailpoi"),
                    c("slope", "veg", "wash", "burn", "trailpoi"))
OCC_MODELS <- lapply(scov_combos, function(x) c(x, BEST_ANNUAL))
OCC_NULL <- FALSE

# Create candidate model set
source("src/multi-season-models/spOccupancy-MS-create-model-formulas.R")
message("Check candidate models:", sep = "\n")
model_specs

# Run second set of candidate models using spOccupancy package
source("src/multi-season-models/spOccupancy-MS-run-candidate-models.R")
# Note: this can take several minutes to run

# View summary table, ranked by WAIC
model_stats %>% arrange(waic)

# Description of columns in summary table:
  # psi: formula for occurrence part of model
  # det: formula for detection part of model
  # max.rhat: maximum value of R-hat across model fixed-effect parameters (want 
    # value < 1.05)
  # min.ESS: minimum value of ESS (effective sample size) across model
    # fixed-effect parameters  (want value > 400)
  # ppc.sites: posterior predictive checks when binning the data across sites
    # (within year). P-values < 0.1 or > 0.9 can indicate that model fails to 
    # adequately represent variation in occurrence or detection across space.
  # ppc.reps: posterior predictive checks when binning the data across 
    # replicates (within year). P-values < 0.1 or > 0.9 can indicate that model 
    # fails to adequately represent variation in detection across replicates.
  # waic: WAIC (Widely Applicable Information Criterion) for comparing models
    # (lower is better)

# Check that r-hat values and ESS look okay for most models.  If not, may 
# need to re-run after increasing N_BATCH or removing site REs. If problem
# persists, may need to remove some covariate combinations from consideration.

#------------------------------------------------------------------------------#
# Select "best" model
#------------------------------------------------------------------------------#

# Identify a model to use for inferences.  Can base this on WAIC or select 
# another model by setting STAT to "model_no" and specifying the 
# "best_index" directly.

# Specify STAT as either: waic or model_no
STAT <- "model_no"   
if (STAT == "model_no") {
  # If STAT == "model_no", specify model of interest by model number in table
  best_index <- 2
} else {
  min_stat <- min(model_stats[,STAT])
  best_index <- model_stats$model_no[model_stats[,STAT] == min_stat] 
}

# Extract output and formulas from best model in 
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

# If all the covariates have sufficient explanatory power, then move on to the 
# next step (looking at results from best model, below). If not, run a model 
# that excludes covariates with little to no explanatory power from the model 
# for occurrence. (One way to assess whether covariates have explanatory power
# is to identify those with f > 0.9)

  # Identify new set(s) of spatial covariates to include in model for inference:
  scov_new <- list(c("roads"), 
                   c("roads", "slope"), 
                   c("boundary", "slope"),
                   c("roads", "veg"),
                   c("roads", "veg", "slope"))
  OCC_MODELS <- lapply(scov_new, function(x) c(x, BEST_ANNUAL))
  # Refine the detection model (if needed)
  # DET_MODELS <- list(c("day", "effort"))
  source("src/multi-season-models/spOccupancy-MS-create-model-formulas.R")
  message("Check candidate models:", sep = "\n")
  model_specs
  
  # Run model(s)
  source("src/multi-season-models/spOccupancy-MS-run-candidate-models.R")
  model_stats %>% arrange(waic)

  # Specify STAT as either: waic or model_no
  STAT <- "model_no"   
  if (STAT == "model_no") {
    # If STAT == "model_no", specify model of interest by model number in table
    best_index <- 2  
  } else {
    min_stat <- min(model_stats[,STAT])
    best_index <- model_stats$model_no[model_stats[,STAT] == min_stat] 
  }
    
  # Extract output and formulas from best model in 
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
  
  # Check if ESS/Rhats are sufficient. If not, might need to run for longer
  summary(best)
  
  # If necessary, increase number of samples and run again:
    # N_CHAINS <- 3
    # N_BURN <- 4000
    # N_THIN <- 20
    # N_BATCH <- 560
    # BATCH_LENGTH <- 25
    # n_samples <- N_BATCH * BATCH_LENGTH
    # message("MCMC specifications will result in a total of ",
    #         (n_samples - N_BURN) * N_CHAINS / N_THIN,
    #         " posterior samples for each parameter.")
    #  
    # source("src/multi-season-models/spOccupancy-MS-run-candidate-models.R")
    # STAT <- "waic"   
    # min_stat <- min(model_stats[,STAT])
    # best_index <- model_stats$model_no[model_stats[,STAT] == min_stat] 
    # best <- out_list[[best_index]]
    # summary(best)
  
  # Save model object to file
  model_filename <- paste0("output/multi-season-models/", PARK, "-", 
                           YEARS[1], "-", YEARS[length(YEARS)],
                           "-", SPECIES, ".rds")
  model_list <- list(model = best, 
                     psi_model = best_psi_model,
                     p_model = best_p_model,
                     data = data_list)
  saveRDS(model_list, file = model_filename)  

#------------------------------------------------------------------------------#
# Evaluate best model and look at estimates
#------------------------------------------------------------------------------#

# Extract names of covariates (with and without "_z" subscripts) from best model
# And for occurrence, extract names of spatial covariates
psi_covs_z <- create_cov_list(best_psi_model)
if (length(psi_covs_z) == 1 & any(psi_covs_z == "1")) {
  psi_covs_z <- character(0)
  psi_covs <- character(0)
  psi_spatcovs_z <- character(0)
  psi_spatcovs <- character(0)
} else {
  psi_covs <- psi_covs_z %>% str_remove_all(pattern = "_z")
  psi_spatcovs_z <- psi_covs_z[!psi_covs_z %in% c("years_z", "visits_z", "traffic_z")]
  psi_spatcovs <- psi_covs[!psi_covs %in% c("years", "visits", "traffic")]  
}
p_covs_z <- create_cov_list(best_p_model)
if (length(p_covs_z) == 1 & any(p_covs_z == "1")) {
  p_covs_z <- character(0)
  p_covs <- character(0)
} else {
  p_covs <- p_covs_z %>% str_remove_all(pattern = "_z")
}

# View parameter estimates
summary(best)

# Can save estimates table to file with code below
# write.csv(estimates,
#           file = paste0("C:/.../",
#                         PARK, "-", SPECIES, "-", 
#                         YEARS[1], YEARS[length(YEARS)], ".csv"),
#           row.names = FALSE)

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
             "variation in occupancy and detection.\n"))
}
if (ppc.rep < 0.1 | ppc.site > 0.9) {
  warning(paste0("PPC indicates that we have not adequately described temporal ",
                 "variation in detection."))
} else {
  cat(paste0("PPC indicates that we have adequately described temporal ",
             "variation in detection.\n"))
}

#------------------------------------------------------------------------------#
# Calculate and visualize predicted probabilities of occurrence, across park
#------------------------------------------------------------------------------#

# Create spatial predictions IF there are spatial covariates in the occurrence
# part of the model. (Note that this can take several minutes to run.)

# If there are time-varying covariates (other than year/trend) in the occurrence 
# part of the model, identify whether we want predictions under average 
# conditions ("averaged") or under observed conditions in the first and last 
# year ("observed"). Note that if we're using "averaged" and years/trend isn't 
# in the model, then predictions from the first and last year will be very 
# similar (but not identical if we're incorporating random effects).
if (any(str_detect(string = psi_covs, 
                   pattern = paste(c("visits", "traffic", "monsoon_ppt", "ppt10"),
                                   collapse = "|")))) {
  ANN_PREDS <- "observed"
}  

if (length(psi_spatcovs) > 0) {
  source("src/multi-season-models/spOccupancy-MS-predictions.R")
  
  # This script creates:
  # best_pred: a list with predictions for each raster cell, MCMC sample
  # preds_mn_firstyr: a raster with mean occurrence probability in each cell 
  # in the first year (across MCMC samples)
  # preds_mn_lastyr: a raster with mean occurrence probability in each cell 
  # in the last year (across MCMC samples)
  # preds_sd_firstyr: a raster with SDs in each cell in the first year
  # preds_sd_lastyr: a raster with SDs in each cell in the last year
  # plot_preds_mn_firstyr: a ggplot object with predicted mean values across 
  # park in the first year
  # plot_preds_mn_lastyr: a ggplot object with predicted mean values across 
  # park in the last year
  # plot_preds_sd_firstyr: a ggplot object with predicted SD values across park 
  # in the first year
  # plot_preds_sd_lastyr: a ggplot object with predicted SD values across park 
  # in the last year
  
  # Plot predicted means in first, last year
  # grid.arrange(plot_preds_mn_firstyr, plot_preds_mn_lastyr, nrow = 2)
  
  # Plot predicted SDs in first, last year
  # grid.arrange(plot_preds_sd_firstyr, plot_preds_sd_lastyr, nrow = 2)
  
  # Can save any of the plots to file (example below):
  plot_save <- plot_preds_mn_lastyr +
    theme_bw(base_size = 8)
  plotname <- paste0("C:/Users/erin/OneDrive/SODN/Mammals/SAGW_20172022_PrelimResults/",
                     PARK, "-", SPECIES, "-OccProbMN-",
                     YEARS[length(YEARS)], ".jpg")
  plot_save1 <- plot_preds_mn_firstyr +
    theme_bw(base_size = 8)
  plotname1 <- paste0("C:/Users/erin/OneDrive/SODN/Mammals/SAGW_20172022_PrelimResults/",
                      PARK, "-", SPECIES, "-OccProbMN-",
                      YEARS[1], ".jpg")
  
  # ggsave(filename = plotname1,
  #        plot = plot_save1,
  #        device = "jpeg",
  #        width = 4,
  #        height = 4,
  #        units = "in",
  #        dpi = 600)
}

#------------------------------------------------------------------------------#
# Make this a plot of mean occurrence probability over time (for years or any
# other annual covariate). Include naive occupancy AND include random effects
# (with different symbols) #####################################################

# Calculate and create figures trends in occurrence probability over time
# (only if "years" is in the model for occurrence)
#------------------------------------------------------------------------------#

# Note: if there are other annual covariates in the model (eg, traffic), this is 
# the predicted trend assuming mean levels of that covariate each year. In other
# words, this is the predicted trend after accounting for all other covariates.

if ("years" %in% psi_covs) {
  trend <- trend_plot_occ(model = best, 
                          data_list = data_list,
                          covariate_table = covariates,
                          raw_occ = TRUE,
                          central_meas = mean,
                          lower_ci = 0.025,
                          upper_ci = 0.975)
  print(trend)
  # Save to file
  plot_save <- trend +
    theme_classic(base_size = 8)
  # plotname <- paste0("C:/.../",
  #                    PARK, "-", SPECIES, "-Trend.jpg")
  
  # ggsave(filename = plotname,
  #        plot = plot_save,
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
# Excluding years (trend) since that was covered in section above. 
psi_continuous <- psi_covs_z[!psi_covs_z %in% c("vegclass2", "vegclass3", "years_z")]
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
# ggsave(filename = paste0("C:/Users/erin/OneDrive/SODN/Mammals/SAGW_20172022_PrelimResults/",
#                          PARK, "-", SPECIES, "-Occ-wash.jpg"),
#        plot = marginal_psi_wash,
#        device = "jpeg",
#        width = 4,
#        height = 4,
#        units = "in",
#        dpi = 600)

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
p_continuous <- p_covs_z[!p_covs_z %in% c("vegclass2", "vegclass3", 
                                          "camera", "lens_2023")]
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
# ggsave(filename = paste0("C:/Users/erin/OneDrive/SODN/Mammals/SAGW_20172022_PrelimResults/",
#                          PARK, "-", SPECIES, "-Det-effort.jpg"),
#        plot = marginal_p_effort,
#        device = "jpeg",
#        width = 4,
#        height = 4,
#        units = "in",
#        dpi = 600)

# If vegetation classes were included as covariates in the model, extract
# detection probabilities for each class
if (sum(str_detect(p_covs, "veg")) > 0) {
  detprobs_veg <- vegclass_estimates(model = best, 
                                     parameter = "det")
  print(detprobs_veg)
}

# If camera and/or lens was included as a covariate in the model, extract 
# detection probabilities for each combination of covariate levels
if (sum(str_detect(p_covs, c("camera|lens_2023"))) > 0) {
  detprob_cat <- det_cat_estimates(model = best,
                                   lower_ci = 0.025,
                                   upper_ci = 0.975)
  print(detprob_cat)
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

#------------------------------------------------------------------------------#
# Clean up
#------------------------------------------------------------------------------#

# Remove weather rasters from local repo
invisible(file.remove(list.files(weather_folder, full.names = TRUE)))
