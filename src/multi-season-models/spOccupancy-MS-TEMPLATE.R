################################################################################
# Template to run and evaluate a suite of multi-season occupancy models for 
# a given park, set of years, and species (using the spOccupancy package)

# Objects that need to be specified to run models (ie, to create 
# src/multi-season-model/PARK/spOccupancy-PARK-SPECIES-YEARS.R)

# ER Zylstra
# Updated 2023-04-07
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

source("src/functions.R")

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

# Select years of interest
YEARS <- 2017:2022

# Look at detection data for various species
detects <- read.csv("output/species-detections-bypark.csv", header = TRUE)
detects <- detects %>%
  dplyr::filter(Park == PARK) %>%
  arrange(desc(propdetect))
# View just those species with a detection rate of 5% (propdetect = proportion 
# of nobs [camera locations * sampling occasion] with species detection)
detects %>% 
  dplyr::filter(propdetect >= 0.05) %>%
  left_join(species, by = c("spp" = "Species_code")) %>%
  select(c(spp, Species, Common_name, nobs, propdetect))

# Select species of interest (ideally with a detection rate of at least 5%)
SPECIES <- "PETA"

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

#------------------------------------------------------------------------------#
# Specify the occurrence portion of candidate models
#------------------------------------------------------------------------------#

# Load dataframe with information about covariates:
covariates <- read.csv("data/covariates/covariates-MS.csv", header = TRUE)

# View those pairs of continuous covariates that are highly correlated
cor_df %>%
  arrange(desc(corr)) %>%
  dplyr::filter(abs(corr) > 0.5)

# See what covariates are available for occurrence part of model
covariates %>%
  dplyr::filter(parameter %in% c("either", "occupancy")) %>%
  dplyr::filter(park %in% c(PARK, "all")) %>%
  select(-c(parameter, park))

# Logical indicating whether a null model for occurrence should be included in 
# the candidate model set
OCC_NULL <- TRUE

# Pick covariates to include in simple candidate models via the short_name 
# column in the covariates dataframe. Note that including "years" as a 
# covariate creates a trend model (logit-linear trend in occurrence probability)
OCC_MODELS1 <- c("years", "visits")

# To combine covariates in a single candidate model, provide a vector of 
# short_names. Compile these vectors into a list.
# e.g., c("aspect", "boundary") would create the following model for occurrence: 
# psi ~ east + north + boundary
OCC_MODELS2 <- list(c("slope2", "years"),
                    # c("veg", "wash", "years"),
                    # c("elev2", "years"),
                    c("boundary", "years"),
                    # c("pois", "years"),
                    # c("roads", "years"),
                    # c("trail", "years"),
                    # c("roadbound", "years"),
                    # c("trailpoi", "years"),
                    c("slope2", "boundary", "years"),
                    c("traffic", "slope2", "years"))

#------------------------------------------------------------------------------#
# Specify the detection portion of candidate models
#------------------------------------------------------------------------------#

# See what covariates are available for detection part of model
covariates %>%
  dplyr::filter(parameter %in% c("either", "detection")) %>%
  dplyr::filter(park %in% c(PARK, "all")) %>%
  select(-c(parameter, park))
  
# Logical indicating whether a null model for detection should be included in 
# the candidate model set
DET_NULL <- FALSE

# Pick covariates to include in simple candidate models via the short_name 
# column in the covariates dataframe
DET_MODELS1 <- c("effort")

# To combine different covariates in a candidate model, provide a vector of 
# short_names (Note: not including camera_2022 in models since that seems to
# cause some problems, likely because that's the last year we have data for. 
# Random yearly effects might be more effective)
DET_MODELS2 <- list(c("day2", "deploy_exp", "effort"))

#------------------------------------------------------------------------------#
# Create (and check) formulas for candidate models
#------------------------------------------------------------------------------#

# Use OCC and DET objects to create formulas for candidate models

# Random effects in occurrence part of model
  # Note: for this stage of the analysis, it's worth considering whether we 
  # include unstructured random effects for both site and time. 
  # Site REs make sense in these implicit dynamic models because we want to 
  # account for non-independence of data from sites over multiple years. 
  # A yearly RE could also make sense because we don't have many covariates that 
  # can explain annual variation in occurrence (or temporal variation beyond 
  # that explained by a linear trend). However, we could run into problems with 
  # convergence if we include both spatial and temporal random effects given the
  # relatively short time series.

  # For now, we can start with an unstructured site effect in all candidate 
  # models, and then evaluate alternative random effect structures at a later 
  # stage. To add an unstructured site RE, we need to add (1 | site) to each 
  # formula for the occurrence part of the model. 
  
  # Specify TIME_RE_OCC and SITE_RE_OCC as either "none" or "unstructured"
  TIME_RE_OCC <- "none" 
  SITE_RE_OCC <- "unstructured"
  
# Random effects in detection part of model  
  # For now, we won't include any random effects
  TIME_RE_DET <- "none"
  SITE_RE_DET <- "none"

source("src/multi-season-models/spOccupancy-MS-create-model-formulas.R")

message("Check candidate models:", sep = "\n")
model_specs

#------------------------------------------------------------------------------#
# Run models and compare fit
#------------------------------------------------------------------------------#

# Set MCMC parameters

N_CHAINS <- 3
N_BURN <- 2000
N_THIN <- 10

# Instead of specifying the total number of samples (like we did for single-
# season models), we'll split the samples into a set of N_BATCH batches, each
# comprised of BATCH_LENGTH samples to improve mixing for the adaptive 
# algorithm. See documentation for the spOccupancy package for more info.

N_BATCH <- 280
BATCH_LENGTH <- 25

n_samples <- N_BATCH * BATCH_LENGTH
message("MCMC specifications will result in a total of ",
        (n_samples - N_BURN) * N_CHAINS / N_THIN, 
        " posterior samples for each parameter.")

# Run candidate models using spOccupancy package
source("src/multi-season-models/spOccupancy-MS-run-candidate-models.R")
  # Note: this can take several minutes to run

# View summary table, ranked by WAIC
  model_stats %>% arrange(waic)

# Description of columns in summary table:
  # psi: formula for occurrence part of model
  # det: formula for detection part of model
  # max.rhat: maximum value of R-hat across model parameters (want value < 1.05)
  # min.ESS: minimum value of ESS (effective sample size) across model
    # parameters (want value > 400)
  # [NOT IN THERE YET] ppc.sites: posterior predictive checks when binning the 
    # data across sites. P-values < 0.1 or > 0.9 can indicate that model fails 
    # to adequately represent variation in occurrence or detection across space.
  # [NOT IN THERE YET] ppc.reps: posterior predictive checks when binning the 
    # data across replicates. P-values < 0.1 or > 0.9 can indicate that model 
    # fails to adequately represent variation in detection over time.
  # waic: WAIC (Widely Applicable Information Criterion) for comparing models
    # (lower is better)

# Check that r-hat values and ESS look okay for most models.  If not, may 
# need to re-run after increasing N_BATCH or removing site REs. If problem
# persists, may need to remove some covariate combinations from consideration.
  
#------------------------------------------------------------------------------#
# Look at results and predictions from "best" model
#------------------------------------------------------------------------------#

# To look at the estimates from any model in the list, use the following,
# replacing "X" with the model_no of interest
# summary(out_list[[X]]) 

# Identify a model to use for inferences.  Can base this on WAIC or can select 
# any model by setting STAT equal to "model_no" and specifying the "best_index" 
# directly.

# Specify STAT as either: waic or model_no
STAT <- "model_no"   

if (STAT == "model_no") {
  # If STAT == "model_no", specify model of interest by model number in table
  best_index <- 4 
} else {
  min_stat <- min(model_stats[,STAT])
  best_index <- model_stats$model_no[model_stats[,STAT] == min_stat] 
}

# Extract output and formulas from best model
best <- out_list[[best_index]]
best_psi_model <- model_specs[best_index, 1]
best_p_model <- model_specs[best_index, 2]

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
# Create table with summary stats that can be saved to file
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

# Can save this table to file with code below
# write.csv(estimates,
#           file = paste0("C:/Users/erin/Desktop/Mammals/",
#                         PARK, "-", SPECIES, "-", 
#                         YEARS[1], YEARS[length(YEARS)], ".csv"),
#           row.names = FALSE)

# Trace plots
plot(best$beta.samples, density = FALSE)
plot(best$alpha.samples, density = FALSE)

# Posterior predictive checks (want Bayesian p-values between 0.1 and 0.9)
# ppc.site <- as.numeric(model_stats$ppc.sites[model_stats$model_no == best_index])
# ppc.rep <- as.numeric(model_stats$ppc.reps[model_stats$model_no == best_index])
# if (ppc.site < 0.1 | ppc.site > 0.9) {
#   warning(paste0("PPC indicates that we have not adequately described spatial ",
#                  "variation in occupancy and/or detection."))
# } else {
#   cat(paste0("PPC indicates that we have adequately described spatial ",
#              "variation in occupancy and detection."))
# } 
# if (ppc.rep < 0.1 | ppc.site > 0.9) {
#   warning(paste0("PPC indicates that we have not adequately described temporal ",
#                  "variation in detection."))
# } else {
#   cat(paste0("PPC indicates that we have adequately described temporal ",
#                  "variation in detection."))
# }

# If there's evidence that spatial variation isn't well explained, plot the 
# difference in the discrepancy measure between the replicate and actual data 
# across each of the sites (identify sites that are causing a lack of fit).
# if (ppc.site < 0.1 | ppc.site > 0.9) {
#   best_ppcs <- ppc.sites[[best_index]]
#   diff_fit <- best_ppcs$fit.y.rep.group.quants[3, ] - best_ppcs$fit.y.group.quants[3, ]
#   
#   # Plot differences
#   par(mfrow = c(1,1))
#   plot(diff_fit, pch = 19, xlab = 'Site ID', ylab = "Replicate - True Discrepancy") 
#   
#   # Identify sites on a map
#   prob_sites <- which(abs(diff_fit) > 0.4)
#   plot(lat~long, data = spatial_covs, las = 1) # all camera locs
#   points(lat~long, data = spatial_covs[prob_sites,], pch = 19, col = "blue")
# }

#------------------------------------------------------------------------------#
# Calculate and visualize predicted probabilities of occurrence, across park
#------------------------------------------------------------------------------#

# Create spatial predictions IF there are spatial covariates in the occurrence
# part of the model. 
# NOTE 1: this script will not work if there are covariates in the model that 
# vary over space and time (eg, precipitation variables). Need to add this.
# NOTE 2: this can take several minutes to run.

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
    grid.arrange(plot_preds_mn_firstyr, plot_preds_mn_lastyr, nrow = 2)
  
  # Plot predicted SDs in first, last year
    grid.arrange(plot_preds_sd_firstyr, plot_preds_sd_lastyr, nrow = 2)
  
  # Can save any of the plots to file (example below):
  plot_save <- plot_preds_mn_lastyr +
    theme_bw(base_size = 8)
  plot_save1 <- plot_preds_mn_firstyr +
    theme_bw(base_size = 8)
  plotname <- paste0("C:/Users/erin/Desktop/Mammals/",
                     PARK, "-", SPECIES, "-OccProbMN-BoundaryWashVeg-",
                     YEARS[length(YEARS)], ".jpg")
  # ggsave(filename = plotname,
  #        plot = plot_save,
  #        device = "jpeg",
  #        width = 4,
  #        height = 4,
  #        units = "in",
  #        dpi = 600)
}

#------------------------------------------------------------------------------#
# Calculate and create figures trends in occurrence probability over time
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
  trend
  # Save to file
  plot_save <- trend +
    theme_classic(base_size = 8)
  plotname <- paste0("C:/Users/erin/Desktop/Mammals/",
                     PARK, "-", SPECIES, "-Trend-(BoundaryWashVeg).jpg")
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

# If vegetation classes were included as covariates in the model, extract
# occurrence probabilities for each class
  if (sum(str_detect(psi_covs, "veg")) > 0) {
    occprobs_veg <- vegclass_estimates(model = best, 
                                       parameter = "occ")
    print(occprobs_veg)
  }

# If there are no covariates in the model (ie, a null model), print overall 
# occurrence probability
  if (psi_n_cont == 0 & length(psi_covs) == 1) {
    if (psi_covs == "1") {
      overall_occ <- mean_estimate(model = best, 
                                   parameter = "occ",
                                   lower_ci = 0.025,
                                   upper_ci = 0.975)
      print(overall_occ)
    }
  }  
  
#------------------------------------------------------------------------------#
# Create figures depicting marginal effects of covariates on detection 
# probability (predicted covariate effects assuming all other covariates held
# constant)
#------------------------------------------------------------------------------#

# Identify continuous covariates in occurrence part of the best model
p_continuous <- p_covs_z[!p_covs_z %in% c("vegclass2", "vegclass3")]
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
  if (p_n_cont == 0 & length(p_covs) == 1) {
    if (p_covs == "1") {
      overall_det <- mean_estimate(model = best, 
                                   parameter = "det",
                                   lower_ci = 0.025,
                                   upper_ci = 0.975)
      print(overall_det)
    }
  }  
