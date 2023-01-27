################################################################################
# Template to run and evaluate a suite of single-season occupancy models for 
# a given park, year, and species (using the spOccupancy package)

# Objects that need to be specified to run models for a particular park,
# year, and species 
# (ie, to create src/single-season-model/YEAR/spOccupancy-PARK-SPECIES-YEAR.R)
  # PARK, YEAR, SPECIES (lines 38, 41, 55)
  # OCC_MODELS, DET_MODELS (starting on line 119)

# ER Zylstra
# Updated 2023-01-26
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(spOccupancy)

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
  filter(Park == PARK & yr == YEAR) %>%
  arrange(desc(propdetect))
# View just those species with a detection rate of 5% 
# (n = camera location * sampling occasion)
detects %>% 
  filter(propdetect >= 0.05) %>%
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
  filter(abs(corr) > 0.5)

# Create names of spatial covariates or covariate groups
# (This helps us avoid having to include "_z" in the variable names, and 
# groups certain variables that will always appear together in a model)
boundary <- "boundary_z"
aspect <- "east_z + north_z"
elevation <- "elev_z"
elevation2 <- "elev_z + I(elev_z ^ 2)"
pois <- "pois_z"
roads <- "roads_z"
slope <- "slope_z"
slope2 <- "slope_z + I(slope_z ^ 2)"
trail <- "trail_z"
if (PARK == "CHIR") {
  burn <- "burn_severity"
}
if (PARK == "SAGW") {
  wash <- "wash_z"
  veg <- "vegclass2 + vegclass3"
}

# Create names of survey covariates
day <- "day_z"
day2 <- "day_z + I(day_z^2)"
deploy <- "deploy_exp"
effort <- "effort_z"

# We're going to run multiple models and store the output in a list.
# To make this easy, we'll create a list of occupancy formulas (specifications
# for covariates in the occurrence part of a model) and a list of detection
# formulas (specifications for covariates in the detection part of the model)

# Create a vector with every combination of covariates you'd like in the
# occurrence part of the model. 
  # Use the following syntax when you'd like to include more than one covariate 
  # (or covariate group): paste(c(variable1, variable2), collapse = " + ")
OCC_MODELS <- c(boundary,
                aspect,
                roads,
                slope,
                wash,
                veg,
                paste(c(veg, boundary), collapse = " + "),
                paste(c(veg, roads), collapse = " + "), 
                paste(c(veg, wash, roads), collapse = " + ")) 

# Create a vector with every combination of covariates you'd like in the
# detection part of the model. 
  # Use the following syntax when you'd like to include more than one covariate 
  # (or covariate group): paste(c(variable1, variable2), collapse = " + ")
DET_MODELS <- c(effort, 
                paste(c(day2, deploy, effort), collapse = " + ")) 

# Create a matrix that contains all combinations for models for occurrence and
# detection
model_specs <- as.matrix(expand.grid(occ = OCC_MODELS, 
                                     det = DET_MODELS,
                                     KEEP.OUT.ATTRS = FALSE))
model_specs <- apply(model_specs, 1:2, function(x) paste0("~", x))

# Create model formulas with R syntax
as.formula.vect <- Vectorize(as.formula)
occ_formulas <- as.formula.vect(model_specs[,1])
det_formulas <- as.formula.vect(model_specs[,2])

#------------------------------------------------------------------------------#
# Run models and compare fit
#------------------------------------------------------------------------------#

# Set MCMC parameters
n_samples <- 5000
n_burn <- 3000
n_thin <- 1
n_chains <- 3

# Running a list of models, doing 4-fold cross validation for each model
  # Not specifying priors, but using defaults which are N(0, var = 2.72)
  # Not specifying initial values -- by default they come from priors
  # Running chains sequentially (n.omp.threads = 1) because vignette states
  # this only speeds things up in spatial models
out_list <- list()
for (i in 1:nrow(model_specs)) {
  out <- PGOcc(occ.formula = occ_formulas[[i]],
               det.formula = det_formulas[[i]], 
               data = data_list, 
               # inits = inits, 
               # priors = priors, 
               n.samples = n_samples, 
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = n_burn, 
               n.thin = n_thin, 
               n.chains = n_chains,
               k.fold = 4,
               k.fold.threads = 4) 
  out_list <- c(out_list, list(out))
}

# Gather results from models (occ/det formulas, Bayesian p-values from posterior 
# predictive checks, WAIC, deviance stat from 4-fold CV)
source("src/single-season-models/spOccupancy-model-stats.R")

# View summary table, with top-ranked WAIC model at top
model_stats %>%
  arrange(waic)
# View summary table, with top-ranked k-fold deviation model at top
model_stats %>%
  arrange(k.fold.dev)

#------------------------------------------------------------------------------#
# Look at results and predictions from "best" model
#------------------------------------------------------------------------------#

# Identify statistic to use for selecting the best model 
# Either WAIC (waic) or deviance from k-fold CV (k.fold.dev)
stat <- "WAIC"
best_index <- model_stats$model_no[model_stats$waic == min(model_stats$waic)] 
best <- out_list[[best_index]]

# Parameter estimates
summary(best)
# Trace plots
plot(best$beta.samples, density = FALSE)
plot(best$alpha.samples, density = FALSE)
par(mfrow = c(1,1))

# Posterior predictive checks
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

# Plot the difference in the discrepancy measure between the replicate 
# and actual data across each of the sites (identify sites that are causing a 
# lack of fit).
best_ppcs <- ppc.sites[[best_index]]
diff_fit <- best_ppcs$fit.y.rep.group.quants[3, ] - best_ppcs$fit.y.group.quants[3, ]
plot(diff_fit, pch = 19, xlab = 'Site ID', ylab = "Replicate - True Discrepancy") 
prob_sites <- which(abs(diff_fit) > 0.4)
# Plot on map
plot(lat~long, data = spatial_covs, las = 1) # all camera locs
points(lat~long, data = spatial_covs[prob_sites,], pch = 19, col = "blue")

# TODO #########################################################
# Create a source script to produce maps with predicted occupancy
# Create a template or source script(s) to create figures with marginal
  # covariate effects
# Think about whether we want to save model output to file?  (maybe not if
  # these SS models don't take long to run?)

