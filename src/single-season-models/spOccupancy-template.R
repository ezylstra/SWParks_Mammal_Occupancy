################################################################################
# Template to run and evaluate a suite of single-season occupancy models for 
# a given park, year, and species (using the spOccupancy package)

# Objects that need to be specified to run models (ie, to create 
# src/single-season-model/YEAR/spOccupancy-PARK-SPECIES-YEAR.R)
  # PARK, YEAR, SPECIES (lines 41, 44, 58)
  # Specifications for occurrence models: OCC_NULL, OCC_MODELS1, OCC_MODELS2
  # Specifications for detection models: DET_NULL, DET_MODELS1, DET_MODELS2
  # Statistic used to select a "best" model for inferences: STAT

# ER Zylstra
# Updated 2023-02-01
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
  
  occm1 <- covariates$formula[covariates$short_name %in% OCC_MODELS1]
  occm2 <- list()
  for (i in 1:length(OCC_MODELS2)) {
    occm2[[i]] <- paste(covariates$formula[covariates$short_name %in% OCC_MODELS2[[i]]],
                        collapse = " + ")
  }
  occ_specs <- c(occm1, unlist(occm2))
  if (OCC_NULL) {
   occ_specs <- c("1", occ_specs) 
  }
  occ_specs <- paste0("~ ", occ_specs) 
  message("Check that these are the candidate models of interest for occurrence.")
  occ_specs

# Logical indicating whether a null model should be included in a candidate set
  DET_NULL <- TRUE

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
  
  detm1 <- covariates$formula[covariates$short_name %in% DET_MODELS1]
  detm2 <- list()
  for (i in 1:length(DET_MODELS2)) {
    detm2[[i]] <- paste(covariates$formula[covariates$short_name %in% DET_MODELS2[[i]]],
                        collapse = " + ")
  }
  det_specs <- c(detm1, unlist(detm2))
  if (DET_NULL) {
    det_specs <- c("1", det_specs) 
  }
  det_specs <- paste0("~ ", det_specs) 
  message("Check that these are the candidate models of interest for detection.")
  det_specs

# Create a matrix that contains all combinations of occurrence and detection 
# covariates for candidate models
model_specs <- as.matrix(expand.grid(occ = occ_specs, 
                                     det = det_specs,
                                     KEEP.OUT.ATTRS = FALSE))

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
set.seed(2023)
out_list <- list()
for (i in 1:nrow(model_specs)) {
  message("Running model ", i, " (of ", nrow(model_specs), ").")
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
               k.fold.threads = 4,
               k.fold.seed = 2023) 
  out_list <- c(out_list, list(out))
}

# Gather results from models (occ/det formulas, Bayesian p-values from posterior 
# predictive checks, WAIC, deviance stat from 4-fold CV)
source("src/single-season-models/spOccupancy-model-stats.R")
  # Takes a minute to run with posterior predictive checks (PPCs)

# View summary table, ranked by WAIC
model_stats %>%
  arrange(waic)
# View summary table, ranked by deviation statistic from k-fold CV
model_stats %>%
  arrange(k.fold.dev)

#------------------------------------------------------------------------------#
# Look at results and predictions from "best" model
#------------------------------------------------------------------------------#

# Identify statistic to use for selecting the best model 
# Either WAIC (waic) or deviance from k-fold CV (k.fold.dev)
# STAT <- "waic"
STAT <- "k.fold.dev"
min_stat <- min(model_stats[,stat])
best_index <- model_stats$model_no[model_stats[,stat] == min_stat] 
best <- out_list[[best_index]]

# View covariate structure
best_psi_model <- model_specs[best_index, 1]
best_p_model <- model_specs[best_index, 2]
message("psi ", best_psi_model)
message("p ", best_p_model)

# Parameter estimates
summary(best)
# Note: this is good for viewing, but will want to use other means to create
# a table for reports/publications

# Trace plots
plot(best$beta.samples, density = FALSE)
plot(best$alpha.samples, density = FALSE)
par(mfrow = c(1,1))

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

#------------------------------------------------------------------------------#
# Calculate predicted probability of occupancy, across park
#------------------------------------------------------------------------------#
# Everything in this section can probably be moved to source script

# Extract layers from spat_raster for covariates in best model
psi_covs_z <- best_psi_model %>%
  str_remove(pattern = "~ ") %>% 
  str_remove_all(pattern = "I[(]") %>%
  str_remove_all(pattern = "[)]") %>%
  str_remove_all(pattern = "\\^2") %>%
  str_split_1(pattern = " [+] ")
psi_covs <- psi_covs_z %>%
  str_remove_all(pattern = "_z")

psi_rasters <- spat_raster[[psi_covs]]

# Standardize values, where needed
zs <- which(str_detect(psi_covs_z, "_z"))
for (z in zs) {
  var_name <- psi_covs[z]
  var_mn <- mean(spatial_covs[, var_name])
  var_sd <- sd(spatial_covs[, var_name])
  psi_rasters[[z]] <- (psi_rasters[[z]] - var_mn) / var_sd
}
# Create quadratics, where needed
quads <- which(duplicated(psi_covs))
for (q in quads) {
  psi_rasters[[q]] <- psi_rasters[[q]] * psi_rasters[[q]]
  names(psi_rasters[[q]]) <- paste0(names(psi_rasters[[q]]), "2")
}
# Add layer for intercept
rast_int <- rast(psi_rasters[[1]])
rast_int <- rast(rast_int, vals = 1)
names(rast_int) <- "int"
psi_rasters <- c(rast_int, psi_rasters)

# Crop and mask by park boundary
park_boundaries <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
park_boundary <- subset(park_boundaries, park_boundaries$UNIT_CODE == PARK)
psi_rasters <- crop(psi_rasters, park_boundary)
psi_rasters <- mask(psi_rasters, park_boundary)

# Convert to a dataframe (one row for each cell, each column = layer)
psi_rasters_df <- as.data.frame(psi_rasters, cell = TRUE)
# Remove row with any covariates equal to NA
psi_rasters_df$nNAs <- apply(psi_rasters_df[, -1], 1, function(x) sum(is.na(x)))
psi_rasters_df <- psi_rasters_df %>%
  dplyr::filter(nNAs == 0) %>%
  select(-nNAs)

# Make predictions with predict.PGOcc() 
# Note: this can take a couple minutes, but is still faster than doing math
best_pred <- predict(best, as.matrix(psi_rasters_df[, -1]))

# Plot predicted occupancy probability (can also take a couple minutes)
plot_dat <- data.frame(cell = psi_rasters_df$cell,
                       mean_psi = apply(best_pred$psi.0.samples, 2, mean),
                       sd_psi = apply(best_pred$psi.0.samples, 2, sd))

preds_mn <- rast(psi_rasters[[1]])
preds_mn[plot_dat[,1]] <- plot_dat[,2]
names(preds_mn) <- "mean"
preds_sd <- rast(psi_rasters[[1]])
preds_sd[plot_dat[,1]] <- plot_dat[,3]
names(preds_sd) <- "sd"

# Could use default plot
# plot(preds_mn)

# Or use tidyterra to plot using ggplot syntax
ggplot() + 
  geom_spatraster(data = preds_mn, mapping = aes(fill = mean)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = '', title = 'Mean occurrence probability') +
  theme_bw()

ggplot() + 
  geom_spatraster(data = preds_sd, mapping = aes(fill = sd)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = '', title = 'SD occurrence probability') +
  theme_bw()

#------------------------------------------------------------------------------#
# Create figures depicting marginal effects of covariates on occurrence 
# probability (predicted covariate effects assuming all other covariates held
# constant)
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
