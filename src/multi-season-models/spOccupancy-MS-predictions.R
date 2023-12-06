################################################################################
# Calculate predicted occurrence probabilities across a park in the first and 
# last year of the study

# In most instances, this will be called from another script (something like:
# src/multi-season-models/PARK/spOccupancy-PARK-SPECIES_YEARS.R)

# ER Zylstra
# Updated 2023-04-07
################################################################################

# Objects created in this script:
  # best_pred is a list with predictions for each raster cell, MCMC sample
  # preds_mn is a raster with mean values in each cell (across MCMC samples)
  # preds_sd is a raster with SDs in each cell (across MCMC samples)

#------------------------------------------------------------------------------#
# Create multi-layer raster with covariate data
#------------------------------------------------------------------------------#

# Identify time-invariant spatial covariates (For the time being, all of our 
# spatial covariates are time-invariant. Leaving this in here in case we 
# add covariates that do vary over space and time)
spattime_covs <- NA
time_invar <- psi_spatcovs[!psi_spatcovs %in% spattime_covs]
time_invar_z <- psi_spatcovs_z[!psi_spatcovs_z %in% spattime_covs]

# Extract layers from park_raster for spatial covariates (that don't vary by
# year) in best model 
if (length(time_invar) > 0) {
  psi_rasters <- park_raster[[time_invar]]
  
  # Standardize values, where needed
  zs <- which(str_detect(time_invar_z, "_z"))
  for (z in zs) {
    var_name <- time_invar[z]
    var_mn <- mean(spatial_covs[, var_name])
    var_sd <- sd(spatial_covs[, var_name])
    psi_rasters[[z]] <- (psi_rasters[[z]] - var_mn) / var_sd
  }
  
  # Create quadratics, where needed
  quads <- which(duplicated(time_invar))
  for (q in quads) {
    psi_rasters[[q]] <- psi_rasters[[q]] * psi_rasters[[q]]
    names(psi_rasters[[q]]) <- paste0(names(psi_rasters[[q]]), "2")
  }
} else {
  # If there are only time-varying spatial covariates, create a raster with 
  # desired geometry (we won't need the actual values)
  psi_rasters <- park_raster[["elev"]]
}

# Crop and mask by park boundary
park_boundaries <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
park_boundary <- subset(park_boundaries, park_boundaries$UNIT_CODE == PARK)
psi_rasters <- crop(psi_rasters, park_boundary)
psi_rasters <- mask(psi_rasters, park_boundary)

#------------------------------------------------------------------------------#
# Predict psi (probability of occupancy) values across the park
#------------------------------------------------------------------------------#

# Convert raster to a dataframe (one row for each cell, each column = layer)
psi_rasters_df <- as.data.frame(psi_rasters, cell = TRUE)

# Remove rows that have any NA covariate values
if (ncol(psi_rasters_df) == 2) {
  psi_rasters_df$nNAs <- 1*is.na(psi_rasters_df[,2])
} else {
  psi_rasters_df$nNAs <- apply(psi_rasters_df[, -1], 1, function(x) sum(is.na(x)))
}
psi_rasters_df <- psi_rasters_df %>%
  dplyr::filter(nNAs == 0) %>%
  select(-nNAs)

# Select years for prediction (will usually select first and last year of study)
pred_years <- YEARS[c(1, length(YEARS))]

# Create 3-dimensional array of data with dimensions: 1 = number of sites; 2 = 
# number of prediction years; 3 = number of predictors, including the intercept.
# (Remember, we're only running this script if we have at >=1 spatial covariate)
X.0 <- array(NA, dim = c(nrow(psi_rasters_df), 
                         length(pred_years), 
                         ncol(best$beta.samples)))
# Make first slice = 1 for the intercept
X.0[, , 1] <- 1

# Extract order of all covariates in the model
cov_order <- occ_estimates$Covariate

# Grab standardized values of annual, non-spatial covariates for prediction 
# years and place them in the appropriate slice of X.0
if ("years_z" %in% cov_order) {
  year_pred <- matrix(rep((pred_years - mean(data_list$occ.covs$years)) / 
                            sd(data_list$occ.covs$years),
                          nrow(psi_rasters_df)),
                      nrow = nrow(psi_rasters_df), ncol = length(pred_years),
                      byrow = TRUE)
  X.0[, , which(cov_order == "years_z")] <- year_pred
}
if ("traffic_z" %in% cov_order) {
  if (ANN_PREDS == "observed") {
    traffic_pred <- matrix(rep(data_list$occ.covs$traffic_z[1, which(YEARS %in% pred_years)],
                               nrow(psi_rasters_df)),
                        nrow = nrow(psi_rasters_df), ncol = length(pred_years),
                        byrow = TRUE)
    X.0[, , which(cov_order == "traffic_z")] <- traffic_pred
  } else {
    X.0[, , which(cov_order == "traffic_z")] <- 0
  }
}
if ("visits_z" %in% cov_order) {
  if (ANN_PREDS == "observed") {
    visits_pred <- matrix(rep(data_list$occ.covs$visits_z[1, which(YEARS %in% pred_years)],
                              nrow(psi_rasters_df)),
                           nrow = nrow(psi_rasters_df), ncol = length(pred_years),
                           byrow = TRUE)
    X.0[, , which(cov_order == "visits_z")] <- visits_pred
  } else {
    X.0[, , which(cov_order == "visits_z")] <- 0
  }
}
if ("monsoon_ppt_z" %in% cov_order) {
  if (ANN_PREDS == "observed") {
    monsoon_pred <- matrix(rep(data_list$occ.covs$monsoon_ppt_z[1, which(YEARS %in% pred_years)],
                               nrow(psi_rasters_df)),
                           nrow = nrow(psi_rasters_df), ncol = length(pred_years),
                           byrow = TRUE)
    X.0[, , which(cov_order == "monsoon_ppt_z")] <- monsoon_pred
  } else {
    X.0[, , which(cov_order == "monsoon_ppt_z")] <- 0
  }
}
if ("ppt10_z" %in% cov_order) {
  if (ANN_PREDS == "observed") {
    ppt10_pred <- matrix(rep(data_list$occ.covs$ppt10_z[1, which(YEARS %in% pred_years)],
                             nrow(psi_rasters_df)),
                         nrow = nrow(psi_rasters_df), ncol = length(pred_years),
                         byrow = TRUE)
    X.0[, , which(cov_order == "ppt10_z")] <- ppt10_pred
  } else {
    X.0[, , which(cov_order == "ppt10_z")] <- 0
  }
}

# If annual, spatial covariates are in the occurrence model and we want to make 
# predictions for the first and last year under observed conditions, then grab
# the standarized values and place them in the appropriate slice of X.0. (Again,
# leaving this in as an example in case we add covariates that vary over space
# and time)
  # if ("monsoon_ppt_z" %in% cov_order) {
  #   if (ANN_PREDS == "observed") {
  #     monsoon_raster <- rast(monsoon_list[which(YEARS %in% pred_years)])
  #     monsoon_raster <- resample(monsoon_raster, psi_rasters, method = "near")
  #     monsoon_df <- as.data.frame(monsoon_raster, cell = TRUE)
  #     monsoon_df <- monsoon_df[monsoon_df$cell %in% psi_rasters_df$cell,]
  #     monsoon_df[, -1] <- (monsoon_df[, -1] - monsoon_ppt_mn) / monsoon_ppt_sd
  #     X.0[, , which(cov_order == "monsoon_ppt_z")] <- as.matrix(monsoon_df[, -1])
  #   } else {
  #     X.0[, , which(cov_order == "monsoon_ppt_z")] <- 0
  #   }
  # }

# Identify which slices of X.0 haven't been filled in yet (if any). The number 
# of slices should equal the number of time-invariant spatial covariates in the 
# model
slicestofill <- which(is.na(X.0[1, 1, ]))

if (length(slicestofill) > 0) {
  if (length(slicestofill) != length(time_invar)) {
    message("Dimensions of X.0 inconsistent with the number of spatial covariates.", 
            " Did not finish creating X.0.")
  } else {
    # Fill in X.0 will values of spatial covariates from psi_rasters_df
    for (i in 1:length(time_invar)) {
      X.0[, , slicestofill[i]] <- psi_rasters_df[, i + 1]
    }
  }
}

# Make predictions with predict.PGOcc() 
  # t.cols: index for which years you want predictions for (of only those years
    # we have data for)
  # ignore.RE: specify whether or not we want to ignore unstructured random 
    # effects in the prediction and just use the fixed effects and any 
    # structured random effects (ignore.RE = TRUE), or include unstructured 
    # random effects for prediction (ignore.RE = FALSE). When ignore.RE = FALSE, 
    # the estimated values of the unstructured random effects are included in 
    # the prediction for both sampled and unsampled sites. For sampled sites, 
    # these effects come directly from those estimated from the model, whereas 
    # for unsampled sites, the effects are drawn from a normal distribution 
    # using our estimates of the random effect variance. Including unstructured 
    # random effects in the predictions will generally improve prediction at 
    # sampled sites, and will lead to nearly identical point estimates at 
    # non-sampled sites, but with larger uncertainty.

ignore.RE <- FALSE
# Are yearly REs in the model?
yrRE <- ifelse(dim(best$X.re)[3] == 2, 1, 0)

# If we want to include unstructured REs in predictions, we need to add slices
# to the X.0 array (putting 0s in there for site, or year if want "averaged" 
# annual effects)
  # Check level names for random effects:
  # best$re.level.names
  # Check summary of random effects for each site, year:
  # summary(best$beta.star.samples)
if (ignore.RE == FALSE) {
  siteRE <- matrix(0, nrow = dim(X.0)[1], ncol = dim(X.0)[2])
  if (yrRE == 1 & ANN_PREDS == "observed") {
    yearRE <- matrix(pred_years, byrow = TRUE,
                     nrow = dim(X.0)[1], ncol = dim(X.0)[2])
    X.0 <- abind(X.0, siteRE, yearRE, along = 3)
    dimnames(X.0)[[3]] <- c(cov_order, "site", "years")
  }
  if (yrRE == 1 & ANN_PREDS == "averaged") {
    yearRE <- matrix(0, nrow = dim(X.0)[1], ncol = dim(X.0)[2])
    X.0 <- abind(X.0, siteRE, yearRE, along = 3)
    dimnames(X.0)[[3]] <- c(cov_order, "site", "years")
  }
  if (yrRE == 0) {
    X.0 <- abind(X.0, siteRE, along = 3)
    dimnames(X.0)[[3]] <- c(cov_order, "site")
  }
}

best_pred <- predict(object = best, 
                     X.0 = X.0,
                     t.cols = which(YEARS %in% pred_years),
                     ignore.RE = ignore.RE,
                     type = "occupancy")

# Plot predicted occupancy probability (can also take a couple minutes)
plot_dat <- data.frame(cell = psi_rasters_df$cell,
                       mean_psi_firstyr = apply(best_pred$psi.0.samples[,,1], 2, mean),
                       mean_psi_lastyr = apply(best_pred$psi.0.samples[,,2], 2, mean),
                       sd_psi_firstyr = apply(best_pred$psi.0.samples[,,1], 2, sd),
                       sd_psi_lastyr = apply(best_pred$psi.0.samples[,,2], 2, sd))

preds_mn_firstyr <- rast(psi_rasters[[1]])
preds_mn_firstyr[plot_dat[,1]] <- plot_dat[,2]
names(preds_mn_firstyr) <- "mean_firstyr"
preds_mn_lastyr <- rast(psi_rasters[[1]])
preds_mn_lastyr[plot_dat[,1]] <- plot_dat[,3]
names(preds_mn_lastyr) <- "mean_lastyr"

preds_sd_firstyr <- rast(psi_rasters[[1]])
preds_sd_firstyr[plot_dat[,1]] <- plot_dat[,4]
names(preds_sd_firstyr) <- "sd_firstyr"
preds_sd_lastyr <- rast(psi_rasters[[1]])
preds_sd_lastyr[plot_dat[,1]] <- plot_dat[,5]
names(preds_sd_lastyr) <- "sd_lastyr"

# Use tidyterra to create plots with the same color scale in both years
minmax_mn <- plot_dat %>%
  select(contains("mean_psi")) %>%
  as.matrix() %>%  
  range
minmax_sd <- plot_dat %>%
  select(contains("sd_psi")) %>%
  as.matrix() %>%
  range

col_scale_mn = scale_fill_gradientn(
  colors = hcl.colors(100, palette = "viridis"),
  limits = minmax_mn,
  na.value = "transparent"
)
col_scale_sd = scale_fill_gradientn(
  colors = hcl.colors(100, palette = "viridis"),
  limits = minmax_sd,
  na.value = "transparent"
)

plot_preds_mn_firstyr <- ggplot() + 
  geom_spatraster(data = preds_mn_firstyr, mapping = aes(fill = mean_firstyr)) + 
  col_scale_mn +
  labs(fill = "", title = paste0("Mean occurrence probability, ", YEARS[1])) +
  theme_bw()
plot_preds_mn_lastyr <- ggplot() + 
  geom_spatraster(data = preds_mn_lastyr, mapping = aes(fill = mean_lastyr)) + 
  col_scale_mn +
  labs(fill = "", title = paste0("Mean occurrence probability, ", YEARS[length(YEARS)])) +
  theme_bw()
# grid.arrange(plot_preds_mn_firstyr, plot_preds_mn_lastyr, nrow = 2)

plot_preds_sd_firstyr <- ggplot() + 
  geom_spatraster(data = preds_sd_firstyr, mapping = aes(fill = sd_firstyr)) + 
  col_scale_sd +
  labs(fill = "", title = paste0("SD occurrence probability, ", YEARS[1])) +
  theme_bw()
plot_preds_sd_lastyr <- ggplot() + 
  geom_spatraster(data = preds_sd_lastyr, mapping = aes(fill = sd_lastyr)) + 
  col_scale_sd +
  labs(fill = "", title = paste0("SD occurrence probability, ", YEARS[length(YEARS)])) +
  theme_bw()
# grid.arrange(plot_preds_sd_firstyr, plot_preds_sd_lastyr, nrow = 2)
