################################################################################
# Calculate predicted occurrence probabilities across a park in the first and 
# last year of the study

# In most instances, this will be called from another script (something like:
# src/multi-season-models/PARK/spOccupancy-PARK-SPECIES_YEARS.R)

# ER Zylstra
# Updated 2023-04-06
################################################################################

# Objects created in this script:
  # best_pred is a list with predictions for each raster cell, MCMC sample
  # preds_mn is a raster with mean values in each cell (across MCMC samples)
  # preds_sd is a raster with SDs in each cell (across MCMC samples)

#------------------------------------------------------------------------------#
# Create multi-layer raster with covariate data
#------------------------------------------------------------------------------#

# Extract layers from spat_raster for spatial covariates in best model 
# (excluding non-spatial covariates like year, traffic, and visits if present)
psi_rasters <- spat_raster[[psi_spatcovs]]

# Standardize values, where needed
zs <- which(str_detect(psi_spatcovs_z, "_z"))
for (z in zs) {
  var_name <- psi_spatcovs[z]
  var_mn <- mean(spatial_covs[, var_name])
  var_sd <- sd(spatial_covs[, var_name])
  psi_rasters[[z]] <- (psi_rasters[[z]] - var_mn) / var_sd
}

# Create quadratics, where needed
quads <- which(duplicated(psi_spatcovs))
for (q in quads) {
  psi_rasters[[q]] <- psi_rasters[[q]] * psi_rasters[[q]]
  names(psi_rasters[[q]]) <- paste0(names(psi_rasters[[q]]), "2")
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

# Predict occurrence in the first and last year of the dataset
pred_years <- YEARS[c(1, length(YEARS))]

# Create 3-dimensional array of data with dimensions: 1 = number of sites; 2 = 
# number of prediction years; 3 = number of predictors, including the intercept
X.0 <- array (1, dim = c(nrow(psi_rasters_df), 
                         length(pred_years), 
                         ncol(best$beta.samples)))

cov_order <- occ_estimates$Covariate

######### TODO: Adapt next section to deal with other annual covariates (eg, visits, traffic)

if (length(cov_order) > 1) {
  if ("years_z" %in% cov_order) {
    year_pred <- matrix(rep((pred_years - mean(data_list$occ.covs$years)) / 
                              sd(data_list$occ.covs$years),
                            nrow(psi_rasters_df)),
                        nrow = nrow(psi_rasters_df), ncol = length(pred_years),
                        byrow = TRUE)
    X.0[, , length(cov_order)] <- year_pred
  }  
  for (i in 2:ncol(psi_rasters_df)) {
    X.0[, , i] <- psi_rasters_df[, i]
  }
}

# Make predictions with predict.PGOcc() 
  # t.cols: index for which years you want predictions for (of only those years
  # we have data for)
  # ignore.RE: specify whether or not we want to ignore unstructured random 
  # effects in the prediction and just use the fixed effects and any structured 
  # random effects (ignore.RE = TRUE), or include unstructured random effects for 
  # prediction (ignore.RE = FALSE). When ignore.RE = FALSE, the estimated values 
  # of the unstructured random effects are included in the prediction for both 
  # sampled and unsampled sites. For sampled sites, these effects come directly 
  # from those estimated from the model, whereas for unsampled sites, the effects 
  # are drawn from a normal distribution using our estimates of the random effect 
  # variance. Including unstructured random effects in the predictions will 
  # generally improve prediction at sampled sites, and will lead to nearly 
  # identical point estimates at non-sampled sites, but with larger uncertainty.

ignore.RE <- FALSE
# If we want to include unstructured REs in predictions, we need to add a slice
# to the X.0 array (putting 0's in there, but I also tried NAs and the results
# were the same)
if (ignore.RE == FALSE) {
  RE <- matrix(0, nrow = dim(X.0)[1], ncol = dim(X.0)[2])
  X.0 <- abind(X.0, RE, along = 3)
  dimnames(X.0)[[3]] <- c(cov_order, "site")
}

best_pred <- predict(object = best, 
                     X.0 = X.0,
                     t.cols = which(YEARs %in% pred_years),
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

# Use tidyterra to create plots using ggplot syntax
plot_preds_mn_firstyr <- ggplot() + 
  geom_spatraster(data = preds_mn_firstyr, mapping = aes(fill = mean_firstyr)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = "", title = paste0("Mean occurrence probability, ", YEARS[1])) +
  theme_bw()
plot_preds_mn_lastyr <- ggplot() + 
  geom_spatraster(data = preds_mn_lastyr, mapping = aes(fill = mean_lastyr)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = "", title = paste0("Mean occurrence probability, ", YEARS[length(YEARS)])) +
  theme_bw()
# grid.arrange(plot_preds_mn_firstyr, plot_preds_mn_lastyr, nrow = 2)

plot_preds_sd_firstyr <- ggplot() + 
  geom_spatraster(data = preds_sd_firstyr, mapping = aes(fill = sd_firstyr)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = "", title = paste0("SD occurrence probability, ", YEARS[1])) +
  theme_bw()
plot_preds_sd_lastyr <- ggplot() + 
  geom_spatraster(data = preds_sd_lastyr, mapping = aes(fill = sd_lastyr)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = "", title = paste0("SD occurrence probability, ", YEARS[length(YEARS)])) +
  theme_bw()
# grid.arrange(plot_preds_sd_firstyr, plot_preds_sd_lastyr, nrow = 2)
