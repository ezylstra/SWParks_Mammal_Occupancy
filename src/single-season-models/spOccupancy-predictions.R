################################################################################
# Calculate predicted probability of occupancy across a park

# In most instances, this will be called from another script (something like:
# src/single-seasons/models/YEAR/spOccupancy-PARK-SPECIES_YEAR.R)

# ER Zylstra
# Updated 2023-02-02
################################################################################

# Objects created in this script:
  # best_pred is a list with predictions for each raster cell, MCMC sample
  # preds_mn is a raster with mean values in each cell (across MCMC samples)
  # preds_sd is a raster with SDs in each cell (across MCMC samples)

#------------------------------------------------------------------------------#
# Create multi-layer raster with covariate data
#------------------------------------------------------------------------------#

# Extract layers from spat_raster for covariates in best model
psi_rasters <- park_raster[[psi_covs]]

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

# Make predictions with predict.PGOcc() 
# Note: this can take a couple minutes, but is still faster than doing math
X.0 <- as.matrix(psi_rasters_df[, -1])
# Make first column = 1 for the intercept
X.0 <- cbind(1, X.0)

# If site REs are in the occurrence model, we will need to decide if we want
# them incorporated into predictions
  
  # Look at RE level names
  # best$re.level.names
  
  # Logical indicating whether to ignore random effects or not (FALSE means
  # that predictions will incorporate random effects)
  ignore.RE <- FALSE
  
  if ("site" %in% colnames(best$X.re) & !ignore.RE) {
    # site = 0 means that random effects will be generated for each raster cell
    # from a normal distribution with the estimated SD
    X.0 <- cbind(X.0, site = 0)
  }

# Generate predictions
best_pred <- predict(object = best, X.0 = X.0, ignore.RE = ignore.RE)

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

# Use tidyterra to create plots using ggplot syntax
plot_preds_mn <- ggplot() + 
  geom_spatraster(data = preds_mn, mapping = aes(fill = mean)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = '', title = 'Mean occurrence probability') +
  theme_bw()

plot_preds_sd <- ggplot() + 
  geom_spatraster(data = preds_sd, mapping = aes(fill = sd)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(fill = '', title = 'SD occurrence probability') +
  theme_bw()
