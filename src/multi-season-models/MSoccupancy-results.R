################################################################################
# Process results from a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-11-18
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(jagsUI)
library(MCMCvis)
library(ggplot2)

#------------------------------------------------------------------------------#
# Load the workspace that contains the JAGS model
#------------------------------------------------------------------------------#

PARK <- "SAGW"
SPECIES <- "LECA"
DATE <- "2022-11-18"
output_file <- paste0("output/models/",
                      tolower(PARK), "-",
                      tolower(SPECIES), "-MS-",
                      DATE,
                      ".Rdata")
load(output_file)

#------------------------------------------------------------------------------#
# First look at model results 
#------------------------------------------------------------------------------#

print(out, digits = 2)

# Trace and density plots
# MCMCtrace(out,pdf = FALSE)
# par(mfrow = c(1,1))

#------------------------------------------------------------------------------#
# Extract posterior samples and set parameters for summary tables, figures
#------------------------------------------------------------------------------#

# Extract posterior samples from jagsUI object and combine into one dataframe
samples <- out$samples
samples <- do.call(rbind, samples)

# Set lower and upper quantiles for credible intervals
lcl <- 0.025
ucl <- 0.975

# Set number of digits for reporting summary statistics
digits <- 2

# Identify a subset of posterior samples we want to use for figures (spatial
# predictions and trends). Keeping 1000s of samples will create memory issues.
nsamp <- 500
subsamples <- floor(seq(1, nrow(samples), length = nsamp))

#------------------------------------------------------------------------------#
# Summaries of posterior distributions for covariate effects (with names)
#------------------------------------------------------------------------------#

# Extract the number of covariates on each parameter
n_cov_psi <- ifelse(exists("cov_psi"), ncol(cov_psi), 0)
n_cov_p <- ifelse(exists("cov_p"), ncol(cov_p), 0)
n_cov_eps <- ifelse(exists("cov_eps"), ncol(cov_eps), 0)
n_cov_gam <- ifelse(exists("cov_gam"), ncol(cov_gam), 0)

# Extract samples for covariate effects
exclude <- paste(c("mean", "PAO", "0", "deviance"), collapse = "|")
covar_samples <- samples[,str_detect(colnames(samples), exclude, negate = TRUE)]
covar_samples <- covar_samples[,order(colnames(covar_samples))]

# Function to calculate the proportion of posterior samples that are > 0 when 
# the median is > 0 or < 0 for medians < 0 (f in summary table)
prop_samples <- function(x) {
  if (median(x) > 0) {
    prop <- sum(x > 0) / length(x)
  } else {
    prop <- sum(x < 0) / length(x)
  }
  return(prop)
} 

# Create summary table
jags_name <- colnames(covar_samples)
parameter <- c(rep("extinction", n_cov_eps), rep("colonization", n_cov_gam), 
               rep("detection", n_cov_p), rep("init occ", n_cov_psi))
covariate <- c(colnames(cov_eps), colnames(cov_gam), 
               colnames(cov_p), colnames(cov_psi))
covariate <- str_remove_all(covariate, "_z|_z2")
covar_summary <- data.frame(jags_name, parameter, covariate)
covar_summary <- covar_summary %>%
  mutate(mean = round(apply(covar_samples, 2, mean), digits),
         sd = round(apply(covar_samples, 2, sd), digits),
         lcl = round(apply(covar_samples, 2, quantile, lcl), digits), 
         median = round(apply(covar_samples, 2, quantile, 0.5), digits),
         ucl = round(apply(covar_samples, 2, quantile, ucl), digits),
         overlap0 = ifelse(lcl * ucl > 0, FALSE, TRUE),
         f = round(apply(covar_samples, 2, prop_samples), digits))
covar_summary    

# Remove unnecessary objects
rm(covar_samples)

#------------------------------------------------------------------------------#
# Summaries of posterior distributions for proportion of sampled sites occupied
#------------------------------------------------------------------------------#

# Extract samples for PAO
pao_samples <- samples[,str_detect(colnames(samples), "PAO")]

# Create summary table
jags_name <- colnames(pao_samples)
year <- min(sitetrans$start_yr):max(sitetrans$end_yr)
pao_summary <- data.frame(jags_name, year)
pao_summary <- pao_summary %>%
  mutate(mean = round(apply(pao_samples, 2, mean), digits),
         sd = round(apply(pao_samples, 2, sd), digits),
         lcl = round(apply(pao_samples, 2, quantile, lcl), digits), 
         median = round(apply(pao_samples, 2, quantile, 0.5), digits),
         ucl = round(apply(pao_samples, 2, quantile, ucl), digits))
pao_summary         

# Remove unnecessary objects
rm(pao_samples)

#------------------------------------------------------------------------------#
# Plotting predicted probability of occupancy in first year
#------------------------------------------------------------------------------#

# Identify covariates in model
psi_covars <- unique(str_remove(colnames(cov_psi), "_z2|_z|[:digit:]"))
# Unzip park rasters
park_folder <- paste0("data/covariates/rasters-", PARK, "/")
park_zip <- paste0("data/covariates/rasters-", PARK, ".zip")
if (length(list.files(park_folder)) == 0) {
  unzip(park_zip)
}

# Select just those rasters that were used in model of initial occupancy
park_rasters <- list.files(path = park_folder, 
                           pattern = ".tif", 
                           full.names = TRUE)
if ("elev" %in% psi_covars) {
  raster_covars <- replace(x = psi_covars,
                           list = which(psi_covars == "elev"), 
                           values = "DEM")
}

# Load rasters into a list
raster_list <- list()
for (i in 1:length(raster_covars)) {
  park_raster_ind <- which(str_detect(park_rasters, raster_covars[i]))
  raster_list[[i]] <- terra::rast(park_rasters[park_raster_ind])
}
names(raster_list) <- psi_covars

# Load park boundary to crop rasters
boundaries <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
boundary <- subset(boundaries, boundaries$UNIT_CODE == PARK)
raster_list <- lapply(raster_list, FUN = crop, ext(boundary))

# Identify categorical covariates (only vegclasses for now)
cat_covar <- "vegclass"

# Extract rasters with continuous covariates
psi_covarsc <- psi_covars[!psi_covars %in% cat_covar]
raster_cont <- raster_list[psi_covarsc]

# Standardize continuous rasters and create quadratic where needed
for (cov in names(raster_cont)) {
  cov_mn <- mean(spatial_covs[,cov])
  cov_sd <- sd(spatial_covs[,cov])
  raster_cont[[cov]] <- (raster_cont[[cov]] - cov_mn)/cov_sd
  if (paste0(cov, "2") %in% covar_summary$covariate) {
    raster_cont <- c(raster_cont, raster_cont[[cov]] ^ 2)
    names(raster_cont)[[length(raster_cont)]] <- paste0(cov, "2")
  }
}
raster_cont <- raster_cont[order(names(raster_cont))]
rast_final <- rast(raster_cont)

# Extract rasters with categorical covariates
if (cat_covar %in% psi_covars) {
  raster_cat <- raster_list[[cat_covar]]
  raster_cat[raster_cat == 4] <- NA
  veg2 <- 1 * (raster_cat == 2)
  names(veg2) <- "veg2"
  veg3 <- 1 * (raster_cat == 3)
  names(veg3) <- "veg3"
  rast_final <- c(rast_final, veg2, veg3)
}

# Create a raster layer with just ones for the intercept
rast_final <- c(rast(rast_final[[1]], vals = 1), rast_final)
names(rast_final)[[1]] <- "int"

# Extract posterior samples for initial occupancy parameters (n = nsamp)
psi_samp <- samples[subsamples, grep("beta_psi", colnames(samples))]

# Convert SpatRaster to a dataframe (one row for each cell, each column = layer)
# note: this removes cells/rows with NAs
rast_final_df <- as.data.frame(rast_final, cell = TRUE)
# Do the math (results in predictions on logit scale for each cell [row] and 
# posterior sample [column])
preds1_df <- as.matrix(rast_final_df[,-1]) %*% t(psi_samp) 
# Convert to probability scale (can take a few seconds)
preds1_df <- exp(preds1_df)/(1 + exp(preds1_df))
# Re-attach cell numbers
preds1_df <- cbind(cell = rast_final_df$cell, preds1_df)

# Calculate the median prediction in each cell and create a raster
preds1_median_df <- cbind(preds1_df[,1], apply(preds1_df[,-1], 1, median))
preds1_median <- rast(rast_final[[1]])
preds1_median[preds1_median_df[,1]] <- preds1_median_df[,2]
# plot(preds1_median)

# Calculate the SD of predictions in each cell and create a raster
preds1_sd_df <- cbind(preds1_df[,1], apply(preds1_df[,-1], 1, sd))
preds1_sd <- rast(rast_final[[1]])
preds1_sd[preds1_sd_df[,1]] <- preds1_sd_df[,2] 
# plot(preds1_sd)

# Remove unnecessary objects
rm(preds1_median_df, preds1_sd_df, raster_list, raster_cont, raster_cat, 
   psi_samp, rast_final_df)

# List of TODOs: 
  # add options to save these rasters and/or plots?
  # think about whether we want to make predictions to highest elevations 
    # where we didn't survey (and where SD of predictions if very high)
  # improve object names throughout

#------------------------------------------------------------------------------#
# Estimate the probability of occupancy across the park in subsequent years
#------------------------------------------------------------------------------#

# List of potential seasonal covariates
seas_covariates <- c("monsoon_ppt") # Will add more in time

# Identify all the rasters we'll need for seasonal covariates in models of 
# extinction and colonization
seas_gam <- str_subset(COVARS_GAM, paste(seas_covariates, collapse = "|"))
seas_eps <- str_subset(COVARS_EPS, paste(seas_covariates, collapse = "|"))
seas_both <- unique(c(seas_gam, seas_eps))
 
# Need to unzip and load rasters with weather data (or eventually any rasters 
# with annually varying covariates)
weather_folder <- "data/covariates/weather-derived-rasters/"
weather_zip <- "data/covariates/weather-derived.zip"
# Unzip weather folder first, if necessary
if (length(list.files(weather_folder)) == 0) {
  unzip(weather_zip, overwrite = TRUE)
}
# List files in weather folder
weather_files <- list.files(weather_folder, full.names = TRUE)

# For each seasonal covariate in the model, create a list of rasters
for (cov in seas_both) {
  cov_files <- weather_files[str_detect(weather_files, cov)]
  # Remove rasters associated with periods outside the years of interest
  # (eg, monsoon rainfall in year x explains transitions between yrs x and x+1)
  cov_yrs <- year[-length(year)]
  cov_yrs_search <- paste0(as.character(cov_yrs), collapse = "|")
  cov_files <- cov_files[str_detect(cov_files, cov_yrs_search)]

  # Identify the mean and SD from surveyed locations for standardizing
  cov_mn <- mean(sitetrans[,cov])
  cov_sd <- sd(sitetrans[,cov])
  
  # Load each raster into a list, resample (so it has the same geometry as 
  # spatial rasters in rast_final created above), and standardize
  cov_list <- list()
  for (i in 1:length(cov_files)) {
    cov_list[[i]] <- rast(cov_files[i])
    cov_list[[i]] <- resample(cov_list[[i]], rast_final, method = "near")
    cov_list[[i]] <- mask(cov_list[[i]], boundary)
    cov_list[[i]] <- (cov_list[[i]] - cov_mn)/cov_sd
  }  
  names(cov_list) <- paste(cov, cov_yrs, sep = "_")
  
  # Name the raster list appropriately (e.g., monsoon_ppt_list)
  assign(paste0(cov, "_list"), cov_list)
}  

# Calculate colonization probability for each cell, year
  # (For now, can accommodate models that have interactions but not quadratics)
  
  # Extract posterior samples for extinction parameters (n = nsamp)
  gam_samp <- samples[subsamples, grep("beta_gam", colnames(samples))]
  
  # Identify spatial covariates in model (if any)
  spat_gam <- str_subset(COVARS_GAM, 
                         paste(seas_covariates, collapse = "|"),
                         negate = TRUE)
  
  # Load rasters associated with spatial covariates into a list 
  # (also, crop and mask each raster to park boundary and standardize)
  if ("elev" %in% spat_gam) {
    spat_gam_r <- replace(x = spat_gam,
                          list = which(spat_gam == "elev"), 
                          values = "DEM")
  }  
  spat_gam_list <- list()
  for (i in 1:length(spat_gam_r)) {
    raster_ind <- which(str_detect(park_rasters, spat_gam_r[i]))
    spat_gam_list[[i]] <- rast(park_rasters[raster_ind])
    spat_gam_list[[i]] <- crop(spat_gam_list[[i]], ext(boundary))
    spat_gam_list[[i]] <- mask(spat_gam_list[[i]], boundary)
    cov_mn <- mean(spatial_covs[,spat_gam[i]])
    cov_sd <- sd(spatial_covs[,spat_gam[i]])
    spat_gam_list[[i]] <- (spat_gam_list[[i]] - cov_mn)/cov_sd
  }
  names(spat_gam_list) <- spat_gam
  
  # Create SpatRaster with a layer for the intercept, a layer for each 
  # spatial covariate, and layers for each year of seasonal covariate data
  gam_list <- list()
  for (i in 1:length(seas_gam)) {
    gam_list <- c(gam_list, get(paste0(seas_gam[[i]], "_list")))
  }
  gam_list <- c(spat_gam_list, gam_list)
  gam_raster <- rast(gam_list)
  gam_raster <- c(rast(gam_raster[[1]], vals = 1), gam_raster)
  names(gam_raster)[[1]] <- "int"
  
  # Convert to a dataframe
  gam_raster_df <- as.data.frame(gam_raster, cell = TRUE) 
  
  # Create season-specific dataframes with a column of 1s for the intercept and
  # one column for each spatial and seasonal covariate
  for (i in 1:length(cov_yrs)) {
    # Extract the columns we need from gam_raster_df
    year_cols <- str_subset(names(gam_raster_df), as.character(cov_yrs[i]))
    gam_yr <- gam_raster_df %>%
      select(int, all_of(spat_gam), all_of(year_cols))
    # Strip the year out of any column name
    colnames(gam_yr) <- str_remove(colnames(gam_yr), paste0("_", cov_yrs[i]))
    
    # Create columns with interactions (if needed)
    n_ints <- str_subset(ls(), "GAM_INT")
    if (length(n_ints) > 0) {
      for (j in 1:length(n_ints)) {
        int <- get(n_ints[j])
        int_ind <- which(colnames(gam_yr) %in% int)
        gam_yr[,paste(int, collapse = "_")] <- gam_yr[,int_ind[1]] * gam_yr[,int_ind[2]]
      }
    }
    
    # Make sure columns are in the same order they appear in the model
    col_order <- covar_summary$covariate[covar_summary$parameter == "colonization"]
    gam_yr <- gam_yr[,c("int", col_order)]
    
    # Do the math (results in predictions of colonization probability on the 
    # logit scale for each cell [row] and posterior sample [column])
    gam_prob <- as.matrix(gam_yr) %*% t(gam_samp)
    # Convert to probability scale
    gam_prob <- exp(gam_prob)/(1 + exp(gam_prob))
    # Attach cell numbers
    gam_prob <- cbind(cell = gam_raster_df$cell, gam_prob)
    # Give this matrix a name that identifies the transition year (eg, gam_1)
    assign(paste0("gam_", i), gam_prob)
  }
  # Remove unnecessary objects
  rm(gam_prob, gam_yr, gam_raster_df, spat_gam_list, gam_list)

# Calculate extinction probability for each cell, year
  # (For now, can accommodate models that have interactions but not quadratics)

  # Extract posterior samples for extinction parameters (n = nsamp)
  eps_samp <- samples[subsamples, grep("beta_eps", colnames(samples))]

  # Identify spatial covariates in model (if any)
  spat_eps <- str_subset(COVARS_EPS, 
                         paste(seas_covariates, collapse = "|"),
                         negate = TRUE)
  
  # Load rasters associated with spatial covariates into a list 
  # (also, crop and mask each raster to park boundary and standardize)
  if ("elev" %in% spat_eps) {
    spat_eps_r <- replace(x = spat_eps,
                          list = which(spat_eps == "elev"), 
                          values = "DEM")
  }  
  spat_eps_list <- list()
  for (i in 1:length(spat_eps_r)) {
    raster_ind <- which(str_detect(park_rasters, spat_eps_r[i]))
    spat_eps_list[[i]] <- rast(park_rasters[raster_ind])
    spat_eps_list[[i]] <- crop(spat_eps_list[[i]], ext(boundary))
    spat_eps_list[[i]] <- mask(spat_eps_list[[i]], boundary)
    cov_mn <- mean(spatial_covs[,spat_eps[i]])
    cov_sd <- sd(spatial_covs[,spat_eps[i]])
    spat_eps_list[[i]] <- (spat_eps_list[[i]] - cov_mn)/cov_sd
  }
  names(spat_eps_list) <- spat_eps
  
  # Create SpatRaster with a layer for the intercept, a layer for each 
  # spatial covariate, and layers for each year of seasonal covariate data
  eps_list <- list()
  for (i in 1:length(seas_eps)) {
    eps_list <- c(eps_list, get(paste0(seas_eps[[i]], "_list")))
  }
  eps_list <- c(spat_eps_list, eps_list)
  eps_raster <- rast(eps_list)
  eps_raster <- c(rast(eps_raster[[1]], vals = 1), eps_raster)
  names(eps_raster)[[1]] <- "int"

  # Convert to a dataframe
  eps_raster_df <- as.data.frame(eps_raster, cell = TRUE) 

  # Create season-specific dataframes with a column of 1s for the intercept and
  # one column for each spatial and seasonal covariate
  for (i in 1:length(cov_yrs)) {
    # Extract the columns we need from eps_raster_df
    year_cols <- str_subset(names(eps_raster_df), as.character(cov_yrs[i]))
    eps_yr <- eps_raster_df %>%
      select(int, all_of(spat_eps), all_of(year_cols))
    # Strip the year out of any column name
    colnames(eps_yr) <- str_remove(colnames(eps_yr), paste0("_", cov_yrs[i]))
    
    # Create columns with interactions (if needed)
    n_ints <- str_subset(ls(), "EPS_INT")
    if (length(n_ints) > 0) {
      for (j in 1:length(n_ints)) {
        int <- get(n_ints[j])
        int_ind <- which(colnames(eps_yr) %in% int)
        eps_yr[,paste(int, collapse = "_")] <- eps_yr[,int_ind[1]] * eps_yr[,int_ind[2]]
      }
    }
    
    # Make sure columns are in the same order they appear in the model
    col_order <- covar_summary$covariate[covar_summary$parameter == "extinction"]
    eps_yr <- eps_yr[,c("int", col_order)]
    
    # Do the math (results in predictions of extinction probability on the 
    # logit scale for each cell [row] and posterior sample [column])
    eps_prob <- as.matrix(eps_yr) %*% t(eps_samp)
    # Convert to probability scale
    eps_prob <- exp(eps_prob)/(1 + exp(eps_prob))
    # Attach cell numbers
    eps_prob <- cbind(cell = eps_raster_df$cell, eps_prob)
    # Give this matrix a name that identifies the transition year (eg, eps_1)
    assign(paste0("eps_", i), eps_prob)
  }
  
  # Remove unnecessary objects
  rm(eps_prob, eps_yr, eps_raster_df, spat_eps_list, eps_list)
  rm(list = str_subset(ls(), "_list"))
  # BUT be sure to keep eps_raster

# Need to make sure we're making predictions of all parameters to the same set 
# of cells (could differ because of covariate availability)
  # Any cells in the gam/eps set that aren't in the preds1_df (init occ prob)?
  eps_noocc <- eps_1[,"cell"][!eps_1[,"cell"] %in% preds1_df[,"cell"]]
  length(eps_noocc) 
  # Any cells in preds_df that aren't in gam/eps sets?
  occ_noeps <- preds1_df[,"cell"][!preds1_df[,"cell"] %in% eps_1[,"cell"]]
  length(occ_noeps)
  # Going forward, using set of cell numbers present in all predictions:
  cells <- sort(intersect(preds1_df[,"cell"], eps_1[,"cell"]))
  preds1_df <- preds1_df[preds1_df[,"cell"] %in% cells,]  
  
# Finally, run through each season and draw values of latent occupancy state (z) 
# for each cell and year

  # Draw values of latent occupancy state in year 1: z[,1]
  # (Occupancy probabilities are in preds1_df)
  z_1 <- matrix(rbinom(length(preds1_df[,-1]), 1, preds1_df[,-1]),
                nrow = nrow(preds1_df[,-1]), ncol = ncol(preds1_df[,-1]))
  # Calculate PAO in year 1
  PAO_1 <- apply(z_1, 2, function(x) sum(x)/length(x))
    # PAO_1 is a vector that has the estimated proportion of cells occupied
    # for each posterior sample (n = nsamp)
  
  # For each subsequent year...
  for (t in 1:length(cov_yrs)) {
    
    # Get the matrix of colonization probabilities for transition from year t to t+1
    gam <- get(paste0("gam_", t))
    gam <- gam[gam[,"cell"] %in% cells,]
    gam <- gam[, -1]
    
    # Get the matrix of extinction probabilities for transition from year t to t+1
    eps <- get(paste0("eps_", t))
    eps <- eps[eps[,"cell"] %in% cells,]
    eps <- eps[, -1]
    
    # Get the matrix of latent states in year t
    if (t == 1) {
      z <- z_1
    } else {
      z <- z_new
    }
    
    # Calculate the probability of occupancy for each cell in year t + 1  
    Ez <- gam * (1 - z) + (1 - eps) * z
    
    # Draw the latent occupancy state from a Bernoulli for each cell in season t+1
    z_new <- matrix(rbinom(length(Ez), 1, Ez), nrow = nrow(Ez), ncol = ncol(Ez))
    
    # Save matrix with occupancy probabilities (Ez)
    assign(paste0("Ez_", t + 1), Ez)
    
    # Calculate estimates of PAO (across entire park) for each year 
    # (will need this to estimate park-wide trends in occupancy)
    assign(paste0("PAO_", t + 1), 
           apply(z_new, 2, function(x) sum(x)/length(x)))
  } 

# Remove unnecessary objects
rm(Ez, z_new, z, gam, eps,
   gam_1, gam_2, gam_3, gam_4, gam_5,
   eps_1, eps_2, eps_3, eps_4, eps_5)  

# Create rasters with summaries of occupancy probabilities in each year
for (t in 1:length(cov_yrs)) {
  Ez <- get(paste0("Ez_", t + 1))
  Ez_stats <- data.frame(cell = cells, 
                         median = apply(Ez, 1, median),
                         sd = apply(Ez, 1, sd))
  preds_median <- rast(preds1_median)
  preds_median[Ez_stats[,1]] <- Ez_stats[,2]
  names(preds_median) <- paste0("median_", cov_yrs[t] + 1)
  assign(paste0("Ez_median_", t + 1), preds_median)
  preds_sd <- rast(preds1_sd)
  preds_sd[Ez_stats[,1]] <- Ez_stats[,3]
  names(preds_sd) <- paste0("sd_", cov_yrs[t] + 1)
  assign(paste0("Ez_sd_", t + 1), preds_sd)  
}

# Merge all rasters with medians of occupancy probability estimates into one,
# multi-layer SpatRaster (Ez_medians) and all rasters with SDs of occupancy 
# probability estimates into one multi-layer SpatRaster (Ez_sds)
Ez_medians <- preds1_median  
names(Ez_medians) <- paste0("median_", year[1])
Ez_sds <- preds1_sd
names(Ez_sds) <- paste0("sd_", year[1])
for (t in 2:(length(cov_yrs)+1)) {
  new_median <- get(paste0("Ez_median_", t))
  Ez_medians <- c(Ez_medians, new_median)
  new_sd<- get(paste0("Ez_sd_", t))
  Ez_sds <- c(Ez_sds, new_sd)
}

# Look at median occupancy probabilities by year
plot(Ez_medians)
# Look at SDs of those estimates by year
plot(Ez_sds)

# Remove unnecessary objects
rm(list = str_subset(ls(), "Ez_median_"))
rm(list = str_subset(ls(), "Ez_sd_"))

#------------------------------------------------------------------------------#
# Estimate trend in occupancy for surveyed locations
#------------------------------------------------------------------------------#

# Extract subset of posterior samples (n = nsamp)
subsamples <- floor(seq(1, nrow(samples), length = nsamp))
pao <- samples[subsamples, grep("PAO", colnames(samples))]

# For each MCMC iteration, estimate a linear trend in logit(occupancy)
# (note: this is a trend for surveyed locations only)
pao_logit <- log(pao / (1 - pao))
year_trend <- 0:(ncol(pao_logit) - 1)
trends <- data.frame(iter = 1:nrow(pao_logit), int = NA, slope = NA)
for (i in 1:nrow(pao_logit)) {
  m <- lm(pao_logit[i,] ~ year_trend)
  trends[i,2:3] <- coef(m)
}

# Summarize trend estimates
summary(trends$slope)
hist(trends$slope, breaks = 50)

# Plot trends on the logit scale (each gray line represents one MCMC iteration)
trends <- trends %>%
  mutate(last_yr = int + max(year_trend) * slope)

ggplot() +
  geom_segment(trends,
               mapping = aes(x = min(year), xend = max(year), 
                             y = int, yend = last_yr),
               size = 0.3, col = "gray") +
  geom_segment(trends,
               mapping = aes(x = min(year), xend = max(year), 
                             y = median(int), yend = median(last_yr)),
               size = 0.8, col = "dodgerblue3") +
  labs(x = "Year", y = "logit(Proportion of sampled sites occupied)")

# Prep data to plot (logit-linear) trends on the probability scale.
# Create a sequence of values that spans yr_trend
yr_predict <- seq(0, max(year_trend), length = 100)
# Create a matrix, with a column of 1s (for the intercept) and yr_predict
yr_predict_matrix <- rbind(1, yr_predict)
# Predict values: logit(occupancy) = beta0 + slope*yr_predict (matrix math)
preds_logit <- as.matrix(trends[,c("int", "slope")]) %*% yr_predict_matrix
# Convert predictions to probability scale. 
preds_prob <- exp(preds_logit)/(1 + exp(preds_logit))
colnames(preds_prob) <- yr_predict
# Calculate the median value across iterations for each value of yr_predict
preds_median <- data.frame(Year = yr_predict + min(year), 
                           Occupancy = apply(preds_prob, 2, median))
# Add column to identify MCMC iteration
preds_prob <- cbind(preds_prob, iter = 1:nrow(preds_prob))
# Convert predictions to long form for ggplot
preds_long <- preds_prob %>%
  as.data.frame %>%
  pivot_longer(cols = !iter,
               names_to = "Year",
               values_to = "Occupancy") %>%
  mutate(Year = as.numeric(Year) + min(year)) %>%
  as.data.frame

# Plot trends on the probability scale (each line represents one MCMC iteration)
ggplot() +
  geom_line(preds_long,
            mapping = aes(x = Year, y = Occupancy, group = iter),
            col = "gray") + 
  geom_line(preds_median,
            mapping = aes(x = Year, y = Occupancy),
            size = 0.8, col = "dodgerblue3") +
  labs(y = "Proportion of sampled sites occupied")

#------------------------------------------------------------------------------#
# Estimate trend in occupancy across the entire park
#------------------------------------------------------------------------------#

# Using the PAO_1 ... PAO_X objects created above
paos <- paste0("PAO_", 1:length(year))
pao_park <- do.call(cbind, mget(paos))

# For each MCMC iteration, estimate a linear trend in logit(occupancy)
pao_park_logit <- log(pao_park / (1 - pao_park))
year_trend <- 0:(ncol(pao_park_logit) - 1)
trends_park <- data.frame(iter = 1:nrow(pao_park_logit), int = NA, slope = NA)
for (i in 1:nrow(pao_park_logit)) {
  m <- lm(pao_park_logit[i,] ~ year_trend)
  trends_park[i,2:3] <- coef(m)
}

# Summarize trend estimates
summary(trends_park$slope)
hist(trends_park$slope, breaks = 50)

# Plot trends on the logit scale (each gray line represents one MCMC iteration)
trends_park <- trends_park %>%
  mutate(last_yr = int + max(year_trend) * slope)

ggplot() +
  geom_segment(trends_park,
               mapping = aes(x = min(year), xend = max(year), 
                             y = int, yend = last_yr),
               size = 0.3, col = "gray") +
  geom_segment(trends_park,
               mapping = aes(x = min(year), xend = max(year), 
                             y = median(int), yend = median(last_yr)),
               size = 0.8, col = "dodgerblue3") +
  labs(x = "Year", y = "logit(Proportion of park occupied)")

# Prep data to plot (logit-linear) trends on the probability scale.
# Create a sequence of values that spans yr_trend
yr_predict <- seq(0, max(year_trend), length = 100)
# Create a matrix, with a column of 1s (for the intercept) and yr_predict
yr_predict_matrix <- rbind(1, yr_predict)
# Predict values: logit(occupancy) = beta0 + slope*yr_predict (matrix math)
preds_park_logit <- as.matrix(trends_park[,c("int", "slope")]) %*% yr_predict_matrix
# Convert predictions to probability scale. 
preds_park_prob <- exp(preds_park_logit)/(1 + exp(preds_park_logit))
colnames(preds_park_prob) <- yr_predict
# Calculate the median value across iterations for each value of yr_predict
preds_park_median <- data.frame(Year = yr_predict + min(year), 
                                Occupancy = apply(preds_park_prob, 2, median))
# Add column to identify MCMC iteration
preds_park_prob <- cbind(preds_park_prob, iter = 1:nrow(preds_park_prob))
# Convert predictions to long form for ggplot
preds_park_long <- preds_park_prob %>%
  as.data.frame %>%
  pivot_longer(cols = !iter,
               names_to = "Year",
               values_to = "Occupancy") %>%
  mutate(Year = as.numeric(Year) + min(year)) %>%
  as.data.frame

# Plot trends on the probability scale (each line represents one MCMC iteration)
ggplot() +
  geom_line(preds_park_long,
            mapping = aes(x = Year, y = Occupancy, group = iter),
            col = "gray") + 
  geom_line(preds_park_median,
            mapping = aes(x = Year, y = Occupancy),
            size = 0.8, col = "dodgerblue3") +
  labs(y = "Proportion of park occupied")

#------------------------------------------------------------------------------#
# Calculate marginal effects for a continuous covariate (example)
#------------------------------------------------------------------------------#

# Estimate how colonization probability would change with a 1-unit (SD) 
# increase in elev

  # Notes: logit is the log of the odds
  # Odds are the probability event happens / probability event doesn't happen
  # So logit(gamma[i]) = log(gamma[i] / (1-gamma[i])) = log(odds)
  # Odds a site gets colonized at mean elevation = exp(beta_gam0)
  # Odds a site gets colonized with 1-SD increase in elevation = 
  # exp(beta_gam0 + beta_gam1 * 1) = exp(beta_gam0)*exp(beta_gam1)
  # So the odds will change by a FACTOR of exp(beta_gam1):
  # Odds[elev+1SD] = Odds[mean elev] * exp(beta_gam1)

# For sagw-leca-MS-2022-10-19.Rdata model, elevation effect is beta_gam[1]
# logit(gamma[i]) = beta_gam0 + beta_gam[1] * elevation[i]
beta_gam <- samples[,"beta_gam[1]"]
change <- exp(beta_gam)
mean(change) # 0.49
quantile(change, probs = c(lcl, ucl)) #0.27, 0.78
# The odds a site is colonized are estimated to be 51% lower for each 
# 1-SD increase in elevation (95% CI = 22-73%) [assuming other covariates held
# constant]

#------------------------------------------------------------------------------#
# Plot marginal effects for a continuous covariate (example)
#------------------------------------------------------------------------------#

# Estimate how detection probability changes with day-of-year

# For sagw-leca-MS-2022-10-19.Rdata model, linear and quadratic effect of day
# are beta_p[2:3]
# logit(p[i]) = beta_p0 + beta_p[2] * doy[i] + beta_p[3] * doy[i] * doy[i]

# Generate a vector of days that span the range that occurred during study
day_nums <- seq(min(surveys$day), max(surveys$day), length = 100)
# Standardize these values (Make sure to grab the mean/SD from the dataframe 
# where values were standardized in MSoccupancy-generic.R. For day number, this 
# was in the occasions dataframe)
day_nums_z <- (day_nums - mean(occasions$mid_yday)) / sd(occasions$mid_yday)
# Create a matrix of covariate values (including the intercept [1])
X_p <- cbind(int = 1, doy = day_nums_z, doy2 = day_nums_z *day_nums_z)
# Create a matrix with posterior samples for the parameters we need
betas_p <- samples[,c("beta_p0", "beta_p[2]", "beta_p[3]")]
# Generate a matrix of predicted values on the logit scale
# Matrix dimensions = 100 x 3000 (3000 predicted values for each day in seq)
pred_logit_p <- X_p %*% t(betas_p)
# Backtransform to the probability scale
pred_prob_p <- exp(pred_logit_p) / (1 + exp(pred_logit_p))
# Calculate mean and credible interval for each value
mean_p <- apply(pred_prob_p, 1, mean)
cri_p <- apply(pred_prob_p, 1, quantile, probs = c(lcl, ucl)) 
# Put objects in dataframe for ggplot
plot_data <- data.frame(mean = mean_p,
                        day_nums = day_nums,
                        lcl = cri_p[1,],
                        ucl = cri_p[2,])

# Plot predictions
  # For day-of-year, transform day numbers we want on x-axis to date 
  # (for easier interpretation)
  day_nums_axis <- seq(min(day_nums), max(day_nums), by = 10)
  days_axis <- parse_date_time(paste(2022, as.character(day_nums_axis)), 
                          orders = "yj")

par(mfrow = c(1, 1))
ggplot() + 
  geom_ribbon(plot_data,
              mapping = aes(x = day_nums, ymin = lcl, ymax = ucl), 
              fill = "gray70") +
  geom_line(plot_data,
            mapping = aes(x = day_nums, y = mean)) +
  scale_y_continuous(limits = c(0.3, 0.9)) +
  scale_x_continuous(breaks = day_nums_axis, labels = format(days_axis, "%m-%d")) + 
  labs(x = "Date", y = "Detection probability") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

#------------------------------------------------------------------------------#
# Calculate marginal effects for a categorical covariate (example)
#------------------------------------------------------------------------------#

# Estimate how vegetation class affects probabilities of occupancy in the first 
# year 

# For sagw-leca-MS-2022-10-13.Rdata model, vegclass1 (low gradient desert) is 
# the reference level, effect of vegclass2 (low hillslope, north-facing) is 
# beta_psi[5], and effect of vegclass3 (medium-high gradient) is beta_psi[6]

# Probability of occupancy in vegclass1 (at mean values of other covariates)
beta_psi0 <- samples[,"beta_psi0"]
vegclass1 <- exp(beta_psi0)/(1 + exp(beta_psi0)) 
mean(vegclass1) # 0.53
quantile(vegclass1, probs = c(lcl, ucl)) #0.25, 0.82
# Probability of initial occupancy = 0.53 (95% CI = 0.25, 0.82)

# Probability of occupancy in vegclass2 (at mean values of other covariates)
beta_veg2 <- samples[,"beta_psi[5]"]
vegclass2L <- beta_psi0 + beta_veg2
vegclass2 <- exp(vegclass2L)/(1 + exp(vegclass2L)) 
mean(vegclass2) # 0.05
quantile(vegclass2, probs = c(lcl, ucl)) #0.00, 0.25
# Probability of initial occupancy = 0.05 (95% CI = 0.00, 0.25)

