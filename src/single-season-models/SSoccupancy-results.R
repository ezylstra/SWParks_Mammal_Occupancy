################################################################################
# Process results from a single-season occupancy analysis

# ER Zylstra
# Updated 2022-12-13
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(jagsUI)
library(MCMCvis)
library(ggplot2)

rm(list = ls())

#------------------------------------------------------------------------------#
# Load the workspace that contains the JAGS model
#------------------------------------------------------------------------------#

PARK <- "SAGW"
SPECIES <- "LECA"
YEAR <- 2020
DATE <- "2022-12-13"
DATE <- str_remove_all(DATE, "-")
output_file <- paste0("output/models/",
                      tolower(PARK), "-",
                      tolower(SPECIES), "-", 
                      YEAR, "-",
                      DATE,".Rdata")
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

# Load park boundary to crop rasters
boundaries <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
boundary <- subset(boundaries, boundaries$UNIT_CODE == PARK)

#------------------------------------------------------------------------------#
# Summaries of posterior distributions for covariate effects (with names)
#------------------------------------------------------------------------------#

# Extract the number of covariates on each parameter
n_cov_psi <- ifelse(exists("cov_psi"), ncol(cov_psi), 0)
n_cov_p <- ifelse(exists("cov_p"), ncol(cov_p), 0)

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
parameter <- c(rep("detection", n_cov_p), rep("init occ", n_cov_psi))
covariate <- c(colnames(cov_p), colnames(cov_psi))
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
# Calculate predicted probability of occupancy, across park
#------------------------------------------------------------------------------#

# Unzip park rasters
park_folder <- paste0("data/covariates/rasters-", PARK, "/")
if (PARK == "ORPI") {
  park_zip1 <- paste0("data/covariates/rasters-ORPI-dist.zip") 
  park_zip2 <- paste0("data/covariates/rasters-ORPI-topo.zip") 
  if (length(list.files(park_folder)) == 0) {
    unzip(park_zip1, overwrite = TRUE)
    unzip(park_zip2, overwrite = TRUE)
  }  
} else {
  park_zip <- paste0("data/covariates/rasters-", PARK, ".zip")
  if (length(list.files(park_folder)) == 0) {
    unzip(park_zip)
  }
}

# Generate list of available rasters
park_rasters <- list.files(path = park_folder, 
                           pattern = ".tif", 
                           full.names = TRUE)

# Load a single raster to obtain geometry, and set all values within the park 
# to 1 (for an intercept layer)
rast_int <- terra::rast(park_rasters[1])
rast_int <- terra::rast(rast_int, vals = 1)
rast_int <- terra::crop(rast_int, boundary)
rast_int <- terra::mask(rast_int, boundary)
names(rast_int) <- "int"
psi_rasters <- rast_int

# If there are covariates in the model, load and process appropriate rasters
if (exists("cov_psi")) {
  # Identify covariates in model
  psi_covars <- unique(str_remove(colnames(cov_psi), "_z2|_z|[:digit:]"))
  if ("elev" %in% psi_covars) {
    raster_covars <- replace(x = psi_covars,
                             list = which(psi_covars == "elev"), 
                             values = "DEM")
  }
  
  # Load rasters into a list
  psi_list <- list()
  for (i in 1:length(raster_covars)) {
    park_raster_ind <- which(str_detect(park_rasters, raster_covars[i]))
    psi_list[[i]] <- terra::rast(park_rasters[park_raster_ind])
  }
  names(psi_list) <- psi_covars
  
  # Crop rasters
  psi_list <- lapply(psi_list, FUN = crop, ext(boundary))
  
  # If vegclass (which is categorical) is in the model, create a layer for each
  # dummy variables in the model (SAGW: vegclasses 2 and 3)
  if ("vegclass" %in% psi_covars) {
    veg_rast <- psi_list[["vegclass"]]
    veg_rast[veg_rast == 4] <- NA
    vegclass2 <- 1 * (veg_rast == 2)
    names(vegclass2) <- "vegclass2"
    vegclass3 <- 1 * (veg_rast == 3)
    names(vegclass3) <- "vegclass3"
    psi_rasters <- c(psi_rasters, vegclass2, vegclass3)
  }
  
  # If there are continuous covariates, append associated layers to psi_rasters
  psi_covarsc <- psi_covars[!psi_covars %in% c("vegclass")]
  if (length(psi_covarsc) > 0) {
    psi_cont <- psi_list[psi_covarsc]
    
    # Standardize continuous rasters and create quadratic where needed
    for (cov in names(psi_cont)) {
      cov_mn <- mean(spatial_covs[,cov])
      cov_sd <- sd(spatial_covs[,cov])
      psi_cont[[cov]] <- (psi_cont[[cov]] - cov_mn)/cov_sd
      if (paste0(cov, "2") %in% covar_summary$covariate) {
        psi_cont <- c(psi_cont, psi_cont[[cov]] ^ 2)
        names(psi_cont)[[length(psi_cont)]] <- paste0(cov, "2")
      }
    }

    psi_rasters <- c(psi_rasters, rast(psi_cont))
    
  }
  
  # Put raster layers in order they appear in model
  cov_order <- covar_summary$covariate[covar_summary$parameter == "init occ"]
  psi_rasters <- psi_rasters[[c("int", cov_order)]]
  
}

# Extract posterior samples for occupancy parameters (n = nsamp)
psi_samp <- samples[subsamples, grep("beta_psi", colnames(samples))]

# Convert SpatRaster to a dataframe (one row for each cell, each column = layer)
# note: this removes cells/rows with NAs
psi_rasters_df <- as.data.frame(psi_rasters, cell = TRUE)
# Do the math (results in predictions on logit scale for each cell [row] and 
# posterior sample [column])
occprob <- as.matrix(psi_rasters_df[,-1]) %*% t(psi_samp) 
# Convert to probability scale (can take a few seconds)
occprob <- exp(occprob)/(1 + exp(occprob))
# Re-attach cell numbers
occprob <- cbind(cell = psi_rasters_df$cell, occprob)
