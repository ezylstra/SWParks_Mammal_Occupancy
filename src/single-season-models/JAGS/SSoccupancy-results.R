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
SPECIES <- "ODHE"
YEAR <- 2022
DATE <- "2023-01-12"
DATE <- str_remove_all(DATE, "-")
output_file <- paste0("output/models-JAGS/",
                      tolower(PARK), "-",
                      tolower(SPECIES), "-", 
                      YEAR, "-slope-",
                      DATE,".Rdata")
load(output_file)

# Check attributes of model (species, park, year, covariates)
cat(model_description) 

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
  } else {
    raster_covars <- psi_covars
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
psi_rasters_df <- as.data.frame(psi_rasters, cell = TRUE)
# Remove row with any covariates equal to NA
psi_rasters_df$nNAs <- apply(psi_rasters_df[,-1], 1, function(x) sum(is.na(x)))
psi_rasters_df <- psi_rasters_df %>%
  filter(nNAs == 0) %>%
  select(-nNAs)
# Do the math (results in predictions on logit scale for each cell [row] and 
# posterior sample [column])
occprob <- as.matrix(psi_rasters_df[,-1]) %*% t(psi_samp) 
# Convert to probability scale (can take a few seconds)
occprob <- exp(occprob)/(1 + exp(occprob))
# Re-attach cell numbers
occprob <- cbind(cell = psi_rasters_df$cell, occprob)

# Draw values of latent occupancy state (z)
z <- matrix(rbinom(length(occprob[,-1]), 1, occprob[,-1]),
            nrow = nrow(occprob[,-1]), ncol = ncol(occprob[,-1]))

# Calculate proportion of area occupied (PAO)
PAO_park <- apply(z, 2, function(x) sum(x)/length(x))
# PAO_park is a vector that has the estimated proportion of cells occupied
# for each posterior sample (n = nsamp)

# Create raster with median, SD occupancy probabilities across the park
# (median and SD)
occ_stats <- data.frame(cell = occprob[,1], 
                        median = apply(occprob[,-1], 1, median),
                        sd = apply(occprob[,-1], 1, sd))
preds_median <- rast(rast_int)
preds_median[occ_stats[,1]] <- occ_stats[,2]
names(preds_median) <- "occupancy_prob_median"
preds_sd <- rast(rast_int)
preds_sd[occ_stats[,1]] <- occ_stats[,3]
names(preds_sd) <- "occupancy_prob_sd"

# Look at estimated occupancy probabilities (and SDs)
par(mfrow = c(2,1))
plot(preds_median, main = "Median")
plot(preds_sd, main = "SD")
par(mfrow = c(1,1))

# Remove unnecessary objects
remove <- c("psi_list", "psi_cont", "veg_rast", "psi_rasters", "psi_rasters_df", 
            "psi_samp", "veg2", "veg3", "z")
remove <- remove[sapply(remove, exists)]
rm(list = remove)

#------------------------------------------------------------------------------#
# Calculate marginal effects for a continuous covariate (example)
#------------------------------------------------------------------------------#

# Estimate how occupancy varies as a function of distance from some
# point-of-interest (POI; includes buildings, parking lots)

# Notes: logit is the log of the odds
# Odds are the probability event happens / probability event doesn't happen
# So logit(psi[i]) = log(psi[i] / (1-psi[i])) = log(odds)
# Odds a site is occupied at the mean distance from a POI = exp(beta_psi0)
# Odds a site is occupied with a 1-SD increase in distance = 
# exp(beta_psi0 + beta_psi1 * 1) = exp(beta_psi0)*exp(beta_psi1)
# So the odds will change by a FACTOR of exp(beta_psi1):
# Odds[dist+1SD] = Odds[mean dist] * exp(beta_psi1)

# For sagw-leca-2020-20221213.Rdata model, pois effect is beta_psi[3]
# logit(psi[i]) = beta_psi0 + beta_psi[3] * dist[i]
beta_psi <- samples[,"beta_psi[3]"]
change <- exp(beta_psi)
mean(change) # 7.53
median(change) # 5.36
quantile(change, probs = c(lcl, ucl)) #1.79, 26.43
# The odds a site is colonized are estimated to be 7.5 times higher (750% higher) 
# for each 1-SD increase in distance (95% CI = 1.8 - 26.4 times higher) 
# [assuming other covariates held constant]

#------------------------------------------------------------------------------#
# Plot marginal effects for a continuous covariate (example)
#------------------------------------------------------------------------------#

# Estimate how detection probability changes with day-of-year

# For sagw-leca-2020-20221213.Rdata model, linear and quadratic effect of day
# are beta_p[1:2]
# logit(p[i]) = beta_p0 + beta_p[1] * doy[i] + beta_p[2] * doy[i] * doy[i]

# Generate a vector of days that span the range that occurred during study
day_nums <- seq(min(surveys$day), max(surveys$day), length = 100)
# Standardize these values (Make sure to grab the mean/SD from the dataframe 
# where values were standardized in SSoccupancy-generic.R. For day number, this 
# was in the occasions dataframe)
day_nums_z <- (day_nums - mean(occasions$mid_yday)) / sd(occasions$mid_yday)
# Create a matrix of covariate values (including the intercept [1])
X_p <- cbind(int = 1, doy = day_nums_z, doy2 = day_nums_z *day_nums_z)
# Create a matrix with posterior samples for the parameters we need
betas_p <- samples[,c("beta_p0", "beta_p[1]", "beta_p[2]")]
# Generate a matrix of predicted values on the logit scale
# Matrix dimensions = 100 x 3000 (3000 predicted values for each day in seq)
pred_logit_p <- X_p %*% t(betas_p)
# Backtransform to the probability scale
pred_prob_p <- exp(pred_logit_p) / (1 + exp(pred_logit_p))
# Calculate mean, median, and credible interval for each value
mean_p <- apply(pred_prob_p, 1, mean)
median_p <- apply(pred_prob_p, 1, median)
cri_p <- apply(pred_prob_p, 1, quantile, probs = c(lcl, ucl)) 
# Put objects in dataframe for ggplot
plot_data <- data.frame(mean = mean_p,
                        median = median_p,
                        day_nums = day_nums,
                        lcl = cri_p[1,],
                        ucl = cri_p[2,])

# Plot predictions
# For day-of-year, transform day numbers we want on x-axis to date 
# (for easier interpretation)
day_nums_axis <- seq(min(day_nums), max(day_nums), by = 10)
days_axis <- parse_date_time(paste(2020, as.character(day_nums_axis)), 
                             orders = "yj")

ggplot() + 
  geom_ribbon(plot_data,
              mapping = aes(x = day_nums, ymin = lcl, ymax = ucl), 
              fill = "gray70") +
  geom_line(plot_data,
            mapping = aes(x = day_nums, y = median)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = day_nums_axis, labels = format(days_axis, "%m-%d")) + 
  labs(x = "Date", y = "Detection probability") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

#------------------------------------------------------------------------------#
# Calculate marginal effects for a categorical covariate (example)
#------------------------------------------------------------------------------#

# Estimate how vegetation class affects probabilities of occupancy

# For sagw-leca-2020-20221213.Rdata model, vegclass1 (low gradient desert) is 
# the reference level, effect of vegclass2 (low hillslope, north-facing) is 
# beta_psi[4], and effect of vegclass3 (medium-high gradient) is beta_psi[5]

# Probability of occupancy in vegclass1 (at mean values of other covariates)
beta_psi0 <- samples[,"beta_psi0"]
vegclass1 <- exp(beta_psi0)/(1 + exp(beta_psi0)) 
mean(vegclass1) # 0.39
median(vegclass1) # 0.39
quantile(vegclass1, probs = c(lcl, ucl)) #0.12, 0.72
# Probability of occupancy = 0.39 (95% CI = 0.12, 0.72)

# Probability of occupancy in vegclass2 (at mean values of other covariates)
beta_veg2 <- samples[,"beta_psi[4]"]
vegclass2L <- beta_psi0 + beta_veg2
vegclass2 <- exp(vegclass2L)/(1 + exp(vegclass2L)) 
mean(vegclass2) # 0.55
median(vegclass2) # 0.56
quantile(vegclass2, probs = c(lcl, ucl)) #0.15, 0.89
# Probability of occupancy = 0.56 (95% CI = 0.15, 0.89)

# Probability of occupancy in vegclass3 (at mean values of other covariates)
beta_veg3 <- samples[,"beta_psi[5]"]
vegclass3L <- beta_psi0 + beta_veg3
vegclass3 <- exp(vegclass3L)/(1 + exp(vegclass3L)) 
mean(vegclass3) # 0.07
median(vegclass3) # 0.05
quantile(vegclass3, probs = c(lcl, ucl)) #0.00, 0.27
# Probability of occupancy = 0.05 (95% CI = 0.00, 0.27)
