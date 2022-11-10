################################################################################
# Process results from a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-11-10
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
DATE <- "2022-10-19"
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

# Identify the number of posterior samples we want to use for figures (spatial
# predictions and trends)
nsamp <- 1000

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

#------------------------------------------------------------------------------#
# Plotting predicted probability of occupancy in first year
#------------------------------------------------------------------------------#

# Essentially, we want to create a heat map with values of psi for each cell

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
subsamples <- floor(seq(1, nrow(samples), length = nsamp))
psi_samp <- samples[subsamples, grep("beta_psi", colnames(samples))]

# Convert SpatRaster to a dataframe (with one row for each cell, each column = layer)
rast_final_df <- as.data.frame(rast_final, cell = TRUE) #130193 rows (removed rows with NAs)
# Do the math (results in predictions on logit scale for each cell [row] and 
# posterior sample [column])
preds_df_logit <- as.matrix(rast_final_df[,-1]) %*% t(psi_samp) 
# Convert to probability scale
preds_df <- exp(preds_df_logit)/(1 + exp(preds_df_logit))
# Re-attach cell numbers
preds_df <- cbind(cell = rast_final_df$cell, preds_df)
# Convert back to a SpatRaster
preds_raster <- rast(rep(rast_final[[1]], nsamp))
names(preds_raster) <- paste("sample", subsamples)
preds_raster[preds_df[,1]] <- preds_df[,-1]
  # Check:
  # plot(preds_raster[[1:2]])
  # psi_samp[1:2,]

# Calculate the median value in each cell (our best estimate of initial occ)
preds_median <- median(preds_raster)
plot(preds_median)
# Calculate the SD across posterior samples in each cell
preds_sd <- stdev(preds_raster, pop = FALSE)
plot(preds_sd)
# Note: pop = TRUE computes the population SD (denom = n), pop = FALSE computes
# the sample standard deviation (denom = n-1)

# TODO: add options to save these rasters and/or plots?

#------------------------------------------------------------------------------#
# Plotting predicted probability of occupancy in last year
#------------------------------------------------------------------------------#

# Identify covariates for extinction and colonization
# gam_covars <- unique(str_remove_all(colnames(cov_gam), "_z2|_z|[:digit:]"))
# eps_covars <- unique(str_remove_all(colnames(cov_eps), "_z2|_z|[:digit:]"))
# Will want to identify interactions where they occur, but this is difficult
# Maybe better to save objects from the wrapper script (eg, EPS_INT1)?
# Could then also bring in COVARS_PSI, COVARS_GAM, COVARS_EPS....
COVARS_GAM <- c("elev", "monsoon_ppt")
COVARS_EPS <- c("elev", "monsoon_ppt")
EPS_INT1 <- c("elev", "monsoon_ppt")

# We already have all the spatial covars unzipped from section above

# List of potential seasonal covariates
seas_covariates <- c("monsoon_ppt") # Will add more in time

# Identify all the rasters we'll need for extinction and colonization
seas_gam <- str_subset(COVARS_GAM, paste(seas_covariates, collapse = "|"))
seas_eps <- str_subset(COVARS_EPS, paste(seas_covariates, collapse = "|"))
seas_both <- unique(c(seas_gam, seas_eps))
 
# Need to unzip and load rasters with weather data (or eventually raster with 
# any annually varying covs)
weather_folder <- "data/covariates/weather-derived-rasters/"
weather_zip <- "data/covariates/weather-derived.zip"
# Unzip weather folder first, if necessary
if (length(list.files(weather_folder)) == 0) {
  unzip(weather_zip, overwrite = TRUE)
}
# List files in weather folder
weather_files <- list.files(weather_folder, full.names = TRUE)

# For each seasonal covariate we'll need, create a list of rasters
for (cov in seas_both) {
  cov_files <- weather_files[str_detect(weather_files, cov)]
  # Remove rasters associated with periods outside the years of interest
  # (eg, monsoon rainfall in year x could explain transitions between years x and x+1)
  cov_yrs <- paste0(as.character(year), collapse = "|")
  cov_files <- cov_files[str_detect(cov_files, cov_yrs)]

  
  ####################################################
  # Everything in the rest of this section needs work.....
  
  # Load each raster and compile into a list
  cov_list <- list()
  for (i in 1:length(cov_files)) {
    cov_list[[i]] <- rast(cov_files[i])
    names(cov_list[[i]]) <- paste(cov, year[i], sep = "_")  ##### is this a good idea?
  }  
  
  # Crop rasters to park boundary
  cov_list <- lapply(cov_list, FUN = crop, ext(boundary))
  
  # Scale values in each raster using means, SDs from sitetrans column (monsoon_ppt)
  
} # Close loop for each seasonal covariate 
  

# Then for eps and col separately...{}

  # add spatial covar raster to list, IF NEEDED (and in the same order as they appear in the model)
  
  # Create a raster layer with just ones for the intercept
  rast_final <- c(rast(rast_final[[1]], vals = 1), rast_final)
  names(rast_final)[[1]] <- "int"

  # Combine all rasters in list to a spatraster
  
  # Convert to a dataframe
  
  # Create season-specific dataframes with intercept, spatial covs, and seasonal cov [1 col per cov]
  
  # If interactions, create new column
  
  # Do linear algebra to get eps[,1], eps[,2], ...gam[,1], gam[,2], etc...
  
# Close eps and gam loops
  
# Then create a seasonal loop
  
  # draw values of z[,1] from preds_df MAKE SURE CELL #s MATCH UP!  
  for (i in 1:nyears) {
    # do element-wise math with z[,i], gam[,i], eps[,i] to get Ez[,i]
    # draw z[,i+1] from Ez[,i]
  }

# General approach:
# For each MCMC iteration:
# Estimate initial occupancy probability for each cell (section above)
# Draw a value latent occupancy state for each cell, z[i,1] from a Bernoulli
# Calculate gam[i,1] and eps[i,1] from beta_gam, beta_eps, and cell covariate values
# Calculate each cell's prob of occupancy in yr 2 based on z[i, 1], gam[i, 1], eps[i, 1]
# Cycle through years.

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
  labs(x = "Year", y = "logit(Proportion of sites occupied)")

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
  labs(y = "Proportion of sites occupied")

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

