################################################################################
# Process results from a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-10-19
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
psi_covars <- unique(str_remove(colnames(cov_psi), "_z2|_z"))
# Unzip park rasters
park_folder <- paste0("data/covariates/rasters-", PARK, "/")
park_zip <- paste0("data/covariates/rasters-", PARK, ".zip")
if (length(list.files(park_folder)) == 0) {
  unzip(park_zip)
}

# Next steps:
# Identify rasters associated with covairates in psi_covars
# Standardize covariates in rasters and create quadratic layers where needed
# Combine rasters with beta posterior distributions
# Calculate the median value (and maybe the SD) for each cell 

#------------------------------------------------------------------------------#
# Plotting predicted probability of occupancy in last year
#------------------------------------------------------------------------------#

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

# Extract 1000 (of 3000) posterior samples for PAO estimates
pao <- samples[seq(1, 3000, by = 3), grep("PAO", colnames(samples))]

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
  mutate(yr2022 = int + 5 * slope)

ggplot() +
  geom_segment(trends,
               mapping = aes(x = 2017, xend = 2022, y = int, yend = yr2022),
               size = 0.3, col = "gray") +
  geom_segment(trends,
               mapping = aes(x = 2017, xend = 2022, 
                             y = median(int), yend = median(yr2022)),
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
preds_median <- data.frame(Year = yr_predict + 2017, 
                           Occupancy = apply(preds_prob, 2, median))
# Add column to identify MCMC iteration
preds_prob <- cbind(preds_prob, iter = 1:nrow(preds_prob))
# Convert predictions to long form for ggplot
preds_long <- preds_prob %>%
  as.data.frame %>%
  pivot_longer(cols = !iter,
               names_to = "Year",
               values_to = "Occupancy") %>%
  mutate(Year = as.numeric(Year) + 2017) %>%
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
  scale_y_continuous(limits = c(0.4, 0.9)) +
  scale_x_continuous(breaks = day_nums_axis, labels = format(days_axis, "%m-%d")) + 
  labs(x = "Date", y = "Probability of occupancy in Year 1") +
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

