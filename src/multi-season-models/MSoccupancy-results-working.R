################################################################################
# Process results from a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-10-19
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(terra)
library(jagsUI)
library(MCMCvis)

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
# Calculate marginal effects for continuous covariate (example)
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
# Plot marginal effects for continuous covariate (example)
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

# Plot predictions
  # For day-of-year, transform day number to date (for easier interpretation)
  day_nums_axis <- seq(min(day_nums), max(day_nums), by = 10)
  days_axis <- parse_date_time(paste(2022, as.character(day_nums_axis)), 
                          orders = "yj")
par(mfrow=c(1,1), mar = c(3, 4, 1, 1))
plot(mean_p ~ day_nums, type = "l", ylim = c(0, 1), axes = FALSE,
     xlab = "", ylab = "")
  axis(1, at = c(par("usr")[1], par("usr")[2]), tck = FALSE, labels = FALSE)
  axis(1, at = day_nums_axis, labels = format(days_axis, "%m-%d"), 
       tcl = -0.25, mgp = c(1.5, 0.4, 0))
  axis(2, at = c(par("usr")[3], par("usr")[4]), tck = FALSE, labels = FALSE)
  axis(2, at = seq(0, 1, by = 0.2), tcl = -0.25, las = 1, mgp = c(1.5, 0.5, 0),
       labels = c(paste0("0.", seq(0, 8, 2)), "1.0"))
  polygon(c(day_nums, rev(day_nums)), c(cri_p[1,], rev(cri_p[2,])), 
          col = rgb(0, 0, 0, 0.2), border = NA)
  mtext("Date (midpoint of occasion)", side = 1, line = 1.8)
  mtext("Probability of occupancy in Year 1", las = 0, side = 2, line = 2.5)

#------------------------------------------------------------------------------#
# Calculate marginal effects for categorical covariate (example)
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





#------------------------------------------------------------------------------#
# TO DO:

# Add code to estimate trends

# Add code to plot predictions

#------------------------------------------------------------------------------#