################################################################################
# Process results from a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-10-13
################################################################################

library(dplyr)
library(stringr)
library(terra)
library(jagsUI)
library(MCMCvis)

#------------------------------------------------------------------------------#
# Load the workspace that contains the JAGS model
#------------------------------------------------------------------------------#

PARK <- "SAGW"
SPECIES <- "LECA"
DATE <- "2022-10-13"
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
# Extract posterior samples, and set parameters for summary tables, figures
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

# Function to calculate the proportion of posterior samples that are > 0 for 
# median values > 0 or < 0 for median values < 0
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
# TO DO:

# Add example code for plotting or calculating marginal effects (in 
# MSoccupancy-sagw-leca.R)

# Add code to estimate trends

# Add code to plot predictions

#------------------------------------------------------------------------------#