################################################################################
# Gather model diagnostic metrics and statistics for model comparisons

# In most instances, this will be called from another script (something like:
# src/single-seasons/models/YEAR/spOccupancy-PARK-SPECIES_YEAR.R)

# ER Zylstra
# Updated 2023-01-26
################################################################################

# This script assumes output from multiple spOccupancy models are stored in 
# out_list

#------------------------------------------------------------------------------#
# Extract a few model diagnostics
#------------------------------------------------------------------------------#

# Assess convergence with Rhat value (maximum across all parameters)
# (would like to see a value <1.05)
max_rhat <- lapply(out_list, function(x) max(c(x$rhat$beta, x$rhat$alpha)))

# Assess effective sample sizes (ESS; minimum across all parameters)
# (would like to see a value > 400)
min_ESS <- lapply(out_list, function(x) min(c(x$ESS$beta, x$ESS$alpha)))

# Posterior predictive checks (how well does our model fit the data?)
  # From vignette: binning the data across sites (group = 1) may help reveal 
  # whether the model fails to adequately represent variation in occurrence and 
  # detection probability across space, while binning the data across replicates 
  # (group = 2) may help reveal whether the model fails to adequately represent 
  # variation in detection probability across the different replicate surveys. 

# Create a function to get a Bayesian p-value from posterior predictive checks
# (want to see values between 0.1 and 0.9)
bayes.p <- function(object, digits = 3) {
  format(round(mean(object$fit.y.rep > object$fit.y), digits), 
         scientific = FALSE,
         nsmall = digits)
}
# Using Freeman-Tukey statistic for now, but could use chi-squared instead 
ppc.sites <- lapply(out_list, FUN = ppcOcc, fit.stat = "freeman-tukey", group = 1)
ppc.reps <- lapply(out_list, FUN = ppcOcc, fit.stat = "freeman-tukey", group = 2)
# These take a few seconds to run

#------------------------------------------------------------------------------#
# Model selection tools
#------------------------------------------------------------------------------#

# Get WAIC values for each model (used to compare models; lower is better)
waic <- lapply(out_list, waicOcc)

# Get deviance stat from k-fold cross validation to compare predictive 
# performance (again, lower is better)
deviances <- unlist(sapply(out_list, "[", "k.fold.deviance"))

# Table with model stats
model_stats <- data.frame(model_no = 1:length(out_list),
                          psi = model_specs[,"occ"],
                          det = model_specs[,"det"],
                          max.rhat = round(unlist(max_rhat), 2),
                          min.ESS = round(unlist(min_ESS)),
                          ppc.sites = unlist(lapply(ppc.sites, FUN = bayes.p)),
                          ppc.reps = unlist(lapply(ppc.reps, FUN = bayes.p)),
                          waic = round(sapply(waic, "[", "WAIC"), 2),
                          k.fold.dev = round(deviances, 2))
