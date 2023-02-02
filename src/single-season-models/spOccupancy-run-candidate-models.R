################################################################################
# Run candidate single-season models (with spOccupancy) and gather model 
# diagnostics and statistics for comparisons

# In most instances, this will be called from another script (something like:
# src/single-seasons/models/YEAR/spOccupancy-PARK-SPECIES_YEAR.R)

# ER Zylstra
# Updated 2023-02-02
################################################################################

# Output from candidate models will be stored in a list (out_list)
# Model summaries/diagnostics will be stored in a dataframe (model_stats)

#------------------------------------------------------------------------------#
# Run models
#------------------------------------------------------------------------------#

# Running candidate models (specified in model_specs). Notes:
  # Doing 4-fold cross validation for each model
  # Not specifying priors, but using defaults which are N(0, var = 2.72)
  # Not specifying initial values -- by default they come from priors
  # Running chains sequentially (n.omp.threads = 1) because vignette states
  # this only speeds things up in spatial models

set.seed(2023)
out_list <- list()
for (i in 1:nrow(model_specs)) {
  message("Running model ", i, " (of ", nrow(model_specs), ").")
  out <- PGOcc(occ.formula = occ_formulas[[i]],
               det.formula = det_formulas[[i]], 
               data = data_list, 
               # inits = inits, 
               # priors = priors, 
               n.samples = N_SAMPLES, 
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = N_BURN, 
               n.thin = N_THIN, 
               n.chains = N_CHAINS,
               k.fold = 4,
               k.fold.threads = 4,
               k.fold.seed = 2023) 
  out_list <- c(out_list, list(out))
}

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
