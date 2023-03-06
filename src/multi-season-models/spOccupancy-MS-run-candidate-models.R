################################################################################
# Run candidate multi-season models (with spOccupancy) and gather model 
# diagnostics and statistics for comparisons

# In most instances, this will be called from another script (something like:
# src/multi-season-models/PARK/spOccupancy-PARK-SPECIES_YEARS.R)

# ER Zylstra
# Updated 2023-03-06
################################################################################

# Output from candidate models will be stored in a list (out_list)
# Model summaries/diagnostics will be stored in a dataframe (model_stats)

#------------------------------------------------------------------------------#
# Run models
#------------------------------------------------------------------------------#

# Running candidate models (specified in model_specs). Notes:
  # Not doing 4-fold cross validation for each model since this takes a while.
    # Could easily go back and do this once selecting a model for inferences (and
    # keeping only that model in a new candidate set)
  # Not specifying priors, but using defaults which are N(0, var = 2.72)
  # Not specifying initial values -- by default they come from priors

ar1 <- ifelse(TIME_RE_OCC == "AR1", TRUE, FALSE)

set.seed(2023)
out_list <- list()
for (i in 1:nrow(model_specs)) {
  message("Running model ", i, " (of ", nrow(model_specs), ").")
  
  if (SITE_RE_OCC == "spatial") {
    
    # set n.omp.threads to something > 1 
    num_cores <- parallel::detectCores() - 2
    if (num_cores > 8) {
      num_cores <- 8
    }
    n.omp.threads <- num_cores
    
    out_list[[i]] <- 
      tryCatch(
        {
        stPGOcc()
        },
        error = function(e){
          message("An error occurred attempting to run model ", i, ":\n", e)
          return(e)
        }
      )

  } else {
  
    out_list[[i]] <- 
      tryCatch(
        {
        tPGOcc(occ.formula = occ_formulas[[i]],
                     det.formula = det_formulas[[i]], 
                     data = data_list, 
                     n.batch = N_BATCH,
                     batch.length = BATCH_LENGTH,
                     ar1 = ar1,
                     n.omp.threads = 1, 
                     verbose = TRUE, 
                     n.burn = N_BURN, 
                     n.thin = N_THIN, 
                     n.chains = N_CHAINS,
                     n.report = 50)
        },
        error = function(e){
          message("An error occurred attempting to run model ", i, ":\n", e)
          return(e)
        }
      )
  }
}

# Print error messages if any models didn't run (can identify problematic 
# models because they're just a list with 2 elements [message, call])
for (i in 1:length(out_list)) {
  if (length(out_list[[i]]) == 2) {
    message("An error occurred attempting to run model ", i, ": ", out_list[[i]][1])
  }
}

#------------------------------------------------------------------------------#
# Extract a few model diagnostics
#------------------------------------------------------------------------------#

# Assess convergence with Rhat value (maximum across all parameters, including 
# random effects). Would like to see a value < 1.05
max_rhat <- lapply(out_list, function(x) 
  if(length(x) == 2) NA else max(unlist(x$rhat)))

# Assess effective sample sizes (ESS; minimum across all parameters, including
# random effects). Would like to see a value > 400.
min_ESS <- lapply(out_list, function(x) 
  if(length(x) == 2) NA else min(unlist(x$ESS)))

# Posterior predictive checks (how well does our model fit the data?)
  # From vignette: binning the data across sites (group = 1) may help reveal 
  # whether the model fails to adequately represent variation in occurrence and 
  # detection probability across space, while binning the data across replicates 
  # (group = 2) may help reveal whether the model fails to adequately represent 
  # variation in detection probability across the different replicate surveys. 

# # Create a function to get a Bayesian p-value from posterior predictive checks
# # (want to see values between 0.1 and 0.9)
# bayes.p <- function(object, digits = 3) {
#   format(round(mean(object$fit.y.rep > object$fit.y), digits), 
#          scientific = FALSE,
#          nsmall = digits)
# }
# # Using Freeman-Tukey statistic for now, but could use chi-squared instead 
# ppc.sites <- lapply(out_list, FUN = ppcOcc, fit.stat = "freeman-tukey", group = 1)
# ppc.reps <- lapply(out_list, FUN = ppcOcc, fit.stat = "freeman-tukey", group = 2)
# # These take a few seconds to run

#------------------------------------------------------------------------------#
# Model selection tools
#------------------------------------------------------------------------------#

# Get WAIC values for each model (used to compare models; lower is better)
waic <- lapply(out_list, function(x) 
  if(length(x) == 2) NA else waicOcc(x))

# Table with model stats
model_stats <- data.frame(model_no = 1:length(out_list),
                          psi = model_specs[,"occ"],
                          det = model_specs[,"det"],
                          max.rhat = round(unlist(max_rhat), 2),
                          min.ESS = round(unlist(min_ESS)),
                          # ppc.sites = unlist(lapply(ppc.sites, FUN = bayes.p)),
                          # ppc.reps = unlist(lapply(ppc.reps, FUN = bayes.p)),
                          waic = round(sapply(waic, "[", "WAIC"), 2))
