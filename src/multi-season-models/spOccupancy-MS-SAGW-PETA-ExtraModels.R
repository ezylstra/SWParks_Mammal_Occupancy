#------------------------------------------------------------------------------#
# Run a few extra models for PETA in SAGW
# TODO: generalize this part of the script later
#------------------------------------------------------------------------------#

model_stats_orig <- model_stats
occ_specs <- c("~ slope_z + I(slope_z^2) + (1 | site)",
               "~ boundary_z + (1 | site)",
               "~ slope_z + I(slope_z^2) + boundary_z + (1 | site)",
               "~ slope_z + I(slope_z^2) + boundary_z + years_z + (1 | site)",
               "~ boundary_z + years_z + boundary_z * years_z + (1 | site)")
det_specs <- "~ effort_z"
model_specs <- as.matrix(expand.grid(occ = occ_specs, 
                                     det = det_specs,
                                     KEEP.OUT.ATTRS = FALSE))

# Create model formulas with R syntax
as.formula.vect <- Vectorize(as.formula)
occ_formulas <- as.formula.vect(model_specs[,1])
det_formulas <- as.formula.vect(model_specs[,2])

# Run candidate models using spOccupancy package
source("src/multi-season-models/spOccupancy-MS-run-candidate-models.R")
# Note: this can take several minutes to run

# View summary table, ranked by WAIC
model_stats %>% arrange(waic)

best_index <- 3 # Slope + Boundary model

# Extract output and formulas from best model
best <- out_list[[best_index]]
best_psi_model <- model_specs[best_index, 1]
best_p_model <- model_specs[best_index, 2]

# Extract covariate names (with and without "_z" subscripts) from best model
psi_covs_z <- create_cov_list(best_psi_model)
psi_covs_z <- psi_covs_z[psi_covs_z != "years_z"]
p_covs_z <- create_cov_list(best_p_model)
p_covs_z <- p_covs_z[p_covs_z != "years_z"]
psi_covs <- psi_covs_z %>% str_remove_all(pattern = "_z")
p_covs <- p_covs_z %>% str_remove_all(pattern = "_z")

# Create table with summary stats that can be saved to file
occ_estimates <- parameter_estimates(model = best, 
                                     parameter = "occ",
                                     lower_ci = 0.025,
                                     upper_ci = 0.975)
det_estimates <- parameter_estimates(model = best, 
                                     parameter = "det",
                                     lower_ci = 0.025,
                                     upper_ci = 0.975)
occ_estimates <- occ_estimates %>%
  rename(Covariate = Parameter) %>%
  mutate(Parameter = "Occurrence", .before = "Covariate")
det_estimates <- det_estimates %>%
  rename(Covariate = Parameter) %>%
  mutate(Parameter = "Detection", .before = "Covariate")

source("src/multi-season-models/spOccupancy-MS-predictions.R")

# Can save any of the plots to file (example below):
plot_save <- plot_preds_mn_lastyr +
  theme_bw(base_size = 8)
plotname <- paste0("C:/Users/erin/Desktop/Mammals/",
                   PARK, "-", SPECIES, "-OccProbMN-SlopeBoundary-",
                   YEARS[length(YEARS)], ".jpg")
ggsave(filename = plotname,
       plot = plot_save,
       device = "jpeg",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)
