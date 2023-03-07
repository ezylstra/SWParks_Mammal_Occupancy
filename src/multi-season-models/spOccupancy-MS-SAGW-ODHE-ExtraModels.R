#------------------------------------------------------------------------------#
# Run a few extra models for ODHE in SAGW
# TODO: generalize this part of the script later
#------------------------------------------------------------------------------#

model_stats_orig <- model_stats
occ_specs <- c("~ wash_z + vegclass2 + vegclass3 + pois_z + (1 | site)",
               "~ wash_z + vegclass2 + vegclass3 + (1 | site)",
               "~ pois_z + years_z + pois_z * years_z + (1 | site)")
det_specs <- "~ day_z + I(day_z^2) + deploy_exp + effort_z + (1 | years)"
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
