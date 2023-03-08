#------------------------------------------------------------------------------#
# Run a few extra models for URCI in SAGW
# TODO: generalize this part of the script later
#------------------------------------------------------------------------------#

model_stats_orig <- model_stats
# occ_specs <- c("~ boundary_z + years_z + (1 | site)",
#                "~ boundary_z + elev_z + I(elev_z^2) + years_z + (1 | site)",
#                "~ boundary_z + wash_z + vegclass2 + vegclass3 + years_z + (1 | site)",
#                "~ boundary_z + wash_z + vegclass2 + vegclass3 + elev_z + I(elev_z^2) + years_z + (1 | site)")
# occ_specs <- c("~ boundary_z + (1 | site)",
#                "~ boundary_z + years_z + (1 | site)",
#                "~ boundary_z + monsoon_ppt_z + (1 | site)",
#                "~ boundary_z + ppt10_z + (1 | site)",
#                "~ wash_z + vegclass2 + vegclass3 + monsoon_ppt_z + (1 | site)",
#                "~ wash_z + vegclass2 + vegclass3 + ppt10_z + (1 | site)")
occ_specs <- c("~ boundary_z + ppt10_z + (1 | site)",
               "~ wash_z + vegclass2 + vegclass3 + ppt10_z + (1 | site)",
               "~ boundary_z + wash_z + vegclass2 + vegclass3 + ppt10_z + (1 | site)")
# det_specs <- c("~ 1", 
#                "~ 1 + (1 | years)")
det_specs <- c("~ 1")
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
