################################################################################
# Develop formulas for occurrence and detection part of candidate models

# In most instances, this will be called from another script (something like:
# src/single-seasons/models/YEAR/spOccupancy-PARK-SPECIES_YEAR.R)

# ER Zylstra
# Updated 2023-02-02
################################################################################

#------------------------------------------------------------------------------#
# Create formulas for occurrence part of candidate models
#------------------------------------------------------------------------------#

occm1 <- covariates$formula[covariates$short_name %in% OCC_MODELS1]
occm2 <- list()
for (i in 1:length(OCC_MODELS2)) {
  occm2[[i]] <- paste(covariates$formula[covariates$short_name %in% OCC_MODELS2[[i]]],
                      collapse = " + ")
}
occ_specs <- c(occm1, unlist(occm2))
if (OCC_NULL) {
  occ_specs <- c("1", occ_specs) 
}
occ_specs <- paste0("~ ", occ_specs) 

#------------------------------------------------------------------------------#
# Create formulas for detection part of candidate models
#------------------------------------------------------------------------------#

detm1 <- covariates$formula[covariates$short_name %in% DET_MODELS1]
detm2 <- list()
for (i in 1:length(DET_MODELS2)) {
  detm2[[i]] <- paste(covariates$formula[covariates$short_name %in% DET_MODELS2[[i]]],
                      collapse = " + ")
}
det_specs <- c(detm1, unlist(detm2))
if (DET_NULL) {
  det_specs <- c("1", det_specs) 
}
det_specs <- paste0("~ ", det_specs) 

#------------------------------------------------------------------------------#
# Combine occurrence and detection formulas to create candidate model 
# specifications
#------------------------------------------------------------------------------#

# Create a matrix that contains all combinations of occurrence and detection 
# covariates for candidate models
model_specs <- as.matrix(expand.grid(occ = occ_specs, 
                                     det = det_specs,
                                     KEEP.OUT.ATTRS = FALSE))

# Create model formulas with R syntax
as.formula.vect <- Vectorize(as.formula)
occ_formulas <- as.formula.vect(model_specs[,1])
det_formulas <- as.formula.vect(model_specs[,2])
