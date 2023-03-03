################################################################################
# Develop formulas for occurrence and detection part of candidate models

# In most instances, this will be called from another script (something like:
# src/multi-season-models/PARK/spOccupancy-PARK-SPECIES_YEARS.R)

# ER Zylstra
# Updated 2023-03-03
################################################################################

#------------------------------------------------------------------------------#
# Create formulas for occurrence part of candidate models
#------------------------------------------------------------------------------#

if (exists("OCC_MODELS1")) {
  occm1 <- covariates$formula[covariates$short_name %in% OCC_MODELS1]
} else {
  occm1 <- NA
}
if (exists("OCC_MODELS2")) {
  occm2 <- list()
  for (i in 1:length(OCC_MODELS2)) {
    occm2[[i]] <- paste(covariates$formula[covariates$short_name %in% OCC_MODELS2[[i]]],
                        collapse = " + ")
  }
  occm2 <- unlist(occm2)
} else {
  occm2 <- NA
}
occ_specs <- c(occm1, occm2)
if (OCC_NULL) {occ_specs <- c("1", occ_specs)}
occ_specs <- occ_specs[!is.na(occ_specs)]
occ_specs <- paste0("~ ", occ_specs) 

# If we want to include unstructured site random effects, add random site
# intercepts using lme4 syntax
if (SITE_RE_OCC == "unstructured") {
  occ_specs <- paste0(occ_specs, " + (1 | site)")
}

# If we want to include unstructured temporal random effects, add random yearly
# intercepts using lme4 syntax
if (TIME_RE_OCC == "unstructured") {
  occ_specs <- paste0(occ_specs, " + (1 | years)")
}

#------------------------------------------------------------------------------#
# Create formulas for detection part of candidate models
#------------------------------------------------------------------------------#

if (exists("DET_MODELS1")) {
  detm1 <- covariates$formula[covariates$short_name %in% DET_MODELS1]
} else {
  detm1 <- NA
}
if (exists("DET_MODELS2")) {
  detm2 <- list()
  for (i in 1:length(DET_MODELS2)) {
    detm2[[i]] <- paste(covariates$formula[covariates$short_name %in% DET_MODELS2[[i]]],
                        collapse = " + ")
  }
  detm2 <- unlist(detm2)
} else {
  detm2 <- NA
}
det_specs <- c(detm1, detm2)
if (DET_NULL) {det_specs <- c("1", det_specs)}
det_specs <- det_specs[!is.na(det_specs)]
det_specs <- paste0("~ ", det_specs) 

# If we want to include unstructured site random effects, add random site
# intercepts using lme4 syntax
if (SITE_RE_DET == "unstructured") {
  det_specs <- paste0(det_specs, " + (1 | site)")
}

# If we want to include unstructured temporal random effects, add random yearly
# intercepts using lme4 syntax
if (TIME_RE_DET == "unstructured") {
  det_specs <- paste0(det_specs, " + (1 | years)")
}

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
