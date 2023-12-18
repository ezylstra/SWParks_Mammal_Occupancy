################################################################################
# Develop formulas for occurrence and detection part of candidate models

# In most instances, this will be called from another script (something like:
# src/multi-season-models/PARK/spOccupancy-PARK-SPECIES_YEARS.R)

# ER Zylstra
# Updated 2023-05-26
################################################################################

#------------------------------------------------------------------------------#
# Create formulas for occurrence part of candidate models
#------------------------------------------------------------------------------#

if (exists("OCC_MODELS")) {
  occm <- list()
  for (i in 1:length(OCC_MODELS)) {
    occm[[i]] <- paste(covariates$formula[covariates$short_name %in% OCC_MODELS[[i]] & 
                                            covariates$park %in% c("all", PARK)],
                       collapse = " + ")
  }
  occm <- unlist(occm)
} else {
  occm <- NA
}

if (OCC_NULL) {occm <- c("1", occm)}
occm <- occm[!is.na(occm)]
occm <- paste0("~ ", occm) 

# Remove any that elements with no covariates (that aren't ~ 1)
occm <- occm[occm != "~ "]

# If we want to include unstructured site random effects, add random site
# intercepts using lme4 syntax
if (SITE_RE_OCC == "unstructured") {
  occm <- paste0(occm, " + (1 | site)")
}

# If we want to include unstructured temporal random effects, add random yearly
# intercepts using lme4 syntax
if (TIME_RE_OCC == "unstructured") {
  occm <- paste0(occm, " + (1 | years)")
}

#------------------------------------------------------------------------------#
# Create formulas for detection part of candidate models
#------------------------------------------------------------------------------#

if (exists("DET_MODELS")) {
  detm <- list()
  for (i in 1:length(DET_MODELS)) {
    detm[[i]] <- paste(covariates$formula[covariates$short_name %in% DET_MODELS[[i]] &
                                            covariates$park %in% c("all", PARK)],
                       collapse = " + ")
  }
  detm <- unlist(detm)
} else {
  detm <- NA
}

if (DET_NULL) {detm <- c("1", detm)}
detm <- detm[!is.na(detm)]
detm <- paste0("~ ", detm) 

# Remove any that elements with no covariates (that aren't ~ 1)
detm <- detm[detm != "~ "]

# If we want to include unstructured site random effects, add random site
# intercepts using lme4 syntax
if (SITE_RE_DET == "unstructured") {
  detm <- paste0(detm, " + (1 | site)")
}

# If we want to include unstructured temporal random effects, add random yearly
# intercepts using lme4 syntax
if (TIME_RE_DET == "unstructured") {
  detm <- paste0(detm, " + (1 | years)")
}

#------------------------------------------------------------------------------#
# Combine occurrence and detection formulas to create candidate model 
# specifications
#------------------------------------------------------------------------------#

# Create a matrix that contains all combinations of occurrence and detection 
# covariates for candidate models
model_specs <- as.matrix(expand.grid(occ = occm, 
                                     det = detm,
                                     KEEP.OUT.ATTRS = FALSE))

# Create model formulas with R syntax
as.formula.vect <- Vectorize(as.formula)
occ_formulas <- as.formula.vect(model_specs[,1])
det_formulas <- as.formula.vect(model_specs[,2])
