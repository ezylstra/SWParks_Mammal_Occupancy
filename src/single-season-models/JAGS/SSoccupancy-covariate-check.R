################################################################################
# Carries out a series of checks on the covariate structure of a single-season
# model specified in SSoccupancy_wrapper.R

# ER Zylstra
# Updated 2023-01-18
################################################################################

# Produce error messages if any quadratic effects are specified incorrectly

  # Initial occupancy probability (PSI)
  if (!all(is.na(PSI_QUADS)) & any(!PSI_QUADS %in% COVARS_PSI)) {
    message("WARNING: Not all covariates listed in PSI_QUADS appear in COVARS_PSI.")
  }
  
  # Detection probability (P)
  if (!all(is.na(P_QUADS)) & any(!P_QUADS %in% COVARS_P)) {
    message("WARNING: Not all covariates listed in P_QUADS appear in COVARS_P.")
  }

# Check that linear models look as expected:
  
  # Initial occupancy probability (PSI)
  psi_formula <- ifelse(is.na(COVARS_PSI), 1, COVARS_PSI)
  if (!all(is.na(PSI_QUADS))) {
    psi_formula <- sort(c(psi_formula, paste0(PSI_QUADS, "2")))
  }
  
  # Detection probability (P)
  p_formula <- ifelse(is.na(COVARS_P), 1, COVARS_P)
  if (!all(is.na(P_QUADS))) {
    p_formula <- sort(c(p_formula, paste0(P_QUADS, "2")))
  }
  
  psi_m <- paste0("psi ~ ", paste(psi_formula, collapse = " + "))
  p_m <- paste0("p ~ ", paste(p_formula, collapse = " + "))

  message("Check that the models below look correct before proceeding!")
  cat(paste(psi_m, p_m, "", sep = "\n")) 
  
  species_descript <- paste0("Species: ", SPECIES)
  park_descript <- paste0("Park: ", PARK)
  yr_descript <- paste0("Year: ", YEAR)
  
  model_description <- paste(species_descript,
                             park_descript,
                             yr_descript,
                             psi_m, p_m, sep = "\n") 
  