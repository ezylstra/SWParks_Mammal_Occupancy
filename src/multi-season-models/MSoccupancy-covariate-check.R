################################################################################
# Carries out a series of checks on the covariate structure of a multi-season
# model specified in MSoccupancy_wrapper.R

# ER Zylstra
# Updated 2022-11-30
################################################################################

# Produce error messages if any quadratic effects or interactions are specified 
# incorrectly

  # Initial occupancy probability (PSI)
  if (!all(is.na(PSI_QUADS)) & any(!PSI_QUADS %in% COVARS_PSI)) {
    message("WARNING: Not all covariates listed in PSI_QUADS appear in COVARS_PSI.")
  }
  
  # Detection probability (P)
  if (!all(is.na(P_QUADS)) & any(!P_QUADS %in% COVARS_P)) {
    message("WARNING: Not all covariates listed in P_QUADS appear in COVARS_P.")
  }
  
  # Extinction probability (EPS)
  if (!all(is.na(EPS_QUADS)) & any(!EPS_QUADS %in% COVARS_EPS)) {
    message("WARNING: Not all covariates listed in EPS_QUADS appear in COVARS_EPS.")
  } 
  if (!exists("EPS_INT1")) {EPS_INT1 <- NA}
  if ((length(EPS_INT1) == 1 & !all(is.na(EPS_INT1))) | length(EPS_INT1) > 2) {
    message("WARNING: EPS_INT1 must contain only two covariate names.")
  }
  if (N_EPS_INTERACTS == 0 & !all(is.na(EPS_INT1))) {
    message("WARNING: EPS_INT1 cannot contain covariates when N_EPS_INTERACTS = 0.")
  }
  if (N_EPS_INTERACTS > 0 & all(is.na(EPS_INT1))) {
    message("WARNING: EPS_INT1 must contain covariates when N_EPS_INTERACTS = 1.")
  }  
  
  # Colonization probability (GAM)
  if (!all(is.na(GAM_QUADS)) & any(!GAM_QUADS %in% COVARS_GAM)) {
    message("WARNING: Not all covariates listed in GAM_QUADS appear in COVARS_GAM.")
  } 
  if (!exists("GAM_INT1")) {GAM_INT1 <- NA}
  if ((length(GAM_INT1) == 1 & !all(is.na(GAM_INT1))) | length(GAM_INT1) > 2) {
    message("WARNING: GAM_INT1 must contain only two covariate names.")
  }
  if (N_GAM_INTERACTS == 0 & !all(is.na(GAM_INT1))) {
    message("WARNING: GAM_INT1 cannot contain covariates when N_GAM_INTERACTS = 0.")
  }
  if (N_GAM_INTERACTS > 0 & all(is.na(GAM_INT1))) {
    message("WARNING: GAM_INT1 must contain covariates when N_GAM_INTERACTS = 1.")
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
  
  # Extinction probability (EPS)
  if (N_EPS_INTERACTS > 0) {
    for (i in 1:N_EPS_INTERACTS) {
      assign(paste0("eps_int",i),paste(get(paste0("EPS_INT", i)), collapse = "*"))
    }
    eps_formula_int <- paste(mget(str_subset(ls(), "eps_int")), collapse = " + ")
  }
  eps_formula <- ifelse(is.na(COVARS_EPS), 1, COVARS_EPS)
  if (!all(is.na(EPS_QUADS))) {
    eps_formula <- sort(c(eps_formula, paste0(EPS_QUADS, "2")))
  }
  if (exists("eps_formula_int")) {eps_formula <- c(eps_formula, eps_formula_int)}
  
  # Colonization probability (GAM)
  if (N_GAM_INTERACTS > 0) {
    for (i in 1:N_GAM_INTERACTS) {
      assign(paste0("gam_int",i),paste(get(paste0("GAM_INT", i)), collapse = "*"))
    }
    gam_formula_int <- paste(mget(str_subset(ls(), "gam_int")), collapse = " + ")
  }
  gam_formula <- ifelse(is.na(COVARS_GAM), 1, COVARS_GAM)
  if (!all(is.na(GAM_QUADS))) {
    gam_formula <- sort(c(gam_formula, paste0(GAM_QUADS, "2")))
  }
  if (exists("gam_formula_int")) {gam_formula <- c(gam_formula, gam_formula_int)}

  psi_m <- paste0("psi ~ ", paste(psi_formula, collapse = " + "))
  p_m <- paste0("p ~ ", paste(p_formula, collapse = " + "))
  eps_m <- paste0("eps ~ ", paste(eps_formula, collapse = " + "))
  gam_m <- paste0("gam ~ ", paste(gam_formula, collapse = " + "))
  
  message("Check that the models below look correct before proceeding!")
  cat(paste(psi_m, p_m, eps_m, gam_m, sep = "\n")) 
  
