################################################################################
# Extract covariates and run a multi-season model in JAGS

# Script is typically called from MSoccupancy-wrapper.R (which has already
# sourced MSoccupancy-prep-data.R)

# ER Zylstra
# Updated 2023-01-18
################################################################################

#------------------------------------------------------------------------------#
# Extract covariate values for each parameter 
# (covariates selected in MSoccupancy-wrapper.R)
#------------------------------------------------------------------------------#

covs_cont <- c(covs_cont, "effort", "day", "monsoon_ppt", "burn_severity_2011")
# Not including camera_new or deploy_exp in this list because we don't want to 
# standardize (0 indicates camera type used prior to 2022 and inexperienced 
# people deploying cameras)

# Initial occupancy (psi)
if (!all(is.na(COVARS_PSI))) {
  # Create vector of standardized covariate names
  covars_psi <- ifelse(COVARS_PSI %in% covs_cont,
                       paste0(COVARS_PSI, "_z"), COVARS_PSI)
  # Extract covariate values
  cov_psi <- spatial_covs %>%
    select(contains(covars_psi))
  # Add quadratics if needed
  if (!all(is.na(PSI_QUADS))) {
    psi_quads <- names(cov_psi)[str_detect(names(cov_psi), 
                                           paste(PSI_QUADS, collapse = "|"))]
    cov_psi[,paste0(psi_quads, "2")] <- cov_psi[,psi_quads] ^ 2
  }
  # Put columns in alphabetical order and convert to a matrix
  cov_psi <- cov_psi %>%
    select(order(colnames(.))) %>%
    as.matrix
}

# Detection probability (p)  
if (!all(is.na(COVARS_P))) { 
  # Create vector of standardized covariate names
  covars_p <- ifelse(COVARS_P %in% covs_cont,
                     paste0(COVARS_P, "_z"), COVARS_P)
  # Extract covariate values
  cov_p <- surveys %>%
    select(contains(covars_p))
  # Add quadratics if needed
  if (!all(is.na(P_QUADS))) {
    p_quads <- names(cov_p)[str_detect(names(cov_p), 
                                       paste(P_QUADS, collapse = "|"))]
    cov_p[,paste0(p_quads, "2")] <- cov_p[,p_quads] ^ 2
  }
  # Put columns in alphabetical order and convert to a matrix
  cov_p <- cov_p %>%
    select(order(colnames(.))) %>%
    as.matrix  
}

# Extinction probability (eps)  
if (!all(is.na(COVARS_EPS))) {
  # Create vector of standardized covariate names
  covars_eps <- ifelse(COVARS_EPS %in% covs_cont,
                       paste0(COVARS_EPS, "_z"), COVARS_EPS)  
  # Extract covariate values  
  cov_eps <- sitetrans %>%
    select(contains(covars_eps))
  # Add quadratics if needed
  if (!all(is.na(EPS_QUADS))) {
    eps_quads <- names(cov_eps)[str_detect(names(cov_eps), 
                                           paste(EPS_QUADS, collapse = "|"))]
    cov_eps[,paste0(eps_quads, "2")] <- cov_eps[,eps_quads] ^ 2
  }
  # Put columns in alphabetical order
  cov_eps <- cov_eps %>%
    select(order(colnames(.)))
  # Add interactions if needed
  if (N_EPS_INTERACTS > 0) {
    eps_interacts <- matrix(NA, nrow = 1, ncol = 2)
    for (i in 1:N_EPS_INTERACTS) {
      eps_interacts[i,] <- get(paste0("EPS_INT", i))
      eps_interacts[i,1] <- ifelse(eps_interacts[i,1] %in% covs_cont,
                                   paste0(eps_interacts[i,1], "_z"),
                                   eps_interacts[i,1])
      eps_interacts[i,2] <- ifelse(eps_interacts[i,2] %in% covs_cont,
                                   paste0(eps_interacts[i,2], "_z"),
                                   eps_interacts[i,2])    
      cov_eps[,paste0(eps_interacts[i,], collapse = "_")] <- 
        cov_eps[,eps_interacts[i,1]] * cov_eps[,eps_interacts[i,2]]
    }
  }
  # Convert to a matrix
  cov_eps <- as.matrix(cov_eps)
}

# Colonization probability (gam)  
if (!all(is.na(COVARS_GAM))) {   
  # Create vector of standardized covariate names
  covars_gam <- ifelse(COVARS_GAM %in% covs_cont,
                       paste0(COVARS_GAM, "_z"), COVARS_GAM)  
  # Extract covariate values  
  cov_gam <- sitetrans %>%
    select(contains(covars_gam))
  # Add quadratics if needed
  if (!all(is.na(GAM_QUADS))) {
    gam_quads <- names(cov_gam)[str_detect(names(cov_gam), 
                                           paste(GAM_QUADS, collapse = "|"))]
    cov_gam[,paste0(gam_quads, "2")] <- cov_gam[,gam_quads] ^ 2
  }
  # Put columns in alphabetical order
  cov_gam <- cov_gam %>%
    select(order(colnames(.)))
  # Add interactions if needed
  if (N_GAM_INTERACTS > 0) {
    gam_interacts <- matrix(NA, nrow = 1, ncol = 2)
    for (i in 1:N_GAM_INTERACTS) {
      gam_interacts[i,] <- get(paste0("GAM_INT", i))
      gam_interacts[i,1] <- ifelse(gam_interacts[i,1] %in% covs_cont,
                                   paste0(gam_interacts[i,1], "_z"),
                                   gam_interacts[i,1])
      gam_interacts[i,2] <- ifelse(gam_interacts[i,2] %in% covs_cont,
                                   paste0(gam_interacts[i,2], "_z"),
                                   gam_interacts[i,2])    
      cov_gam[,paste0(gam_interacts[i,], collapse = "_")] <- 
        cov_gam[,gam_interacts[i,1]] * cov_gam[,gam_interacts[i,2]]
    }
  }
  # Convert to a matrix
  cov_gam <- as.matrix(cov_gam)
}

# Calculate the number of covariates for each parameter  
if (exists("cov_psi")) {n_cov_psi <- ncol(cov_psi)}
if (exists("cov_p")) {n_cov_p <- ncol(cov_p)}
if (exists("cov_eps")) {n_cov_eps <- ncol(cov_eps)}
if (exists("cov_gam")) {n_cov_gam <- ncol(cov_gam)}

#------------------------------------------------------------------------------#
# Package things up for JAGS
#------------------------------------------------------------------------------#

# z (latent occupancy for each site & season) is what we're interested in
# In JAGS, z will be stored in a matrix (n_sites[i] * n_seasons[t])

# In order to create initial values for z, we'll need to summarize detections 
# over occasions at each site in each season. To do this:
# Create an array with detection data (row = site, col = occ, slice = season)
# Then use an apply function, summarizing over columns

n_sites <- max(surveys$site_index)
n_seasons <- max(surveys$season_index)
max_n_occasions <- max(surveys$occ_index)
season_occ <- data.frame(season_index = 1:n_seasons)
season_occ$yr <- surveys$yr[match(season_occ$season_index, surveys$season_index)]
for (i in 1:n_seasons) {
  season_occ$n_occasions[i] <- 
    ifelse(sum(surveys$season_index == i) == 0, NA, 
           max(surveys$occ_index[surveys$season_index == i]))
}

# Create an array with detection data (row = site, col = occ, slice = season)
y_array <- array(NA, dim = c(n_sites, max_n_occasions, n_seasons))
for (i in 1:n_seasons) {
  if (is.na(season_occ$n_occasions[i])) {next}
  y_array[,1:season_occ$n_occasions[i],i] <-
    as.matrix(dh_df[,grep(season_occ$yr[i], colnames(dh_df))])
}

# Function to create a matrix with information about known latent states, z[i,t]
# JAGS won't try to estimate z when site is known to be occupied
known_state_occ <- function(y_array){
  state <- apply(y_array, 
                 c(1,3), 
                 function(x) ifelse(sum(is.na(x)) == length(x), 
                                    NA, max(x, na.rm=T)))
  state[state==0] <- NA
  return(state)
}
# check:
# known_state_occ(y_array)

# Function to create initial values for unknown latent states, z[i, t]
inits_state_occ <- function(y_array){
  state <- apply(y_array, 
                 c(1,3), 
                 function(x) ifelse(sum(is.na(x)) == length(x), 
                                    NA, 
                                    max(x, na.rm=T)))
  # Initial value of 1 whenever occupancy state is unknown
  state[state==1] <- 2
  state[is.na(state) | state==0] <- 1
  state[state==2] <- NA
  return(state)
}  
# check:
# inits_state_occ(y_array)

# Bundle data for JAGS that will be needed in any model, regardless of covariate 
# structure
jags_data <- list(y = surveys$det,
                  n_sites = n_sites,
                  n_seasons = n_seasons,
                  n_obs = nrow(surveys),
                  n_sitetrans = nrow(sitetrans),
                  site = surveys$site_index,
                  season = surveys$season_index,
                  site_ec = sitetrans$site_index,
                  trans_ec = sitetrans$trans_index,
                  z = known_state_occ(y_array))

# Vector of parameters to monitor (needed for all models)
params <- c("mean_psi", "beta_psi0", "mean_p", "beta_p0",
            "mean_eps", "beta_eps0", "mean_gam", "beta_gam0", "PAO")

# List of initial values (needed for all models)
inits_list <- list(mean_psi = runif(1, 0, 1),
                   mean_p = runif(1, 0, 1),
                   mean_eps = runif(1, 0, 1),
                   mean_gam = runif(1, 0, 1),
                   z = inits_state_occ(y_array))

# Add in data, parameters, inits associated with covariates
if (exists("cov_psi")) {
  jags_data <- c(jags_data,
                 list(cov_psi = cov_psi),
                 list(n_cov_psi = n_cov_psi))
  params <- c(params, "beta_psi")
  inits_list <- c(inits_list,
                  list(beta_psi = runif(n_cov_psi, -2, 2)))
}
if (exists("cov_p")) {
  jags_data <- c(jags_data,
                 list(cov_p = cov_p),
                 list(n_cov_p = n_cov_p))
  params <- c(params, "beta_p")
  inits_list <- c(inits_list,
                  list(beta_p = runif(n_cov_p, -2, 2)))
}
if (exists("cov_eps")) {
  jags_data <- c(jags_data,
                 list(cov_eps = cov_eps),
                 list(n_cov_eps = n_cov_eps))
  params <- c(params, "beta_eps")
  inits_list <- c(inits_list,
                  list(beta_eps = runif(n_cov_eps, -2, 2)))
}
if (exists("cov_gam")) {
  jags_data <- c(jags_data,
                 list(cov_gam = cov_gam),
                 list(n_cov_gam = n_cov_gam))
  params <- c(params, "beta_gam")
  inits_list <- c(inits_list,
                  list(beta_gam = runif(n_cov_gam, -2, 2)))
}

inits <- function(){inits_list}

# Identify correct JAGS model
jags_model <- "JAGS/JAGS_MS"
if (exists("cov_psi")) {jags_model <- paste0(jags_model, "_", "psi")}
if (exists("cov_p")) {jags_model <- paste0(jags_model, "_", "p")}
if (exists("cov_eps")) {jags_model <- paste0(jags_model, "_", "eps")}
if (exists("cov_gam")) {jags_model <- paste0(jags_model, "_", "gam")}
jags_model <- paste0(jags_model, ".txt")

#------------------------------------------------------------------------------#
# Run model in JAGS
#------------------------------------------------------------------------------#

nc <- 3      # Number of chains
na <- 5000   # Number of iterations to run in the adaptive phase
nb <- 10000  # Number of iterations to discard (burn-in)
ni <- 30000  # Number of iterations per chain (including burn-in)
nt <- 20     # Thinning rate
# Note: with this number of iterations, the model can take many minutes to run

out <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags_model,
            n.chains = nc,
            n.adapt = na,
            n.burnin = nb,
            n.iter = ni,
            n.thin = nt,
            parallel = TRUE)

# print(out)

# Trace and density plots
# MCMCtrace(out,pdf = FALSE)
