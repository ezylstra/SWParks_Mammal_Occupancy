################################################################################
# Run a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-10-19
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(jagsUI)
library(MCMCvis)

#------------------------------------------------------------------------------#
# Load all photo, location, events, species data 
#------------------------------------------------------------------------------#

source("src/photo-data/format-mammal-data.R")

# dat = information about each photo (date, time, species, location)
# events = information about each camera deployment (dates, location, duration)
# event_mat = camera location x day matrix with 1/0 indicating whether camera
#             was deployed or not
# locs = information about each camera location (park, lat/long, name)
# species = table with species observed (species code, common name, # of obs)

#------------------------------------------------------------------------------#
# Specify model parameters (objects in all caps)
#------------------------------------------------------------------------------#

# Select park of interest ("CHIR", "ORPI", or "SAGW")
PARK <- "SAGW"

# Select species of interest (select option from species$Species_code)
SPECIES <- "LECA"

# Specify years to include (earliest possible year = 2016 for ORPI, 2017 other parks)
YEARS <- 2017:2022

# Covariates for initial occupancy (psi)
  # Options (all parks): elev, boundary, pois, roads, trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_PSI <- c("pois", "elev", "vegclass")

  # If you want to include quadratic effects for any of the elements in 
  # COVARS_PSI, list them by name. If no quadratic effects, PSI_QUADS <- NA
  PSI_QUADS <- c("elev", "pois")

# Covariates for detection probability (p)
  # Options (all parks): effort, day, camera_new, deploy_exp, elev, boundary, 
  # pois, roads, trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_P <- c("effort", "day", "camera_new", "deploy_exp")

  # If you want to include quadratic effects for any of the elements in 
  # COVARS_P, list them by name. If no quadratic effects, P_QUADS <- NA
  P_QUADS <- c("day")

# Covariates for extinction probability (eps)
  # Options (all parks): monsoon_ppt, elev
  # Set as NA if you don't want any covariates
  COVARS_EPS <- c("elev", "monsoon_ppt")

  # If you want to include quadratic effects for any of the elements in 
  # COVARS_EPS, list them by name. If no quadratic effects, EPS_QUADS <- NA
  EPS_QUADS <- NA  

  # To include interactions between covariates:
  # Specify the number of desired interactions, and create a vector for each 
  # interaction with the names of the two variables to include (note: these 
  # variables must appear in COVARS_EPS)
  N_EPS_INTERACTS <- 1
  EPS_INT1 <- c("elev", "monsoon_ppt")

# Covariates for colonization probability (gam)
  # Options (all parks): monsoon_ppt, elev
  # Set as NA if you don't want any covariates
  COVARS_GAM <- c("elev", "monsoon_ppt")
  
  # If you want to include quadratic effects for any of the elements in 
  # COVARS_GAM, list them by name. If no quadratic effects, GAM_QUADS <- NA
  GAM_QUADS <- NA  
  
  # To include interactions between covariates:
  # Specify the number of desired interactions, and create a vector for each 
  # interaction with the names of the two variables to include (note: these 
  # variables must appear in COVARS_GAM)
  N_GAM_INTERACTS <- 0
  # GAM_INT1 <- c("monsoon_ppt", "elev")

#------------------------------------------------------------------------------#  
# Run multi-season model and save results
#------------------------------------------------------------------------------#
  
source("src/multi-season-models/MSoccupancy-generic.R")

# Create filename to store results
# For now, using park, species, and date in filename but could add more 
# descriptors if needed (e.g., something about covariates):
date <- Sys.Date()
model_filename <- paste0("output/models/",
                         tolower(PARK), "-",
                         tolower(SPECIES), "-MS-",
                         date,
                         ".Rdata")

# List of all objects we want to save if they exist
obj_to_save <- c("out", "surveys", "spatial_covs", "sitetrans", "occasions",
                 "n_cov_psi", "n_cov_p", "n_cov_eps", "n_cov_gam",
                 "COVARS_PSI", "COVARS_P", "COVARS_EPS", "COVARS_GAM", 
                 ls()[str_detect(ls(), "INT")])
# Remove items from the list if they don't exist
obj_to_save <- obj_to_save[sapply(obj_to_save, exists)]

# Save to file
save(list = obj_to_save, file = model_filename)
