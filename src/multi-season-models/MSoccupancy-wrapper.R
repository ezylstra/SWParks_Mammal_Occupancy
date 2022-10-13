################################################################################
# Run a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-10-13
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

# Park of interest ("CHIR", "ORPI", or "SAGW")
PARK <- "SAGW"

# Species of interest (select option from species$Species_code)
SPECIES <- "LECA"

# Specify years to include (earliest possible year = 2016)
YEARS <- 2017:2022

# Covariates for initial occupancy (psi)
  # Options (all parks): elev, boundary, pois, roads, trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_PSI <- c("pois", "elev", "vegclass")

  # Indicate (by number) if you want to include quadratic effects for any of the
  # elements in covariates_psi (NA if none; c(1,3) would include linear and 
  # quadratic effects of the 1st and 3rd covariates in the covariates_psi vector
  # above)
  PSI_QUADS <- c(2,1)

# Covariates for detection probability (p)
  # Options (all parks): effort, day, camera_new, elev, boundary, pois, roads, 
  # trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_P <- c("effort", "day", "camera_new")

  # Indicate (by number) if you want to include quadratic effects for any of the
  # elements in covariates_p
  P_QUADS <- c(2)

# Covariates for extinction probability (eps)
  # Options (all parks): monsoon_ppt, elev
  # Set as NA if you don't want any covariates
  COVARS_EPS <- c("elev", "monsoon_ppt")

  # Indicate (by number) if you want to include quadratic effects for any of the
  # elements in covariates_eps
  EPS_QUADS <- NA  

  # To include interactions between covariates:
  # Specify the number of desired interactions, and create a vector with the 
  # names of the two variables to include (note: these variables must appear in 
  # covariates_eps)
  N_EPS_INTERACTS <- 1
  EPS_INT1 <- c("elev", "monsoon_ppt")

# Covariates for colonization probability (gam)
  # Options (all parks): monsoon_ppt, elev
  # Set as NA if you don't want any covariates
  COVARS_GAM <- c("elev", "monsoon_ppt")
  
  # Indicate (by number) if you want to include quadratic effects for any of the
  # elements in covariates_gam
  GAM_QUADS <- NA  
  
  # To include interactions between covariates:
  # Specify the number of desired interactions, and create a vector with the 
  # names of the two variables to include (note: these variables must appear in 
  # covariates_gam)
  N_GAM_INTERACTS <- 0
  # GAM_INT1 <- c("monsoon_ppt", "elev")

#------------------------------------------------------------------------------#  
# Run multi-season model and save results
#------------------------------------------------------------------------------#
  
source("src/multi-season-models/MSoccupancy-generic.R")
  
# Save JAGS model and other objects to file in output/models
date <- Sys.Date()
  
# For now, using park, species, and date to identify the model, but could
# add more descriptors to the filename:
model_filename <- paste0("output/models/",
                         tolower(PARK), "-",
                         tolower(SPECIES), "-MS-",
                         date,
                         ".Rdata")

obj_to_save <- c("out", "surveys", "spatial_covs", "sitetrans")
if (!all(is.na(COVARS_PSI))) {obj_to_save <- c(obj_to_save, "cov_psi")}
if (!all(is.na(COVARS_P))) {obj_to_save <- c(obj_to_save, "cov_p")}
if (!all(is.na(COVARS_EPS))) {obj_to_save <- c(obj_to_save, "cov_eps")}
if (!all(is.na(COVARS_GAM))) {obj_to_save <- c(obj_to_save, "cov_gam")}
save(list = obj_to_save, file = model_filename)

