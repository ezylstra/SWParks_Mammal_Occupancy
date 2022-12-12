################################################################################
# Run a multi-season occupancy analysis

# ER Zylstra
# Updated 2022-12-12
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

# Select year of interest
YEAR <- 2020

# Covariates for occupancy (psi)
  # Options (all parks): elev, boundary, pois, roads, trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_PSI <- c("pois", "elev", "slope")

  # If you want to include quadratic effects for any of the elements in 
  # COVARS_PSI, list them by name. If no quadratic effects, PSI_QUADS <- NA
  PSI_QUADS <- c("elev", "slope")

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

  source("src/single-season-models/SSoccupancy-covariate-check.R")

#------------------------------------------------------------------------------#  
# Run single-season model and save results
#------------------------------------------------------------------------------# 

source("src/single-season-models/SSoccupancy-generic.R")

# Create filename to store results
# For now, using park, species, yr, and date in filename but could change 
# descriptors if needed (e.g., something about covariates):
date <- str_remove_all(Sys.Date(), "-")
model_filename <- paste0("output/models/",
                         tolower(PARK), "-",
                         tolower(SPECIES), "-",
                         YEAR, "-", date, ".Rdata")

# List of all objects we want to save if they exist
obj_to_save <- c("out", "surveys", "spatial_covs", "occasions",
                 "cov_psi", "cov_p", "COVARS_PSI", "COVARS_P")
# Remove items from the list if they don't exist
obj_to_save <- obj_to_save[sapply(obj_to_save, exists)]

# Save to file
save(list = obj_to_save, file = model_filename)
