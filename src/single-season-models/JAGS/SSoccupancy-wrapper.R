################################################################################
# Run a single-season occupancy analysis

# ER Zylstra
# Updated 2023-01-18
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

# Prep detection and covariate data
source("src/single-season-models/JAGS/SSoccupancy-prep-data.R")

#------------------------------------------------------------------------------#
# Evaluate potential covariates
#------------------------------------------------------------------------------#

# Brief descriptions of each covariate:
  # boundary = distance to nearest park boundary (m)
  # east = eastness (aspect) Ranges from 1 (east-facing) to -1 (west-facing)
  # elev = elevation (m)
  # north = northness (aspect) Ranges from 1 (north-facing) to -1 (south-facing)
  # pois = distance to nearest point of interest (m; incl buildings, lots)
  # roads = distance to nearest road (m; includes local roads, vehicular trails)
  # slope = slope (degrees)
  # trail = distance to nearest trail (m)
  # wash [SAGW only] = distance to nearest wash (m) 
  # vegclass [SAGW only] = vegetation class where 1 = low gradient desert with 
    # high-cover mixed cactus and 2-15% tree cover; 2 = Low hillslope and 
    # mountain foothills, rocky, often north facing, cooler, wetter; 
    # 3 = Medium-high gradient, contrasting topography (hilly), often Jojoba 
    # dominant; 4 = developed (no cameras located in this vegclass)
  # burn_severity_2011 [CHIR only] = severity of 2011 burn (integer values, 0:4 
    # with 0 = unburned to 4 = high burn severity)
  # effort = proportion of days during sampling occasion that camera was running
  # day = day of the year (1:366)
  # deploy_exp = experience of crew deploying camera (0 = all novices; 1 = at
    # least one experience person present; 2 = at least one expert present)

# Look at pairwise correlations among continuous covariates. May want to 
# avoid including a pair of covariates in a model for occupancy if they are 
# highly correlated: |r| > 0.6 or 0.7.

message("Identify pairs of covariates are highly correlated. 
(May want to avoid including pairs with |r| > 0.6 or 0.7 in the same model)")
cor_df %>% arrange(desc(abs(corr))) %>% filter(abs(corr) >= 0.6)

#------------------------------------------------------------------------------#
# Specify covariates (objects in all caps)
#------------------------------------------------------------------------------#

# Covariates for occupancy (psi)
  # Options (all parks): elev, boundary, pois, roads, trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_PSI <- c("pois", "elev", "vegclass")

  # If you want to include quadratic effects for any of the elements in 
  # COVARS_PSI, list them by name. If no quadratic effects, PSI_QUADS <- NA
  PSI_QUADS <- c("elev")

# Covariates for detection probability (p)
  # Options (all parks): effort, day, deploy_exp, elev, boundary, pois, roads, 
  # trail, east, north, slope
  # Additional options (SAGW): wash, vegclass
  # Additional options (CHIR): burn_severity_2011
  # Set as NA if you don't want any covariates
  COVARS_P <- c("effort", "day", "deploy_exp")

  # If you want to include quadratic effects for any of the elements in 
  # COVARS_P, list them by name. If no quadratic effects, P_QUADS <- NA
  P_QUADS <- c("day")

  source("src/single-season-models/JAGS/SSoccupancy-covariate-check.R")

#------------------------------------------------------------------------------#  
# Run a single-season model and save results
#------------------------------------------------------------------------------# 

source("src/single-season-models/JAGS/SSoccupancy-run-model.R")  

# Create filename to store results
# For now, using park, species, yr, and date in filename but could change 
# descriptors if needed (e.g., something about covariates):
date <- str_remove_all(Sys.Date(), "-")
model_filename <- paste0("output/models-JAGS/",
                         tolower(PARK), "-",
                         tolower(SPECIES), "-",
                         YEAR, "-", date, ".Rdata")

# List of all objects we want to save if they exist
obj_to_save <- c("out", "surveys", "spatial_covs", "occasions",
                 "cov_psi", "cov_p", "COVARS_PSI", "COVARS_P", 
                 "model_description")
# Remove items from the list if they don't exist
obj_to_save <- obj_to_save[sapply(obj_to_save, exists)]

# Save to file
save(list = obj_to_save, file = model_filename)
