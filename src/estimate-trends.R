################################################################################
# Estimate trends based on a multi-season dynamic model
# SAGW, LECA (black-tailed jackrabbit)

# ER Zylstra
# Updated 2022-07-07
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(jagsUI)

# Load photo, location, events, species data 
source("src/format-mammal-data.R")

# dat = information about each photo (date, time, species, location)
# events = information about each camera deployment (dates, location, duration)
# event_mat = camera location x day matrix with 1/0 indicating whether camera
#             was deployed or not
# locs = information about each camera location (park, lat/long, name)
# species = table with species observed (species code, common name, # of obs)

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

park <- "SAGW"
species <- "LECA"

# Will eventually need some character string to identify the model with 
# particular covariates, information contained in the name of the model file.
# For now, will just use "MS-test" for model created in sagw-leca-multiseason.R

# Load JAGS model 
model_file <- paste0("output/models/",
                     tolower(park), "-",
                     tolower(species), "-",
                     "MS-test.rds")
jags_model <- readRDS(file = model_file)

# TO DO:
# extract posterior samples for PAO





