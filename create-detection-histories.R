################################################################################
# SODN -- Camera trap data, 2016-2022
# Creating encounter histories

# ER Zylstra
# 2022-06-02
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

# rm(list = ls())

#-------------------------------------------------------------------------------#
# Run script to import, format data
#-------------------------------------------------------------------------------#

source("format-mammal-data.R")

#-------------------------------------------------------------------------------#
# Subset data by species and park
#-------------------------------------------------------------------------------#

spp <- "LECA"

park <- "SAGW"

datsub <- dat %>% 
  filter(Species_code == spp & Park == park) %>%
  select(LocationID, StdLocName, obsdate, yr, mon, yday, o_day, obstime, POINT_X, POINT_Y)

eventsub <- events %>%
  filter(Park == park) %>%
  select(StdLocName, d_date, r_date, d_yr, r_yr, duration, d_day, r_day)

#-------------------------------------------------------------------------------#
# 
#-------------------------------------------------------------------------------#

# To do:

# Extract rows from events_mat for selected park
# Trim yearly deployments
  # 28 day minimum [4, 7-day occasions] so most/all cameras deployed over entire period
  # Select only those columns that fall within trimmed deployments.
  # Change values to NA when camera wasn't operational.
# Create daily detection histories
  # In matrix above, change values to 0 when camera was operational.
  # Change 0s to 1s when species was detected
# Create detection histories for multi-day "occasions"
  # Can have 1 or 0 even if camera wasn't operational the entire time
# Create survey covariate that represents "effort" (prop of days camera operational)

# Create single-season models in JAGS
# Think about how we want to specify multi-season models


