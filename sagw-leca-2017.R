################################################################################
# Single-season occupancy analysis
# SAGW, LECA (black-tailed jackrabbit), 2017

# ER Zylstra
# 2022-06-15
################################################################################

# Note: Will probably want to create a template script for single-season analyses
# This will serve as a first case study

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

# rm(list = ls())

# Load photo observations, camera location info, event matrix, and species list
source("format-mammal-data.R")

# Load sampling occasion data
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

# Identify park, species, and year of interest
park <- "SAGW"
species <- "LECA"
year <- 2017

# Use day number (day 1 = 01 Jan 2016) for column names in event_mat 
colnames(event_mat) <- 1:ncol(event_mat)

# Extract sampling occasion info for selected park and year
occasions <- occasions %>%
  filter(Park == park & yr == year)

# Select photo observations for park, species and year
# Retain a maximum of one observation per day at each location
obs <- dat %>% 
  filter(Park == park & Species_code == species & yr == year) %>%
  select(StdLocName, obsdate, yr, o_day) %>%
  arrange(StdLocName, obsdate) %>%
  distinct

# Extract information about camera locations in selected park
locs_park <- locs %>%
  filter(UnitCode == park) %>%
  select(StdLocName, POINT_X, POINT_Y) %>%
  rename(loc = StdLocName, long = POINT_X, lat = POINT_Y)

# Extract rows from events matrix that correspond to locations in selected park
event_mat <- event_mat[rownames(event_mat) %in% locs_park$loc,]

# Convert occasion start/end dates to day numbers
occasions$start_day <- 
  as.numeric(date(occasions$start)) - as.numeric(as.Date("2015-12-31"))
occasions$end_day <- 
  as.numeric(date(occasions$end)) - as.numeric(as.Date("2015-12-31"))

# Create a list of days included in sampling occasions
occ_days <- NULL
for (i in 1:nrow(occasions)) {
  occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
}

# Extract columns from events matrix that correspond to sampling occasions
event_mat <- event_mat[,colnames(event_mat) %in% occ_days]

# Convert the events matrix into detection histories (ddh = daily dh)
ddh <- event_mat

# Change 0s to NA (NA = camera wasn't deployed)
ddh[ddh == 0] <- NA

# Change 1s to 0 (0 indicates that the species wasn't detected)
ddh[ddh == 1] <- 0

# Replace 0s with 1s when the species was detected
for (i in 1:nrow(obs)) {
  ddh[rownames(ddh) == obs$StdLocName[i], 
      colnames(ddh) == as.character(obs$o_day[i])] <- 1
}
# checks:
sum(ddh == 1, na.rm = TRUE)
sum(obs$o_day %in% occ_days)

# Create a function to aggregate daily detection data during each occasion
  # NA if the camera was not operational throughout entire occasion (all values = NA)
  # 1 if species was detected one or more times (regardless if there are NAs)
  # 0 if species was never detected
  paNA <- function(x) {
    if (sum(is.na(x)) == length(x)) {NA} else 
      if (sum(x, na.rm = TRUE) == 0) {0} else {1} 
  }
  
# Create a function to calculate the proportion of a sampling period a camera
# was operational
  propNA <- function(x) {
    (occasions$duration[1] - sum(is.na(x))) / occasions$duration[1]
  }

# Summarize detection data (dh) and effort during each occasion 
dh <- effort <-  matrix(NA, 
                        nrow = nrow(ddh), 
                        ncol = ncol(ddh) / occasions$duration[1],
                        dimnames = list(rownames(ddh), NULL))

for (i in 1:ncol(dh)) {
  multiday <- ddh[,colnames(ddh) %in% occasions$start_day[i]:occasions$end_day[i]]
  dh[,i] <- apply(multiday, 1, paNA)
  effort[,i] <- apply(multiday, 1, propNA)
}
  # Might need to replace 0 values in effort with NA

#-------------------------------------------------------------------------------#
# Spatial covariates
#-------------------------------------------------------------------------------#
# Will expand this section as more covariates become available

spatial_covs <- locs_park 

# Ensure the order is the same as what's in the detection history matrix
spatial_covs <- spatial_covs[match(rownames(dh), spatial_covs$loc),]

