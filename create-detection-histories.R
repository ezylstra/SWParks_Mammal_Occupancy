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
# Subset photo data by species and park
#-------------------------------------------------------------------------------#
# Note: will eventually want to generalize this, so we could pick > 1 spp

spp <- "LECA"

park <- "SAGW"

datsub <- dat %>% 
  filter(Species_code == spp & Park == park) %>%
  select(LocationID, StdLocName, obsdate, yr, mon, yday, o_day, obstime, POINT_X, POINT_Y)

#-------------------------------------------------------------------------------#
# Define sampling periods, occasions for selected park
#-------------------------------------------------------------------------------#

# Length of secondary occasions (within a year), in days
occ_length <- 7

# Add colnames to event_mat 
days_df <- data.frame(daynum = 1:max(dat$o_day), 
                      date = seq(as.Date("2016-01-01"), max(dat$obsdate), by = 1))
days_df$yr <- year(days_df$date)
colnames(event_mat) <- days_df[,1]

# Subset events data
eventsub <- events %>%
  filter(Park == park) %>%
  select(StdLocName, d_date, r_date, d_yr, r_yr, duration, d_day, r_day)

# Extract rows from events matrix that correspond to locations in selected park
dh <- event_mat[rownames(event_mat) %in% eventsub$StdLocName,]

# Calculate the number of cameras that are deployed each day
days_df$n_cameras <- colSums(dh)
  
# Calculate proportion of cameras in selected park that are deployed each day
days_df$prop_deploy <- days_df$n_cameras / nrow(dh)

# Identify those dates when > threshold proportion of cameras deployed
threshold <- 0.60 # Can adjust this at anyt time
days_df$mostdeployed <- 1*(days_df$prop_deploy > threshold)
  
# Look at consecutive days with sufficient number of cameras deployed
data.frame(unclass(rle(days_df$mostdeployed)))
  # filter(days_df, yr == 2017 & prop_deploy > 0)

# Identify the first day in each year that the threshold was met
firstday_yr <- days_df %>%
  filter(mostdeployed == 1) %>%
  group_by(yr) %>%
  summarize(day1 = min(daynum[mostdeployed == 1])) %>%
  as.data.frame()

# Generate occasion number (within a year)
days_df$occasion <- NA
for (i in days_df$daynum[days_df$mostdeployed == 1]) {
  days_df$occasion[i] <- ceiling((i - firstday_yr$day1[firstday_yr$yr == days_df$yr[i]] + 1) / 7)
}

# Exclude any "occasions" shorter than occ_length
occasions <- days_df %>%
  filter(!is.na(occasion)) %>%
  group_by(yr, occasion) %>%
  summarize(duration = length(occasion)) %>%
  as.data.frame()
occasions$keep <- 1*(occasions$duration == 7) 
days_df <- left_join(days_df, select(occasions, -duration))
days_df$keep[is.na(days_df$keep)] <- 0

# Extract just those columns from dh that correspond to new sampling occasions
dh <- dh[, colnames(dh) %in% days_df$daynum[days_df$keep == 1]]

#-------------------------------------------------------------------------------#
# Create detection histories 
#-------------------------------------------------------------------------------#

# Change 0s in dh to NA (NA = camera wasn't operational)
dh[dh == 0] <- NA

# Change 1s in dh to 0 (0 will indicate that the species wasn't detected)
dh[dh == 1] <- 0

# Replace 0s with 1s when the species was detected


# Next steps:

# Create detection histories for multi-day "occasions"
  # Can have 1 or 0 even if camera wasn't operational the entire time
# Create survey covariate that represents "effort" (prop of days camera operational)


