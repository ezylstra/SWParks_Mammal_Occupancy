################################################################################
# SODN -- Camera trap data, 2016-2022
# Creating detection histories

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
# Note: will likely want to make this more general (eg, pick > 1 spp)

spp <- "LECA"

park <- "SAGW"

# Retain a maximum of one observation per day at each location
datsub <- dat %>% 
  filter(Species_code == spp & Park == park) %>%
  select(LocationID, StdLocName, obsdate, yr, mon, yday, o_day, POINT_X, POINT_Y) %>%
  distinct

#-------------------------------------------------------------------------------#
# Define sampling periods, occasions for selected park
#-------------------------------------------------------------------------------#

# Length of secondary occasions, in days
occ_length <- 7
  # May want to evaluate whether 7 days is the best choice (Iannarilli et al. 2019?)

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
ddh <- event_mat[rownames(event_mat) %in% eventsub$StdLocName,]

# Calculate the number of cameras that are deployed each day
days_df$n_cameras <- colSums(ddh)
  
# Calculate proportion of cameras in selected park that are deployed each day
days_df$prop_deploy <- days_df$n_cameras / nrow(ddh)

# Identify those dates when > threshold proportion of cameras deployed
threshold <- 0.60
days_df$abovethresh <- 1*(days_df$prop_deploy > threshold)
  
# Look at consecutive days with sufficient number of cameras deployed
data.frame(unclass(rle(days_df$abovethresh)))
  # filter(days_df, yr == 2017 & prop_deploy > 0)

# Identify the first day in each year that the threshold was met
firstday_yr <- days_df %>%
  filter(abovethresh == 1) %>%
  group_by(yr) %>%
  summarize(day1 = min(daynum)) %>%
  as.data.frame

# Generate occasion number (within a year)
days_df$occasion <- NA
for (i in days_df$daynum[days_df$abovethresh == 1]) {
  days_df$occasion[i] <- 
    ceiling((i - firstday_yr$day1[firstday_yr$yr == days_df$yr[i]] + 1) / occ_length)
}

# Exclude any occasions shorter than occ_length
occasions <- days_df %>%
  filter(!is.na(occasion)) %>%
  group_by(yr, occasion) %>%
  summarize(duration = length(occasion),
            start = min(date),
            end = max(date)) %>%
  as.data.frame
occasions$keep <- 1*(occasions$duration == occ_length) 
days_df <- left_join(days_df, select(occasions, c(yr, occasion, keep)))
days_df$keep[is.na(days_df$keep)] <- 0
  # occasions
  # days_df[760:800,]

# Extract just those columns from dh that correspond to sampling occasions
ddh <- ddh[, colnames(ddh) %in% days_df$daynum[days_df$keep == 1]]

#-------------------------------------------------------------------------------#
# Create detection histories 
#-------------------------------------------------------------------------------#

# Change 0s in dh to NA (NA = camera wasn't operational)
ddh[ddh == 0] <- NA

# Change 1s in dh to 0 (0 will indicate that the species wasn't detected)
ddh[ddh == 1] <- 0

# Replace 0s with 1s when the species was detected
for (i in 1:nrow(datsub)) {
  ddh[rownames(ddh) == datsub$StdLocName[i], 
     colnames(ddh) == as.character(datsub$o_day[i])] <- 1
}
# checks:
table(ddh)
days_df[days_df$daynum %in% datsub$o_day[!datsub$o_day %in% colnames(ddh)],]
  # Detections that aren't included in detection histories are either:
    # on days when a low proportion of cameras were deployed (< threshold)
    # or, during occasions that were < occ_length

# Create a function to aggregate daily detection data during each occasion
  # NA if the camera was not operation during entire occasion (all values = NA)
  # 1 if species was detected one or more times (regardless if there are NAs)
  # 0 if species was never detected
  paNA <- function(x) {
    if (sum(is.na(x)) == length(x)) {NA} else 
      if (sum(x, na.rm = TRUE) == 0) {0} else {1} 
  }

# Summarize detection data (dh) and effort (eff) during each occasion 
occ_starts <- seq(1, ncol(ddh), by = occ_length)
occ_ends <- occ_starts + occ_length - 1
dh <- eff <-  matrix(NA, 
                     nrow = nrow(ddh), 
                     ncol = ncol(ddh) / occ_length,
                     dimnames = list(rownames(ddh), NULL))
for (i in 1:ncol(dh)) {
  multiday <- ddh[,occ_starts[i]:occ_ends[i]]
  dh[,i] <- apply(multiday, 1, paNA)
  eff[,i] <- apply(multiday, 1, function(x) {(occ_length - sum(is.na(x))) / occ_length})
}

#-------------------------------------------------------------------------------#
# Temporal covariates
#-------------------------------------------------------------------------------#
# Will expand this section as more covariates become available

temporal_covs <- occasions %>%
  filter(keep == 1) %>%
  select(-c(duration, keep)) %>%
  rowwise %>%
  mutate(mid = mean.Date(c(start,end)),
         y2018 = 1*(yr == 2018),
         y2020 = 1*(yr == 2020),
         y2021 = 1*(yr == 2021),
         y2022 = 1*(yr == 2022),
         trend = yr - 2017) %>%
  as.data.frame

#-------------------------------------------------------------------------------#
# Spatial covariates
#-------------------------------------------------------------------------------#
# Will expand this section as more covariates become available

spatial_covs <- locs %>%
  filter(UnitCode == park) %>%
  select(UnitCode, StdLocName, POINT_X, POINT_Y) %>%
  rename(park = UnitCode, loc = StdLocName, long = POINT_X, lat = POINT_Y)

# Ensure the order is the same as what's in the detection history matrix
spatial_covs <- spatial_covs[match(rownames(ddh), spatial_covs$loc),]

