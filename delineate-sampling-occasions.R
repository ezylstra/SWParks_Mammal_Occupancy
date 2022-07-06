################################################################################
# SODN -- Camera trap data, 2016-2022
# Delineate sampling occasions

# ER Zylstra
# 2022-06-02
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

# rm(list = ls())

source("format-mammal-data.R")

  # dat = information about each photo (date, time, species, location)
  # events = information about each camera deployment (dates, location, duration)
  # event_mat = camera location x day matrix with 1/0 indicating whether camera
  #             was deployed or not
  # locs = information about each camera location (park, lat/long, name)
  # species = table with species observed (species code, common name, # of obs)

# Create a dataframe with information about each day of the study
days_df <- data.frame(daynum = 1:max(dat$o_day), 
                      date = seq(as.Date("2016-01-01"), max(dat$obsdate), by = 1))
days_df$yr <- year(days_df$date)

# Use day number (day 1 = 01 Jan 2016) for column names in event_mat 
colnames(event_mat) <- days_df[,1]

# Set the length of sampling occasions, in days
occ_length <- 7
  # May want to evaluate if 7 days is the best choice (Iannarilli et al. 2019?)

# Set the minimum proportion of cameras in a park that need to be operational
# for that date to be included in a sampling occasion
threshold <- 0.60

# For now, just delineating occasions at SAGW.
# Will eventually edit the script to loop through all parks.
n <- 0

for (park in c("SAGW")) {
  n <- n + 1

  # Extract events data for selected park
  event_park <- events %>%
    filter(Park == park) %>%
    select(StdLocName, d_date, r_date, d_yr, r_yr, duration, d_day, r_day)
  
  # Extract rows from events matrix that correspond to locations in selected park
  event_mat_park <- event_mat[rownames(event_mat) %in% event_park$StdLocName,]
  
  # Calculate the number of cameras that are deployed each day
  days_park <- days_df
  days_park$Park <- park
  days_park$n_cameras <- colSums(event_mat_park)
  
  # Calculate proportion of cameras in selected park that are deployed each day
  days_park$prop_deploy <- days_park$n_cameras / nrow(event_mat_park)
  
  # Identify those dates when the proportion of cameras deployed >= threshold
  days_park$at_thresh <- 1*(days_park$prop_deploy >= threshold)
  
  # Look at consecutive days with sufficient number of cameras deployed
  # data.frame(unclass(rle(days_park$at_thresh)))
  # filter(days_df, yr == 2017 & prop_deploy > 0)
  
  # Identify the first day in each year that the threshold was met
  firstday_yr <- days_park %>%
    filter(at_thresh == 1) %>%
    group_by(yr) %>%
    summarize(day1 = min(daynum)) %>%
    as.data.frame
  
  # Generate occasion number (within a year)
  days_park$occasion <- NA
  for (i in days_park$daynum[days_park$at_thresh == 1]) {
    days_park$occasion[i] <- 
      ceiling((i - firstday_yr$day1[firstday_yr$yr == days_df$yr[i]] + 1) / occ_length)
  }
  
  # Create a dataframe with information about each occasion
  occasions <- days_park %>%
    filter(!is.na(occasion)) %>%
    group_by(Park, yr, occasion) %>%
    summarize(duration = length(occasion),
              start = min(date),
              end = max(date)) %>%
    as.data.frame
  
  # Exclude any occasions shorter than occ_length
  occasions$keep <- 1*(occasions$duration == occ_length) 
  days_park <- left_join(days_park, select(occasions, c(yr, occasion, keep)))
  days_park$keep[is.na(days_df$keep)] <- 0
  # occasions
  # days_park[760:800,]
  
  # Remove occasions that won't be used in analysis from the dataframe
  occasions <- occasions %>% 
    filter(keep == 1) %>%
    select(-keep)

  # Append information about sampling occasions to occasions_allparks
  if (n == 1) {
    occasions_allparks <- occasions
  } else {
    occasions_allparks <- rbind(occasions_allparks, occasions)
  }
  
  # Export days_park dataframe as csv?
  # write.csv(days_park,
  #           file = paste0(park, "_days.csv"),
  #           row.names = FALSE)
  
}

# Export occasions dataframe as csv
write.csv(occasions_allparks, 
          file = "data/occasions/occasions-all-parks.csv",
          row.names = FALSE)

