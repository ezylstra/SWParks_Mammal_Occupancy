################################################################################
# SODN -- Camera trap data, 2016-2022
# Delineate sampling occasions

# ER Zylstra
# Updated: 2022-10-20
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

# rm(list = ls())

source("src/photo-data/format-mammal-data.R")

  # dat = information about each photo (date, time, species, location)
  # events = information about each camera deployment (dates, location, duration)
  # event_mat = camera location x day matrix with 1/0 indicating whether camera
  #             was deployed or not
  # locs = information about each camera location (park, lat/long, name)
  # species = table with species observed (species code, common name, # of obs)

# Create a dataframe with information about each day of the study
# Day number (daynum): day 1 = 01 Jan 2016)
days_df <- data.frame(daynum = 1:max(events$r_day), 
                      date = seq(as.Date("2016-01-01"), 
                                 max(events$r_date), 
                                 by = 1))
days_df$yr <- year(days_df$date)

# Set the length of sampling occasions, in days
occ_length <- 7
  # May want to evaluate if 7 days is the best choice (Iannarilli et al. 2019?)

# Set the maximum number of sampling occasions in a year (because batteries 
# will start dying)
occ_max <- 6

# Identify the minimum proportion of cameras in a park that need to be 
# operational for that date to be included in a sampling occasion
threshold <- 0.60

# List of parks to include
parks <- c("CHIR", "ORPI", "SAGW")

# Range of photo observations for each park
park_photos <- dat %>%
  filter(Park %in% parks) %>%
  group_by(Park) %>%
  summarize(first_yr = min(yr),
            last_yr = max(yr)) %>%
  data.frame

n <- 0

for (park in parks) {

  n <- n + 1

  # Extract events data for selected park
  event_park <- events %>%
    filter(Park == park) %>%
    select(StdLocName, d_date, r_date, d_yr, r_yr, duration, d_day, r_day)
    
  # Limit events to only those years when we have photo data
  yr_min <- park_photos$first_yr[park_photos$Park == park]
  yr_max <- park_photos$last_yr[park_photos$Park == park]
  event_park <- event_park %>%
    filter(d_yr %in% yr_min:yr_max)
  
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
  # filter(days_park, yr == 2017 & prop_deploy > 0)
  
  # Identify "streaks" of days when the threshold is met
  days_park1 <- days_park %>%
    filter(at_thresh == 1)
  days_park1$streak <- 1
  for (i in 2:nrow(days_park1)) {
    days_park1$streak[i] <- ifelse(days_park1$daynum[i] - days_park1$daynum[i-1] == 1,
                                   days_park1$streak[i-1],
                                   days_park1$streak[i-1] + 1)
  }
  
  # Starting a new "streak" at CHIR in 2019, when all cameras were re-deployed 
  # in summer (between 17 and 21 June)
  if (park == "CHIR") {
    days_park1$streak <- ifelse(days_park1$date >= "2019-06-21",
                                days_park1$streak + 1,
                                days_park1$streak)
  }

  # Identify the first day in each streak that the threshold was met
  firstday <- days_park1 %>%
    group_by(streak) %>%
    summarize(day1 = min(daynum),
              yr = year(date[daynum == min(daynum)])) %>%
    data.frame
  
  # Generate occasion number (within a streak)
  days_park$streak <- days_park1$streak[match(days_park$daynum, days_park1$daynum)]
  days_park$occasion <- NA
  for (i in days_park1$daynum) {
    days_park$occasion[i] <- 
      ceiling((i - firstday$day1[firstday$streak == days_park$streak[i]] + 1) / occ_length)
  }
  
  # Create a dataframe with information about each occasion
  occasions <- days_park %>%
    filter(!is.na(occasion)) %>%
    group_by(Park, streak, occasion) %>%
    summarize(duration = length(occasion),
              start = min(date),
              end = max(date)) %>%
    as.data.frame  
  
  # Exclude any occasions shorter than occ_length, and then retain a maximum
  # of occ_max occasions per streak/year
  occasions$full_duration <- 1*(occasions$duration == occ_length)
  occasions_max <- occasions %>%
    group_by(streak) %>%
    summarize(max_occ = max(occasion[full_duration == 1]),
              yr = year(min(start))) %>%
    data.frame
  occasions_max$max_keep <- ifelse(occasions_max$max_occ > occ_max, 
                                   occ_max, 
                                   occasions_max$max_occ)
  occasions <- left_join(occasions, 
                         select(occasions_max, c(streak, max_keep, yr)), 
                         by = "streak")
  occasions$keep <- ifelse(occasions$occasion > occasions$max_keep, 0, 1)

  # Remove occasions that won't be used in analysis from the dataframe
  # (Assuming only one streak starts per year, we can use yr instead of streak)
  occasions <- occasions %>% 
    filter(keep == 1) %>%
    select(-c(streak, full_duration, max_keep, keep))

  # Add occasion ID and convert occasion start/end dates to day numbers
  occasions <- occasions %>%
    mutate(yr_occ = paste0(occasions$yr, "_", occasions$occasion),
           start_day = as.numeric(date(start)) - as.numeric(as.Date("2015-12-31")),
           end_day = as.numeric(date(end)) - as.numeric(as.Date("2015-12-31")))

  # Append information about sampling occasions to occasions_allparks
  if (n == 1) {
    occasions_allparks <- occasions
  } else {
    occasions_allparks <- rbind(occasions_allparks, occasions)
  }
}

# Export occasions dataframe as csv
# write.csv(occasions_allparks, 
#           file = "data/occasions/occasions-all-parks.csv",
#           row.names = FALSE)

