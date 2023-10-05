################################################################################
# SODN -- Camera trap data, 2016-2022
# Delineate sampling occasions and calculate the number of species detections
# during sampling occasions at each park

# ER Zylstra
# Updated: 2023-10-05
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(sf)
library(terra)

# rm(list = ls())

source("src/photo-data/format-mammal-data.R")

  # dat = information about each photo (date, time, species, location)
  # events = information about each camera deployment (dates, location, duration)
  # event_mat = camera location x day matrix with 1/0 indicating whether camera
  #             was operational or not
  # locs = information about each camera location (park, lat/long, name)
  # species = table with species observed (species code, common name, # of obs)

#------------------------------------------------------------------------------#
# Delineate sampling occasions
#------------------------------------------------------------------------------#

# Create a dataframe with information about each day of the study
# Day number (daynum): day 1 = 01 Jan 2016)
days_df <- data.frame(daynum = 1:max(events$end_day), 
                      date = seq(as.Date("2016-01-01"), 
                                 max(events$active_end), 
                                 by = 1))
days_df$yr <- year(days_df$date)

# Set the length of sampling occasions, in days
occ_length <- 7
  # May want to evaluate if 7 days is the best choice (Iannarilli et al. 2019?)

# Set the maximum number of sampling occasions in a year
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

  # Append information about sampling occasions to occ_all
  if (n == 1) {
    occ_all <- occasions
  } else {
    occ_all <- rbind(occ_all, occasions)
  }
}

# Export occasions dataframe as csv
write.csv(occ_all,
          file = "data/occasions/occasions-all-parks.csv",
          row.names = FALSE)

#------------------------------------------------------------------------------#
# Calculate number of species detections
#------------------------------------------------------------------------------#

n <- 0

for (park in parks) {
  
  n <- n + 1
  
  # Extract sampling occasion info for selected park and years
  park_yrs <- park_photos %>%
    filter(Park == park)
  occasions <- occ_all %>%
    filter(Park == park) %>%
    filter(yr %in% park_yrs$first_yr:park_yrs$last_yr)
  
  # Create a list of days included in sampling occasions
  occ_days <- NULL
  for (i in 1:nrow(occasions)) {
    occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
  }
  
  # Extract photo observations for park and occasion days
  # Retain only one species observation per day at a given location
  obs <- dat %>% 
    filter(Park == park & o_day %in% occ_days) %>%
    select(Park, StdLocName, Species_code, yr, o_day) %>%
    rename(spp = Species_code) %>%
    distinct
  
  # Attach occasion to each row and then distinct again
  for (i in 1:nrow(obs)) {
    obs$yr_occ[i] <- occasions$yr_occ[obs$o_day[i] >= occasions$start_day & 
                                        obs$o_day[i] <= occasions$end_day]
  }
  # Remove replicate observations of species at a location during each occasion
  obs <- obs %>%
    select(-o_day) %>%
    distinct() %>%
    mutate(detect = 1)
  
  # Extract rows from events matrix that correspond to locations in selected park
  locs_park <- locs$StdLocName[locs$UnitCode == park]
  event_mat_park <- event_mat[rownames(event_mat) %in% locs_park,]
  
  # Extract columns from events matrix that correspond to sampling occasions
  event_mat_park <- event_mat_park[,colnames(event_mat_park) %in% occ_days]
  
  # Summarize event data by occasion 
  event_occ <- matrix(NA, 
                      nrow = nrow(event_mat_park), 
                      ncol = nrow(occasions),
                      dimnames = list(rownames(event_mat_park), occasions$yr_occ))
  
  # Create matrix with 1/0 indicating whether camera was operational during that 
  # occasion
  for (i in 1:ncol(event_occ)) {
    multiday <- event_mat_park[,colnames(event_mat_park) %in% 
                                 occasions$start_day[i]:occasions$end_day[i]]
    event_occ[,i] <- apply(multiday, 1, sum)
    event_occ[event_occ > 1] <- 1
  }
  # check: 
  # table(c(event_occ), useNA = "always")
  
  # Convert to long form
  event_occ_long <- event_occ %>%
    as.data.frame() %>%
    mutate(StdLocName = rownames(.)) %>%
    pivot_longer(cols = !last_col(),
                 names_to = "yr_occ",
                 values_to = "obs") %>%
    data.frame()
  
  # Merge observation and species detection information
  detects <- expand.grid(Park = park,
                         StdLocName = sort(locs_park),
                         yr_occ = occasions$yr_occ,
                         spp = unique(obs$spp))
  detects <- left_join(detects, event_occ_long, by = c("StdLocName", "yr_occ"))
  detects <- left_join(detects, obs[,c("StdLocName", "spp", "yr_occ", "detect")],
                       by = c("StdLocName", "spp", "yr_occ"))
  detects$detect[is.na(detects$detect)] <- 0
  detects$yr <- as.numeric(str_sub(detects$yr_occ, 1, 4))
  
  # Summarize by park, yr, and species
  # nobs is number of "observations" (location * occasion when camera operational)
  # ndetects is number of detections of that species (for given year)
  detects_yr <- detects %>%
    group_by(Park, yr, spp) %>%
    summarize(.groups = "keep",
              nobs = sum(obs),
              ndetects = sum(detect)) %>%
    mutate(propdetect = round(ndetects/nobs,2)) %>%
    data.frame()
  
  # Summarize by park and species
  # nobs is number of "observations" (location * occasion when camera operational)
  # ndetects is number of detections of that species
  detects_allyrs <- detects %>%
    group_by(Park, spp) %>%
    summarize(.groups = "keep",
              nobs = sum(obs),
              ndetects = sum(detect)) %>%
    mutate(propdetect = round(ndetects/nobs,2)) %>%
    data.frame()
  detects_allyrs
  
  # Create dataframe with information on species detections at all parks
  if (n == 1) {
    spp_detections_yr <- detects_yr
    spp_detections <- detects_allyrs
  } else {
    spp_detections_yr <- rbind(spp_detections_yr, detects_yr)
    spp_detections <- rbind(spp_detections, detects_allyrs)
  }  
}

# Look at species detection rates by park
spp_detections %>%
  arrange(spp, Park) %>%
  pivot_wider(id_cols = spp, 
              names_from = Park, 
              names_sort = TRUE,
              values_from = propdetect) %>%
  data.frame()

# Export species detections dataframes as csvs
write.csv(spp_detections,
          file = "output/species-detections-bypark.csv",
          row.names = FALSE)
write.csv(spp_detections_yr,
          file = "output/species-detections-byparkyr.csv",
          row.names = FALSE)

#------------------------------------------------------------------------------#
# Calculate number of species detections, 2017-2022 for sharing
#------------------------------------------------------------------------------#

n <- 0

for (park in parks) {
  
  n <- n + 1
  
  # Extract sampling occasion info for selected park and years
  park_yrs <- park_photos %>%
    filter(Park == park)
  occasions <- occ_all %>%
    filter(Park == park) %>%
    filter(yr %in% 2017:2022)
  
  # Create a list of days included in sampling occasions
  occ_days <- NULL
  for (i in 1:nrow(occasions)) {
    occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
  }
  
  # Extract photo observations for park and occasion days
  # Retain only one species observation per day at a given location
  obs <- dat %>% 
    filter(Park == park & o_day %in% occ_days) %>%
    select(Park, StdLocName, Species_code, yr, o_day) %>%
    rename(spp = Species_code) %>%
    distinct
  
  # Attach occasion to each row and then distinct again
  for (i in 1:nrow(obs)) {
    obs$yr_occ[i] <- occasions$yr_occ[obs$o_day[i] >= occasions$start_day & 
                                        obs$o_day[i] <= occasions$end_day]
  }
  # Remove replicate observations of species at a location during each occasion
  obs <- obs %>%
    select(-o_day) %>%
    distinct() %>%
    mutate(detect = 1)
  
  # Extract rows from events matrix that correspond to locations in selected park
  locs_park <- locs$StdLocName[locs$UnitCode == park]
  event_mat_park <- event_mat[rownames(event_mat) %in% locs_park,]
  
  # Extract columns from events matrix that correspond to sampling occasions
  event_mat_park <- event_mat_park[,colnames(event_mat_park) %in% occ_days]
  
  # Summarize event data by occasion 
  event_occ <- matrix(NA, 
                      nrow = nrow(event_mat_park), 
                      ncol = nrow(occasions),
                      dimnames = list(rownames(event_mat_park), occasions$yr_occ))
  
  # Create matrix with 1/0 indicating whether camera was operational during that 
  # occasion
  for (i in 1:ncol(event_occ)) {
    multiday <- event_mat_park[,colnames(event_mat_park) %in% 
                                 occasions$start_day[i]:occasions$end_day[i]]
    event_occ[,i] <- apply(multiday, 1, sum)
    event_occ[event_occ > 1] <- 1
  }
  # check: 
  # table(c(event_occ), useNA = "always")
  
  # Convert to long form
  event_occ_long <- event_occ %>%
    as.data.frame() %>%
    mutate(StdLocName = rownames(.)) %>%
    pivot_longer(cols = !last_col(),
                 names_to = "yr_occ",
                 values_to = "obs") %>%
    data.frame()
  
  # Merge observation and species detection information
  detects <- expand.grid(Park = park,
                         StdLocName = sort(locs_park),
                         yr_occ = occasions$yr_occ,
                         spp = unique(obs$spp))
  detects <- left_join(detects, event_occ_long, by = c("StdLocName", "yr_occ"))
  detects <- left_join(detects, obs[,c("StdLocName", "spp", "yr_occ", "detect")],
                       by = c("StdLocName", "spp", "yr_occ"))
  detects$detect[is.na(detects$detect)] <- 0
  detects$yr <- as.numeric(str_sub(detects$yr_occ, 1, 4))
  
  # Summarize by park and species
  # nobs is number of "observations" (location * occasion when camera operational)
  # ndetects is number of detections of that species
  detects_allyrs <- detects %>%
    group_by(Park, spp) %>%
    summarize(.groups = "keep",
              nobs = sum(obs),
              ndetects = sum(detect)) %>%
    mutate(propdetect = round(ndetects/nobs,2)) %>%
    data.frame()
  detects_allyrs
  
  # Create dataframe with information on species detections at all parks
  if (n == 1) {
    spp_detections <- detects_allyrs
  } else {
    spp_detections <- rbind(spp_detections, detects_allyrs)
  }  
}

spp_detections <- spp_detections %>%
  rename(Species_code = spp) %>%
  left_join(., species[, c(1, 3)], by = "Species_code") %>%
  rename(Species = Common_name) %>%
  select(-Species_code) %>%
  mutate(Park = ifelse(Park == "CHIR", "Chiricahua NM",
                       ifelse(Park == "ORPI", "Organ Pipe Cactus NM",
                              "Saguaro NP, Tucson Mtn District"))) %>%
  relocate(Species, .after = Park) %>%
  arrange(Park, desc(ndetects))

# Export species detections dataframe as csvs
write.csv(spp_detections,
          file = "output/species-detections-bypark-20172022.csv",
          row.names = FALSE)


