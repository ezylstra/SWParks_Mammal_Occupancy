################################################################################
# SODN -- Camera trap data, 2016-2022
# Delineate sampling occasions and calculate the number of species detections
# during sampling occasions at a given park

# ER Zylstra
# Updated: 2023-10-12
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(sf)
library(terra)

# Specify park (CHIR, ORPI, or SAGW)
PARK <- "SAGW"

source("src/photo-data/format-mammal-data.R")

  # dat = information about each photo (date, time, species, location)
  # events = information about each camera deployment (dates, location, duration)
  # event_mat = camera location x day matrix with 1/0 indicating whether camera
  #             was operational or not
  # locs = information about each camera location (park, lat/long, name)
  # species = table with species observed (species code, common name)

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

# Range of photo observations
yr_min <- min(dat$yr)
yr_max = max(dat$yr)

# Extract columns from events data and limit events to only those years when
# we have photo data
events <- events %>%
  select(Park, LocationName, loc, d_yr, active_start, active_end, 
         operational, start_day, end_day) %>%
  filter(d_yr %in% yr_min:yr_max)

# Calculate the number of cameras that are deployed each day
days_park <- days_df
days_park$Park <- PARK
days_park$n_cameras <- colSums(event_mat)
  
# Calculate proportion of cameras in selected park that are deployed each day
days_park$prop_deploy <- days_park$n_cameras / nrow(event_mat)
  
# Identify those dates when the proportion of cameras deployed >= threshold
days_park$at_thresh <- 1*(days_park$prop_deploy >= threshold)
  
# Look at consecutive days with sufficient number of cameras deployed
# data.frame(unclass(rle(days_park$at_thresh)))

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
if (PARK == "CHIR") {
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
            end = max(date),
            .groups = "keep") %>%
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

# Export occasions dataframe as csv
occasions_file <- paste0("data/occasions/occasions-", PARK, ".csv")
write.csv(occasions,
          file = occasions_file,
          row.names = FALSE)

#------------------------------------------------------------------------------#
# Calculate number of species detections
#------------------------------------------------------------------------------#

# Create a list of days included in sampling occasions
occ_days <- NULL
for (i in 1:nrow(occasions)) {
  occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
}
  
# Extract photo observations for park and occasion days
# Retain only one species observation per day at a given location
obs <- dat %>% 
  filter(o_day %in% occ_days) %>%
  select(Park, LocationName, loc, Species_code, yr, o_day) %>%
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

# Extract columns from events matrix that correspond to sampling occasions
event_mat <- event_mat[,colnames(event_mat) %in% occ_days]
  
# Summarize event data by occasion 
event_occ <- matrix(NA, 
                    nrow = nrow(event_mat), 
                    ncol = nrow(occasions),
                    dimnames = list(rownames(event_mat), occasions$yr_occ))
  
# Create matrix with 1/0 indicating whether camera was operational during that 
# occasion
for (i in 1:ncol(event_occ)) {
  multiday <- event_mat[,colnames(event_mat) %in% 
                          occasions$start_day[i]:occasions$end_day[i]]
  event_occ[,i] <- apply(multiday, 1, sum)
  event_occ[event_occ > 1] <- 1
}
# check: 
# table(c(event_occ), useNA = "always")
  
# Convert to long form
event_occ_long <- event_occ %>%
  as.data.frame() %>%
  mutate(loc = rownames(.)) %>%
  pivot_longer(cols = !last_col(),
               names_to = "yr_occ",
               values_to = "obs") %>%
  data.frame()
  
# Merge observation and species detection information
detects <- expand.grid(Park = PARK,
                       loc = sort(locs$loc),
                       yr_occ = occasions$yr_occ,
                       spp = unique(obs$spp))
detects <- left_join(detects, event_occ_long, by = c("loc", "yr_occ"))
detects <- left_join(detects, obs[,c("loc", "spp", "yr_occ", "detect")],
                     by = c("loc", "spp", "yr_occ"))
detects$detect[is.na(detects$detect)] <- 0
detects$yr <- as.numeric(str_sub(detects$yr_occ, 1, 4))
  
# Summarize by year and species
# nobs is number of "observations" (location * occasion when camera operational)
# ndetects is number of detections of that species (for given year)
detects_yr <- detects %>%
  group_by(Park, yr, spp) %>%
  summarize(.groups = "keep",
            nobs = sum(obs),
            ndetects = sum(detect)) %>%
  mutate(propdetect = round(ndetects / nobs, 2)) %>%
  data.frame()

# Summarize by species
# nobs is number of "observations" (location * occasion when camera operational)
# ndetects is number of detections of that species
detects_allyrs <- detects %>%
  group_by(Park, spp) %>%
  summarize(.groups = "keep",
            nobs = sum(obs),
            ndetects = sum(detect)) %>%
  mutate(propdetect = round(ndetects / nobs, 2)) %>%
  arrange(desc(propdetect)) %>%
  mutate(first_yr = yr_min,
         last_yr = yr_max) %>%
  data.frame()
detects_allyrs

# Export species detections dataframes as csvs
allyrs_file <- paste0("output/species-detections-", PARK, ".csv")
write.csv(detects_allyrs,
          file = allyrs_file,
          row.names = FALSE)
byyr_file <- paste0("output/species-detections-byyr-", PARK, ".csv")
write.csv(detects_yr,
          file = byyr_file,
          row.names = FALSE)
