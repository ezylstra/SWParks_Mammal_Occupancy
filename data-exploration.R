################################################################################
# SODN -- Camera trap data, 2016-2022
# Data exploration

# ER Zylstra
# 2022-05-19
################################################################################

library(dplyr)
library(lubridate)
library(stringr)

#-------------------------------------------------------------------------------#
# Run script to import, format data
#-------------------------------------------------------------------------------#
source("format-mammal-data.R")

# rm(list = ls())

#-------------------------------------------------------------------------------#
# Visualize when and where cameras were deployed
#-------------------------------------------------------------------------------#
table(events$Park, events$d_yr)

events$d_yday <- yday(events$d_date)
events$r_yday <- yday(events$r_date)

events_py <- as.data.frame(group_by(events, Park, d_yr) %>% 
  summarize(n_events = length(d_yday), 
            start_date = min(d_yday),
            end_date = max(r_yday)))

events_py



# TO DO:
# Summarize events - when, length, #detections by species, when detections occur
# (use to figure out when batteries likely died...)
# Summarize detections for common species - by park, year
# (time between detections for common species - 7-day occasion reasonable?)


# Replace the long location name (StdLocName) with a SiteID?
# Remove unknown species (or maybe just some of them?)

