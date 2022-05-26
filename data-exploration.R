################################################################################
# SODN -- Camera trap data, 2016-2022
# Data exploration

# ER Zylstra
# 2022-05-19
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)

# rm(list = ls())

#-------------------------------------------------------------------------------#
# Run script to import, format data
#-------------------------------------------------------------------------------#

source("format-mammal-data.R")

#-------------------------------------------------------------------------------#
# Visualize when and where cameras were deployed
#-------------------------------------------------------------------------------#

# More generally
events <- events %>% 
  group_by(Park) %>% 
  mutate(locnum = as.numeric(as.factor(StdLocName))) %>%
  as.data.frame()

# Plot SAGW events
ggplot() +
  geom_segment(filter(events, Park == "SAGW"),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 1, color = 'dodgerblue3') +
  labs(x = 'Date', y = 'Camera number')



# TO DO:
# Summarize events - when, length, #detections by species, when detections occur
# (use to figure out when batteries likely died...)
# Summarize detections for common species - by park, year
# (time between detections for common species - 7-day occasion reasonable?)


# Replace the long location name (StdLocName) with a SiteID?
# Remove unknown species (or maybe just some of them?)

