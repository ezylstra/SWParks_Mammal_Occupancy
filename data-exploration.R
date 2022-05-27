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

# Create a park-specific index for camera location (1-65)
events <- events %>% 
  group_by(Park) %>% 
  mutate(locnum = as.numeric(as.factor(StdLocName))) %>%
  as.data.frame()

# Plot events at three main parks, excluding 2016 (different locs sampled at CHIR)
ggplot() + 
  geom_segment(filter(events, Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 0.5, color = 'dodgerblue3') +
  labs(x = 'Date', y = 'Camera number') + 
  facet_grid(rows = vars(Park))
ggsave("output/SamplingEvents_3Parks.jpg", 
       width = 6.5, height = 6.5, 
       units = "in")

# Plot events at four smaller parks
ggplot() + 
  geom_segment(filter(events, !Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 0.5, color = 'dodgerblue3') +
  labs(x = 'Date', y = 'Camera number') + 
  facet_grid(rows = vars(Park))
ggsave("output/SamplingEvents_OtherParks.jpg", 
       width = 6.5, height = 6.5, 
       units = "in")




# TO DO:
# Summarize events - when, length, #detections by species, when detections occur
# (use to figure out when batteries likely died...)
# Summarize detections for common species - by park, year
# (time between detections for common species - 7-day occasion reasonable?)


# Replace the long location name (StdLocName) with a SiteID?
# Remove unknown species (or maybe just some of them?)

