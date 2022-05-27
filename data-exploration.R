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
library(gridExtra)

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
               size = 0.5, color = "dodgerblue3") +
  labs(x = 'Date', y = "Camera number") + 
  facet_grid(rows = vars(Park))
# Save plot in output/folder
  # ggsave("output/SamplingEvents_3Parks.jpg", 
  #        width = 6.5, height = 6.5, 
  #        units = "in")

# Plot events at four smaller parks
ggplot() + 
  geom_segment(filter(events, !Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 0.5, color = "dodgerblue3") +
  labs(x = 'Date', y = "Camera number") + 
  facet_grid(rows = vars(Park))
# Save plot in output/folder
  # ggsave("output/SamplingEvents_OtherParks.jpg", 
  #        width = 6.5, height = 6.5, 
  #        units = "in")

#-------------------------------------------------------------------------------#
# Visualize deployments and photos
#-------------------------------------------------------------------------------#

# SAGW 2018
g1 <- ggplot(filter(dat, Park == "SAGW" & yr == 2018), aes(obsdate)) + 
  geom_histogram(binwidth = 1) +
  labs(x = "", y = "Number of photos") + 
  coord_cartesian(xlim = as.Date(c("2018-01-01", "2018-03-31"))) + 
  annotate("text", x = as.Date("2018-01-01"), y = Inf, hjust=0, vjust=2, label = "2018") +
  theme(text=element_text(size = 8),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = margin(0.2, 0.2, 0, 0.2, unit = "cm"))
g2 <- ggplot() + 
  geom_segment(filter(events, Park == "SAGW" & d_yr == 2018),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 0.5, color = "dodgerblue3") +
  labs(x = "", y = "Camera number") + 
  coord_cartesian(xlim = as.Date(c("2018-01-01", "2018-03-31"))) + 
  theme(text=element_text(size = 8),
        axis.title.x=element_blank(),
        plot.margin = margin(0.1, 0.2, 0.2, 0.2, unit = "cm"))

g12 <- grid.arrange(g1, g2, nrow = 2)

g14 <- grid.arrange(g12, g12, ncol = 2)

# Just testing to see how I can merge multiple plots.
ggsave("output/Test.jpg", 
       g14,
       width = 10, height = 12,
       units = "cm")

# TO DO:
# Summarize events - when, length, #detections by species, when detections occur
# (use to figure out when batteries likely died...)
# Summarize detections for common species - by park, year
# (time between detections for common species - 7-day occasion reasonable?)
