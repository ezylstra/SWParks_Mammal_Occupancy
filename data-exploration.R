################################################################################
# SODN -- Camera trap data, 2016-2022
# Data exploration

# ER Zylstra
# 2022-05-19
################################################################################

#-------------------------------------------------------------------------------#
# Run script to import, format data
#-------------------------------------------------------------------------------#
source("format-mammal-data.R")

# rm(list = ls())

# TO DO:
# Summarize events - when, length, #detections by species, when detections occur
# (use to figure out when batteries likely died...)
# Summarize detections for common species - by park, year
# (time between detections for common species - 7-day occasion reasonable?)


# Replace the long location name (StdLocName) with a SiteID?
# Remove unknown species (or maybe just some of them?)

