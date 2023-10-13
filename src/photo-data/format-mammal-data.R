################################################################################
# SODN -- Camera trap data, 2016-2022
# Import and format data

# ER Zylstra
# 2023-10-12
################################################################################

# Need to load these packages if not calling this script via source()

library(dplyr)
library(lubridate)
library(stringr)
library(sf)
library(terra)

# Need to specify park if not calling this script via source()
# PARK <- "SAGW"

#------------------------------------------------------------------------------#
# Import data
#------------------------------------------------------------------------------#

# Species list
species_list <- read.csv("data/mammals/PROTECTED_SpeciesList.csv")

# Observations
dat <- read.csv(paste0("data/mammals/PROTECTED_Detections_", PARK, ".csv"))

# Camera locations
locs_ann <- vect(paste0("data/mammals/PROTECTED_CameraLocations_Annual_",
                        PARK, ".shp"))
names(locs_ann) <- c("UnitCode", "StdLocName", "LocationName", "DeployDate",
                     "StdLocName_Flag", "LocationName_Flag", "DeployDate_Flag",
                     "geometry_Flag")

  # Logical indicating whether to save a new shapefile with central location
  # for each camera across years (should only need to do this once per year)
  centroid_save <- FALSE

# Deployment schedule
events <- read.csv(paste0("data/mammals/PROTECTED_Events_", PARK, ".csv"))

# Experience level of personnel deploying cameras
deploys <- read.csv("data/covariates/deployment-personnel.csv")

#------------------------------------------------------------------------------#
# Notes on "Flag" columns
#------------------------------------------------------------------------------#

# Flag columns will have an acceptance rating and may or may not have a
# qualification code in brackets. 

# Acceptance ratings:
# A: accepted. Meets data quality standards. (default)
# AE: accepted, estimated. Data have been estimated or corrected to ensure that 
  # data meet quality standard. 
# R: rejected. Data do not meet quality standards and should probably be removed
  # before analysis.
# P: provisional/preliminary. Data subject to change based on quality assurance 
  # or quality control processes.

#------------------------------------------------------------------------------#
# Format species list
#------------------------------------------------------------------------------#

# Will remove all species except for medium to large mammals that we want to 
# include in occupancy models. 
exclude <- c("Harris's antelope squirrel", "Merriam's kangaroo rat", 
             "round-tailed ground squirrel", "unknown animal", 
             "unknown kangaroo rat", "unknown rodent", 
             "unknown woodrat", "western white-throated woodrat")

species <- species_list %>%
  filter(!Common_Name %in% exclude) %>%
  rename(Species_code = Accepted_Code,
         Species = Scientific_Name,
         Common_name = Common_Name, 
         Nativeness = Nativity) %>%
  select(Species_code, Species, Common_name, TSN, Family, Nativeness, Protected)

#------------------------------------------------------------------------------#
# Format events data
#------------------------------------------------------------------------------#

# Notes about sampling "events" (may or may not be relevant with new files)
# Occasionally cameras were immediately re-deployed for continuous sampling
# Sometimes cameras left out for >1 yr. Not sure how long they collected photos.

# Checked for flagged data
count(events, events[, grep("Flag", colnames(events))])
# Remove any events that have one or more Flags = R (Reject)
events <- events %>%
  mutate(across(ends_with("Flag"), function(x) ifelse(x == "R", 1, 0))) %>%
  mutate(reject_sum = rowSums(select(., ends_with("Flag")))) %>%
  filter(reject_sum == 0)

# Only keep necessary columns and remove any events that aren't associated with 
# the focal park:
events <- events %>%
  select(-c(StdLocName, CrewRetrieve)) %>%
  filter(UnitCode == PARK)

# Convert deployment, retrieval, active dates to date objects, and check that 
# active start/ends are always within deployment dates
events <- events %>%
  mutate(d_date = ymd(DeployDate),
         r_date = ymd(RetrievalDate),
         active_start = ymd(ActiveStart),
         active_end = ymd(ActiveEnd),
         actst_check = ifelse(active_start < d_date, 1, 0),
         actend_check = ifelse(active_end > r_date, 1, 0))
if(sum(events$actst_check) > 0 | sum(events$actend_check) > 0) {
  stop("One or more active dates fall outside of deployment window.\n")
}
# If no issues, remove checks columns:
events <- select(events, -c(DeployDate, RetrievalDate, ActiveStart, ActiveEnd, 
                            actst_check, actend_check))

  #-- Fix known issues in events dataset --------------------------------------#
  # DON'T KNOW IF WE'LL NEED ANY OF THIS WITH NEW DATASET  
  #
  # # At a few locations at ORPI in 2021, two cameras were deployed at the same 
  # # location simultaneously (removing event information for 1 of the cameras)
  # events <- events %>%
  #   filter(!(UnitCode == "ORPI" & LocationName == "V101_16W" & 
  #              year(d_date) == 2021 & CameraName == "SODN_040")) %>%
  #   filter(!(UnitCode == "ORPI" & LocationName ==  "V102_107W" & 
  #              year(d_date) == 2021 & CameraName == "SODN_134")) %>%
  #   filter(!(UnitCode == "ORPI" & LocationName ==  "V103_06W" & 
  #              year(d_date) == 2021 & CameraName == "SODN_167"))     
  # 
  # # Change the retrieval date for CHIR camera 502-003 deployed in 2019
  # events <- events %>%
  #   mutate(r_date = if_else(UnitCode == "CHIR" & LocationName == "V502_003" & 
  #                             year(d_date) == 2019,
  #                           parse_date_time("2021-05-12", 
  #                                           orders = "%Y-%m-%d"),
  #                           r_date))  
  # 
  #----------------------------------------------------------------------------#

# Create new year columns
events <- events %>%
  mutate(d_yr = year(d_date),
         r_yr = year(r_date)) 

# Reorder columns and sort events dataframe by location and date
events <- events %>%
  relocate(c(d_date, r_date, active_start, active_end), .after = LocationName) %>%
  arrange(UnitCode, LocationName, d_date)

# Restrict events dataframe to 2016 forward for ORPI, 2017 forward for other
# parks (CHIR has events listed in 2016 but no corresponding photos. Sampling
# methods were different in 2016, so probably ok not to track those photos down)
events <- events %>%
  filter(!(UnitCode == "ORPI" & d_yr < 2016)) %>%
  filter(!(UnitCode != "ORPI" & d_yr < 2017))

# Calculate length of deployment, in days
events$duration <- as.double(difftime(as.POSIXct(events$r_date), 
                                      as.POSIXct(events$d_date), 
                                      units = 'days'))

# Calculate length of time camera was operational
events$operational <- as.double(difftime(as.POSIXct(events$active_end), 
                                         as.POSIXct(events$active_start), 
                                         units = 'days'))

# Summarize/Visualize
# summary(events$operational)
# hist(events$operational, breaks = 25)
# table(events$operational < events$duration)

# Look at events when camera was operational for < 15 days
# events %>% filter(operational < 15) %>%
#   select(d_date, active_start, r_date, active_end, LocationName, BatteryStatus,
#          TotalPics, operational) %>%
#   arrange(d_date)

#------------------------------------------------------------------------------#
# Format and organize information about deployment personnel
#------------------------------------------------------------------------------#  

# Identify deployment personnel that are experts or experienced
exp2 <- deploys$personnel[deploys$expertise == "expert"]
exp1 <- deploys$personnel[deploys$expertise == "experienced"]

# Create deploy_exp variable for each camera deployment with: 
  # 2 = at least one expert present 
  # 1 = at least one experienced person present
  # 0 = all novices
events <- events %>%
  mutate(deploy_exp = ifelse(str_detect(CrewDeploy, 
                                        paste(exp2, collapse = "|")),
                             2, ifelse(str_detect(CrewDeploy, 
                                                  paste(exp1, collapse = "|")),
                                       1, 0)),
         deploy_exp = ifelse(is.na(deploy_exp), 0, deploy_exp))

#------------------------------------------------------------------------------#
# Format and organize mammal observation data
#------------------------------------------------------------------------------#

# Checked for flagged data
count(dat, dat[, grep("Flag", colnames(dat))])
# Remove any detections that have one or more Flags = R (Reject)
dat <- dat %>%
  mutate(across(ends_with("Flag"), function(x) ifelse(x == "R", 1, 0))) %>%
  mutate(reject_sum = rowSums(select(., ends_with("Flag")))) %>%
  filter(reject_sum == 0)

# Exclude photo observations of mammals that aren't in our species list
dat <- filter(dat, Accepted_Code %in% species$Species_code)

# Remove unnecessary columns (including all columns with species info except 
# Accepted_Code) and remove information that isn't associated with the focal 
# park:
dat <- dat %>%
  select(UnitCode, LocationName, ImageDate, Accepted_Code) %>%
  rename(Species_code = Accepted_Code) %>%
  filter(UnitCode == "SAGW")

# Create new date-, time-related columns
dat$datetime <- parse_date_time(dat$ImageDate, orders = c("%m/%d/%Y %H:%M:%S"))
dat <- dat %>%
  mutate(obsdate = date(datetime),
         yr = year(datetime),
         mon = month(datetime),
         yday = yday(datetime),
         time24 = hour(datetime) + minute(datetime) / 60 + second(datetime) / 3600) %>%
  select(-ImageDate)

# Finally, as an extra check, remove detections that occur outside active dates
# (might be able to remove this eventually). Note that the events data now 
# include an ImageDate_Flag that marks these instances as "R" and they are being 
# removed above. But OK to leave in for now as a double-check.
dat <- dat %>%
  left_join(events[, c("UnitCode", "LocationName", 
                       "active_start", "active_end", "d_yr")], 
            by = c("UnitCode", "LocationName", "yr" = "d_yr")) %>%
  filter(obsdate >= active_start & obsdate <= active_end) %>%
  select(-c(active_start, active_end))

#------------------------------------------------------------------------------#
# Attach spatial data to detections
#------------------------------------------------------------------------------#

# Reproject to use the same crs as other objects used in the project 
# (EPSG:4269; lon/lat NAD83)
locs_ann <- terra::project(locs_ann, "EPSG:4269")

# Checked for flagged data
locs_ann_df <- as.data.frame(locs_ann)
locs_ann_df$reject <- 
  (rowSums(locs_ann_df[,endsWith(names(locs_ann_df), "Flag")] == "R") >= 1)
if (sum(locs_ann_df$reject) > 0) {
  stop("One or more camera locations has been flagged.\n")
}

# Calculate the centroid of annual deployment locations for each camera
ann_sf <- sf::st_as_sf(locs_ann)
centroids_sf <- ann_sf %>%
  group_by(UnitCode, LocationName) %>%
  summarize(geometry = st_union(geometry),
            .groups = "keep") %>%
  st_centroid()

# Creating a simple location name (loc) to match things up between dat and locs. 
# Most (but not all) entries will look like: PARK_XXX_XXX. If there is a leading 
# "V" or "W" on LocationName in locs file (like there was in old version), 
# should remove first.
centroids_sf <- centroids_sf %>%
  mutate(LocationName = ifelse(str_sub(LocationName, 1, 1) %in% c("V", "W"),
                               str_sub(LocationName, 2, nchar(LocationName)),
                               LocationName)) %>%
  mutate(loc = paste0(UnitCode, "_", LocationName))

# Save shapefile
if (centroid_save) {
  centroid_file <- paste0("data/mammals/PROTECTED_CameraLocations_Centroids_",
                          PARK, ".shp")
  writeVector(vect(centroids_sf),
              centroid_file,
              overwrite = TRUE) 
}

# Create dataframe with location info and attach coordinates to dat
coords <- st_coordinates(centroids_sf)
locs <- as.data.frame(centroids_sf) %>%
  cbind(coords) %>%
  select(-geometry) %>%
  rename(longitude = X,
         latitude = Y) %>%
  relocate(loc, .before = "UnitCode") %>%
  arrange(loc)
dat <- dat %>%
  mutate(LocationName = ifelse(str_sub(LocationName, 1, 1) %in% c("V", "W"),
                               str_sub(LocationName, 2, nchar(LocationName)),
                               LocationName)) %>%
  mutate(loc = paste0(UnitCode, "_", LocationName)) %>%
  left_join(locs[, c("loc", "longitude", "latitude")], by = "loc")

# Check that each detection has coordinates
# sum(is.na(dat$longitude)) == 0

#------------------------------------------------------------------------------#
# Create a matrix with period each camera was operational
#------------------------------------------------------------------------------# 

# Summarize sampling events by park and year (and compare to detections data)
# table(events$UnitCode, events$d_yr)
# table(dat$UnitCode, dat$yr)

# Create short location name in events
events <- events %>%
  mutate(LocationName = ifelse(str_sub(LocationName, 1, 1) %in% c("V", "W"),
                               str_sub(LocationName, 2, nchar(LocationName)),
                               LocationName)) %>%
  mutate(loc = paste0(UnitCode, "_", LocationName))

# Check that all locations in detections dataset appear in the events dataset
detlocs <- sort(unique(dat$loc))  
detlocs[!detlocs %in% events$loc] # Yes, all appear in events

# Convert active start/end dates to integers, setting Jan 1 2016 equal to 1
events$start_day <- as.numeric(events$active_start) - as.numeric(as.Date("2015-12-31"))
events$end_day <- as.numeric(events$active_end) - as.numeric(as.Date("2015-12-31"))

# Create location by date matrix 
# (1 indicates date was during a sampling event, 0 otherwise)
  
  # Create matrix with all 0s
  eventlocs <- sort(unique(events$loc))
  event_mat <- matrix(0, nrow = length(eventlocs), ncol = max(events$end_day)) 
  rownames(event_mat) <- eventlocs
  colnames(event_mat) <- 1:max(events$end_day)
  
  # Replace 0s with 1s during each sampling event (ie, active window)
  for (i in 1:nrow(event_mat)) {
    temp <- events[events$loc == eventlocs[i], ]
    
    for (j in 1:nrow(temp)) {
      event_mat[i, temp$start_day[j]:temp$end_day[j]] <- 1
    }
  }

# Create character strings with location and date for each mammal detection
dat$o_day <- as.numeric(dat$obsdate) - as.numeric(as.Date("2015-12-31"))
dat$locday <- paste(dat$loc, dat$o_day, sep = "_")

# Create character strings with location and date for day during each active window
eventvec <- as.character(vector())
for (i in 1:length(eventlocs)) {
  eventvec <- append(eventvec, paste(eventlocs[i], which(event_mat[i,] == 1), sep = "_"))
}
  # check:
  # head(eventvec); head(events[,c("loc", "start_day", "end_day")])
  
# Check that all detections fall within active windows
summary(dat$locday %in% eventvec) 
  # Yes, all detections fall within active dates

#------------------------------------------------------------------------------#
# Final formatting
#------------------------------------------------------------------------------# 

# Remove locday column from dat and rename UnitCode to Park
dat <- dat %>% 
  select(-locday) %>% 
  rename(Park = UnitCode)

# Rename UnitCode in other dataframes
events <- rename(events, Park = UnitCode)
locs <- rename(locs, Park = UnitCode)

#------------------------------------------------------------------------------#
# Remove objects that are no longer needed
#------------------------------------------------------------------------------# 

rm(ann_sf, centroids_sf, coords, deploys, locs_ann, locs_ann_df, species_list,
   temp, centroid_save, detlocs, eventlocs, eventvec, exclude, exp1, exp2, i, j)
