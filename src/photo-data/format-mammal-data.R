################################################################################
# SODN -- Camera trap data, 2016-2022
# Import and format data

# ER Zylstra
# 2023-10-04
################################################################################

# Need to load these packages if not calling this script via source()

library(dplyr)
library(lubridate)
library(stringr)
library(sf)
library(terra)

#------------------------------------------------------------------------------#
# Import data
#------------------------------------------------------------------------------#

# Species list
species_list <- read.csv("data/mammals/PROTECTED_SpeciesList.csv")

# Observations
dat <- read.csv("data/mammals/PROTECTED_Detections.csv")

# Camera locations
locs_ann <- vect("data/mammals/PROTECTED_CameraLocations_Annual.shp")
names(locs_ann) <- c("UnitCode", "StdLocName", "LocationName", "DeployDate",
                     "StdLocName_Flag", "LocationName_Flag", "DeployDate_Flag",
                     "geometry_Flag")

# Deployment schedule
events <- read.csv("data/mammals/PROTECTED_Events.csv")

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
# include in occupancy models. Note that there are 2 entries for UNCA, so we'll 
# remove the one with Common_Name = "unknown canid" and TSN = NA.
exclude <- c("Harris's antelope squirrel", "Merriam's kangaroo rat", 
             "round-tailed ground squirrel", "unknown animal", 
             "unknown canid", "unknown kangaroo rat", "unknown rodent", 
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
  rowwise() %>%
  mutate(reject_sum = sum(c_across(StdLocName_Flag:ActiveEnd_Flag) == "R")) %>%
  filter(reject_sum == 0) %>%
  select(-c(contains("Flag"), reject_sum)) %>%
  data.frame()

# Only keep necessary columns and remove any events that aren't associated with 
# the 3 focal parks:
events <- events %>%
  select(-c(StdLocName, CrewRetrieve)) %>%
  filter(UnitCode %in% c("CHIR", "ORPI", "SAGW"))

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

# checks
  # count(events, CrewDeploy, deploy_exp)

#------------------------------------------------------------------------------#
# Format and organize mammal observation data
#------------------------------------------------------------------------------#

# Remove rows with StudyAreaName == "OTHER" (these are not actually photos,
# but instead are information from an extra column in the csv when there's
# an entry for a lost or damaged camera)
dat <- filter(dat, StudyAreaName != "OTHER")

# Remove empty or unnecessary columns:
dat <- select(dat,-c(StudyAreaName, StudyAreaID, 
                     UTM_E, UTM_N, UTMZone, 
                     FileName, VisitID, ImgID, ImageNum, Highlight,
                     SpeciesID, DetailText, Individuals))

# Rename StudyAreaAbbr column and change MOCA entries to MOCC to be consistent
# with other data sources
dat <- dat %>%
  rename(Park = StudyAreaAbbr) %>%
  mutate(Park = replace(Park, Park == "MOCA", "MOCC"))

# Exclude photo observations of non-mammals or rodents
birds <- c("bird species", "unknown bird", "unknown sparrow", "unknown raven",
           "unknown woodpecker", "unknown jay", "unknown kingbird", 
           "zone-tailed hawk", "red-tailed hawk", "turkey vulture",
           "great-horned owl", "western screech-owl", "barn owl", 
           "lesser nighthawk", "greater roadrunner", "wild turkey", 
           "gilded flicker", "northern flicker", "gila woodpecker",
           "great-tailed grackle", "common raven", "northern mockingbird",  
           "mexican jay", "california scrub-jay", "curve-billed thrasher",
           "eurasian collared dove", "white-winged dove", "mourning dove", 
           "northern cardinal", "gambel's quail", "montezuma quail",
           "cactus wren", "canyon wren", "bewick's wren", "rock wren",
           "dark-eyed junco", "black-throated sparrow",  
           "rufous-crowned sparrow", "canyon towhee", "spotted towhee")

herps <- c("unknown reptile", "gophersnake", "clark's spiny lizard", 
           "eastern collared lizard", "western whiptail", "desert tortoise", 
           "side-blotched lizard", "common side-blotched lizard")

mammals_to_exclude <- c("human", "unknown rodent", "unknown kangaroo rat", 
                        "unknown deer mouse", "unknown woodrat", 
                        "desert kangaroo rat", "merriam's kangaroo rat", 
                        "white-throated woodrat", "cliff chipmunk", 
                        "round-tailed ground squirrel",
                        "harris's antelope squirrel")

other <- c("datasheet", "camera lost", "vehicle", "",
           "tarantula", "unknown animal")

dat <- dat %>%
  mutate(CommonName_lower = tolower(CommonName)) %>%
  filter(! CommonName_lower %in% c(birds,
                                   herps, 
                                   mammals_to_exclude,
                                   other)) %>%
  select(-CommonName_lower) 

# Species
  # Create separate table with species information
  species <- count(dat, 
                   CommonName, 
                   paste0(dat$Genus, " ", dat$Species), 
                   ShortName)
  colnames(species) <- c("Common_name", "Species", "Species_code", "n")
   # 35 different "species" -- includes Unknowns
  species[species$n %in% range(species$n[-which(grepl("Unk", 
                                                      species$Common_name))]),]
    # Number of observations per known species ranges from 
    # 10 (Arizona gray squirrel) to >13,000 (WT deer) 
  # Remove Species and Common_name variables from observations dataframe 
  # and rename column with species codes
  dat <- dat %>%
    select(-c(Species, Genus, CommonName)) %>%
    rename(Species_code = ShortName)

# Extract year from ImgPath (used to organize/file photos)
summary(n.backslashes <- str_count(dat$ImgPath, "\\\\"))  
  # ImgPath always has 7 backslashes (ie, 8 character strings)
n.strings <- mean(n.backslashes) + 1
dat <- cbind(dat, str_split_fixed(dat$ImgPath, "\\\\", n.strings)[,6])
names(dat)[ncol(dat)] <- "FY_filepath"
  # Check:
  count(dat, FY_filepath)  #7 FYs (16-22)

# Remove ImgPath column
dat <- select(dat, -ImgPath)

# Format date and time
  # Create new date-time column
  dat$datetime <- parse_date_time(dat$ImageDate, 
                                  orders = c("%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y"))
  # Create new date-, time-related columns
  dat <- dat %>%
    mutate(obsdate = date(datetime),
           yr = year(datetime),
           mon = month(datetime),
           yday = yday(datetime),
           obstime = hour(datetime) + minute(datetime) / 60 + second(datetime) / 3600)
  
  # Remove ImageDate column
  dat <- select(dat, -ImageDate)
  
#------------------------------------------------------------------------------#
# Attach spatial data to observations
#------------------------------------------------------------------------------#

# Some lat/longs in dat are wrong (for 7 ORPI cameras -- will be fixed in a 
# future iteration), and some lat/longs are missing for smaller parks.  
# For now, we'll use lat/longs from the locs csv

# Will create a simple location name (loc_short) to match things up between 
# dat and locs. Most entries will look like: PARK_XXX_XXX
# Note: do want to keep StdLocName column in locs because that's exactly what 
# appears in the events dataframe

dat <- dat %>%
  mutate(loc_short = paste0(Park, "_", LocationName))
  
locs <- locs %>% 
  # Remove GICL locations
  filter(UnitCode != "GICL") %>%
  # Remove locations that don't appear in events dataframe 
  filter(StdLocName %in% events$StdLocName) %>%
  # Remove leading "WBC_", "V", and "W" from MarkerName where they appear
  mutate(loc_short_np = ifelse(str_detect(MarkerName, "WBC"), 
                               str_replace(MarkerName, fixed ("WBC_"), ""),
                               ifelse(str_sub(MarkerName, 1, 1) %in% c("V", "W"),
                                      str_sub(MarkerName, 2, -1),
                                      MarkerName))) %>%
  # Add the park code as a prefix
  mutate(loc_short = paste0(UnitCode, "_", loc_short_np)) %>%
  select(-loc_short_np)
  
# Check that each location has a unique set of lat/longs:
# dim(unique(locs[,c("POINT_X", "POINT_Y", "loc_short")]))
# dim(unique(locs[,c("POINT_X", "POINT_Y")]))

# Attach correct lat/longs to photo observations dataframe
dat <- left_join(dat, select(locs, c(loc_short, StdLocName, POINT_X, POINT_Y)),
                 by = "loc_short")
# Check that every photo observation has lat/longs:
# summary(dat[, c("POINT_X", "POINT_Y")])

# Remove LatitudeDD, LongitudeDD columns and rename POINT_X, POINT_Y
dat <- dat %>%
  select(-c(LatitudeDD, LongitudeDD)) %>%
  rename(longitude = POINT_X,
         latitude = POINT_Y)
  
#------------------------------------------------------------------------------#
# Linking events file to observations 
#------------------------------------------------------------------------------# 

# Summarize sampling events by park and year (and compare to photo obs data)
table(events$Park, events$d_yr)
table(dat$Park, dat$yr)

# Check that all cameras in photo obs dataset appear in the events dataset
obslocs <- sort(unique(dat$StdLocName))
obslocs[!obslocs %in% events$StdLocName]  # Yes, all appear in events
  
# Convert deployment/retrieval dates to integers, setting Jan 1 2016 equal to 1
events$d_day <- as.numeric(events$d_date) - as.numeric(as.Date("2015-12-31"))
events$r_day <- as.numeric(events$r_date) - as.numeric(as.Date("2015-12-31"))

# Create location by date matrix 
# (1 indicates date was during a sampling event, 0 otherwise)
  
  # Create matrix with all 0s
  eventlocs <- sort(unique(events$StdLocName))
  event_mat <- matrix(0, nrow = length(eventlocs), ncol = max(events$r_day)) 
  rownames(event_mat) <- eventlocs
  colnames(event_mat) <- 1:max(events$r_day)
  
  # Replace 0s with 1s during each sampling event
  for (i in 1:nrow(event_mat)) {
    temp <- events[events$StdLocName == eventlocs[i], ]
    
    for (j in 1:nrow(temp)) {
      event_mat[i, temp$d_day[j]:temp$r_day[j]] <- 1
    }
  }
  # Checks at 1 random location in each of 3 bigger parks
  # (only 1s during sampling event and 0s outside of event?):
  # # CHIR
  # events[events$StdLocName == eventlocs[28], 
  #        c("StdLocName", "d_date", "r_date", "d_day","r_day")]
  # sum(event_mat[28, 1:629] == 1); sum(event_mat[28, 1:629] == 0) 
  # sum(event_mat[28, 630:699] == 1); sum(event_mat[28, 630:699] == 0) 
  # # ORPI
  # events[events$StdLocName == eventlocs[111], 
  #        c("StdLocName", "d_date", "r_date", "d_day","r_day")]
  # sum(event_mat[111, 1:112] == 1); sum(event_mat[111, 1:112] == 0) 
  # sum(event_mat[111, 113:168] == 1); sum(event_mat[111, 113:168] == 0)
  # # SAGW
  # events[events$StdLocName == eventlocs[190], 
  #        c("StdLocName", "d_date", "r_date", "d_day","r_day")]  
  # sum(event_mat[190, 1883:2210] == 1); sum(event_mat[190, 1883:2210] == 0) 
  # sum(event_mat[190, 2211:2247] == 1); sum(event_mat[190, 2211:2247] == 0)
  
# Create character strings with location and date for each mammal observation
dat$o_day <- as.numeric(dat$obsdate) - as.numeric(as.Date("2015-12-31"))
dat$locdate <- paste(dat$StdLocName, dat$o_day, sep = "_")

# Create character strings with location and date for day during each sampling event
eventvec <- as.character(vector())
for (i in 1:length(eventlocs)) {
  eventvec <- append(eventvec, paste(eventlocs[i], which(event_mat[i,] == 1), sep = "_"))
}
  # check:
  head(eventvec); head(events[,c("StdLocName", "d_day", "r_day")])
  
# Check that all photo dates were during listed sampling events
summary(dat$locdate %in% eventvec) 
  # Yes, all photo dates occur during known event

#-----------------------------------------------------------------------------------#
# Remove objects that are no longer needed
#-----------------------------------------------------------------------------------# 

rm(deploys, mowe_add, temp, birds, eventlocs, eventvec, exp1, exp2, herps, i, j,
   mammals_to_exclude, mowe_locs, n.backslashes, n.strings, obslocs, other)
