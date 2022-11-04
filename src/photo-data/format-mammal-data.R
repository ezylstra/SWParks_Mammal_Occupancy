################################################################################
# SODN -- Camera trap data, 2016-2022
# Import and format data

# ER Zylstra
# 2022-04-28
################################################################################

# Need to load these packages if not calling this script via source()

library(dplyr)
library(lubridate)
library(stringr)

rm(list = ls())

#------------------------------------------------------------------------------#
# Import data
#------------------------------------------------------------------------------#

# Observations
dat <- read.csv("data/mammals/MAMMALS_ALL_2022-08-12.csv")

# Camera locations
locs <- read.csv("data/mammals/SODN_Wildlife_Locations_XY_Revised_20220502.csv")[,2:9]

# Deployment schedule
events <- read.csv("data/mammals/SODN_Wildlife_Events_20221025.csv")[,c(3:20, 28:29)]

# Experience level of personnel deploying cameras
deploys <- read.csv("data/covariates/deployment-personnel.csv")

#------------------------------------------------------------------------------#
# Format events data
#------------------------------------------------------------------------------#

# Notes about sampling "events":
# Occasionally cameras were immediately re-deployed for continuous sampling
# Occasionally two cameras were deployed at the same location simultaneously
# Sometimes cameras left out for >1 yr. Not sure how long they collected photos.

# Only keep necessary columns
events <- select(events, c(StdLocName, ProtocolVersion, DeployDate, 
                           RetrievalDate, CameraName, MountMethod, 
                           BatteryStatus, CameraSensitivity, 
                           DelaySec, ImagePer, TotalPics, CrewDeploy))

# Add column to identify Park and remove information about parks that aren't in 
# photo observations dataset (GICL and National Wildlife refuges)
events <- events %>%
  mutate(Park = str_split_fixed(events$StdLocName, "_", 4)[,2]) %>%
  filter(!Park %in% c("GICL", "LCNWR", "SBNWR"))

# Create new deployment date-time column and remove DeployDate
events <- events %>%
  mutate(d_datetime = parse_date_time(events$DeployDate, 
                                      orders = c("%m/%d/%Y %H:%M", 
                                                 "%m/%d/%Y %H:%M:%S"))) %>%
  select(-DeployDate)

# Create new retrieval date-time column and remove RetrievalDate
events <- events %>%
  mutate(r_datetime = parse_date_time(events$RetrievalDate, 
                                      orders = c("%m/%d/%Y %H:%M", 
                                                 "%m/%d/%Y %H:%M:%S"))) %>%
  select(-RetrievalDate)

  #-- Fix known issues in events dataset --------------------------------------#

  # Change the retrieval date for CHIR camera 502-003 deployed in 2019
  events <- events %>%
    mutate(r_datetime = if_else(StdLocName == "Wildlife_CHIR_V502_003" & 
                                  year(d_datetime) == 2019,
                                parse_date_time("2021-05-12",
                                                orders = "%Y-%m-%d"),
                                r_datetime))  

  # Change the retrieval date for all cameras deployed in 2019 at MOCC
  events <- events %>%
    mutate(r_datetime = if_else(Park == "MOCC" & year(d_datetime) == 2019,
                                       parse_date_time("2021-10-14",
                                                       orders = "%Y-%m-%d"),
                                       r_datetime))

  # Add two deployments for all cameras as MOWE
  mowe_locs <- events$StdLocName[events$Park == "MOWE"]
  mowe_add <- data.frame(matrix(NA, 
                                ncol = ncol(events), 
                                nrow = 2 * length(mowe_locs)))
  colnames(mowe_add) <- colnames(events)
  mowe_add <- mowe_add %>%
    mutate(StdLocName = rep(mowe_locs, 2),
           Park = "MOWE",
           d_datetime = parse_date_time(rep(c("2018-08-07", "2019-05-29"), 
                                            each = length(mowe_locs)),
                                        orders = "%Y-%m-%d"),
           r_datetime = parse_date_time(rep(c("2019-05-29", "2021-10-14"), 
                                            each = length(mowe_locs)),
                                        orders = "%Y-%m-%d"))                                                               
  events <- rbind(events, mowe_add)
  #----------------------------------------------------------------------------#

# Create new date and year columns
events <- events %>%
  mutate(d_date = date(d_datetime),
         d_yr = year(d_datetime),
         r_date = date(r_datetime),
         r_yr = year(r_datetime)) 

# Reorder columns and sort events dataframe by location and date
events <- relocate(events, c(Park, d_date, r_date), .after = StdLocName)
events <- arrange(events, StdLocName, d_datetime)

# Restrict events dataframe to 2016 forward
events <- filter(events, d_yr > 2015)

# Calculate length of deployment, in days
events$duration <- as.double(difftime(as.POSIXct(events$r_datetime), 
                                      as.POSIXct(events$d_datetime), 
                                      units = 'days'))
  # Summarize/Visualize 
  # summary(events$duration)
  # hist(events$duration, breaks = 25)

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
  # There are sampling events in CHIR in 2016 with no corresponding observations
  # Sampling methods were different that year, so probably ok not to track 
  # these down

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
  # checks at random locations 
  # (only 1s during sampling event and 0s outside of event?):
  events[events$StdLocName == eventlocs[206], 
         c("StdLocName", "d_date", "r_date", "d_day","r_day")]
  sum(event_mat[206, 1:390] == 1); sum(event_mat[206, 1:390] == 0) 
  sum(event_mat[206, 391:424] == 1); sum(event_mat[206, 391:424] == 0) 
  
  events[events$StdLocName == eventlocs[111], 
         c("StdLocName", "d_date", "r_date", "d_day","r_day")]
  sum(event_mat[111, 1:115] == 1); sum(event_mat[111, 1:115] == 0) 
  sum(event_mat[111, 810:853] == 1); sum(event_mat[111, 810:853] == 0) 
  
# Create character strings with location and date for each mammal observation
dat$o_day <- as.numeric(dat$obsdate) - as.numeric(as.Date("2015-12-31"))
dat$locdate <- paste(dat$StdLocName, dat$o_day, sep = "_")

# Create character strings with location and date for day during each sampling event
eventvec <- as.character(vector())
for (i in 1:length(eventlocs)) {
  eventvec <- append(eventvec, paste(eventlocs[i], which(event_mat[i,] == 1), sep = "_"))
}
  # check:
  head(eventvec); head(events[,c(1, 13:14, 18, 17)])
  
# Check that all photo dates were during listed sampling events
summary(dat$locdate %in% eventvec) 
  # Yes, all photo dates occur during known event

#-----------------------------------------------------------------------------------#
# Remove objects that are no longer needed
#-----------------------------------------------------------------------------------# 

rm(list = setdiff(ls(), c("dat", "locs", "events", "event_mat", "species")))
