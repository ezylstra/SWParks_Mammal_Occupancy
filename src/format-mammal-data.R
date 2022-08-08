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
dat <- read.csv("data/mammals/MAMMALS_ALL_2022-08-01.csv")

#orpi_replace <- read.csv("data/mammals/ORPI_2020_102_54W_adjustedDates.csv")[,c(1:2,4:14,16:19)]
#  colnames(orpi_replace)[c(1:13,16:17)] <- colnames(dat)[c(1:13,15:16)]
  
# Camera locations
locs <- read.csv("data/mammals/SODN_Wildlife_Locations_XY_Revised_20220502.csv")[,2:9]

# Deployment schedule
events <- read.csv("data/mammals/SODN_Wildlife_Events_Revised_20220512.csv")[,3:26]

#------------------------------------------------------------------------------#
# Replace data from one ORPI camera in 2020 with data that has corrected date 
#------------------------------------------------------------------------------#

# ORPI date adjustment should no longer be necessary, because of new csv

# Remove datasheet photo
#orpi_replace <- filter(orpi_replace, Common_name != 'Datasheet')

# Fix a few species names/codes
#orpi_replace <- mutate(orpi_replace, 
#                       Species = ifelse(Common_name == "Unknown fox", 
#                                        "Caninae sp.",
#                                        paste(Genus, Species, sep = " ")),
#                       Species_code = ifelse(Common_name == "Unknown fox",
#                                             "UNFO", 
#                                             Species_code))
#orpi_replace <- select(orpi_replace, !Genus)

# Check that number of 2020 observations at ORPI_102_54W in dat is the same as replacement file
#nrow(filter(dat, str_detect(ImgPath, "ORPI_102_54W") & FieldSeason == 2020))
#nrow(orpi_replace)

# Remove original rows in dat and replace with those from new file with correct dates
#orpi_replace$Highlight <- as.character(orpi_replace$Highlight)
#dat <- filter(dat, !(str_detect(ImgPath, "ORPI_102_54W") & FieldSeason == 2020))
#dat <- bind_rows(dat, orpi_replace)

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
  dat$datetime <- parse_date_time(dat$ImageDate, c("%m/%d/%Y %I:%M:%S %p", 
                                                   "%m/%d/%Y"))
  # Create new date column
  dat$obsdate <- date(dat$datetime)
  # Create year variable (numeric)
  dat$yr <- year(dat$datetime)
  # Create month variable (numeric)
  dat$mon <- month(dat$datetime)
  # Create day-of-year variable (numeric)
  dat$yday <- yday(dat$datetime)
  # Convert time to a decimal value in [0,24)
  dat$obstime <- hour(dat$datetime) + minute(dat$datetime) / 60 + 
                 second(dat$datetime) / 3600
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
  
#-----------------------------------------------------------------------------------#
# Linking events file to observations 
# TODO: UPDATE THIS SECTION WHEN WE GET NEW EVENTS DATAFILE
#-----------------------------------------------------------------------------------# 
  
# Notes about sampling "events":
# Occasionally (esp at smaller parks), cameras were immediately re-deployed for continuous sampling
# Occasionally two cameras were deployed at the same location simultaneously
# Sometimes cameras left out for > 1 yr. Not sure how long they actually collected photos.
  
# Only keep necessary columns
events <- select(events, c(StdLocName, ProtocolVersion, DeployDate, RetrievalDate, 
                           CameraName, MountMethod, BatteryStatus, CameraSensitivity, 
                           DelaySec, ImagePer, TotalPics))

# Add column to identify Park
events$Park <- str_split_fixed(events$StdLocName, "_", 4)[,2]

# Format deployment date
  # Create new date-time column
  events$d_datetime <- parse_date_time(events$DeployDate, "%m/%d/%Y %H:%M")
  # Create new date column
  events$d_date <- date(events$d_datetime)
  # Create new year column
  events$d_yr <- year(events$d_datetime)
  # Remove DeployDate column
  events <- select(events, -DeployDate)

# Format retrieval date and time
  # Create new date-time colu
  events$r_datetime <- parse_date_time(events$RetrievalDate, "%m/%d/%Y %H:%M")
  # Create new date column
  events$r_date <- date(events$r_datetime)
  # Create new year column
  events$r_yr <- year(events$r_datetime)
  # Remove RetrievalDate column
  events <- select(events, -RetrievalDate)  
  
# Reorder columns and sort events dataframe by location and date
  events <- relocate(events, c(Park, d_date, r_date), .after = StdLocName)
  events <- arrange(events, StdLocName, d_datetime)

# Restrict events dataframe to 2016 forward
events <- filter(events, d_yr > 2015)

# Summarize sampling events by park and year (and compare to mammal observation data)
table(events$Park, events$d_yr)
table(dat$Park, dat$yr)
  # There are sampling events in CHIR in 2016 with no corresponding observations
  # Sampling methods were different that year, so probably ok not to track 
  # these down

# Calculate length of deployment, in days
events$duration <- as.double(difftime(as.POSIXct(events$r_datetime), 
                                      as.POSIXct(events$d_datetime), 
                                      units = 'days'))
  # Summarize/Visualize 
  summary(events$duration)
  # hist(events$duration, breaks = 25)

# View events with duration < 1 day
filter(events, duration < 1)
head(filter(events, Park == 'TONT'))
# Remove these events (all at TONT when cameras immediately re-deployed)
events <- filter(events, duration > 1)

# Check that cameras in our observations dataset all appear in the events dataset
obslocs <- sort(unique(dat$StdLocName))
obslocs[!obslocs %in% events$StdLocName]  # Yes, all appear in events
  
# Convert deployment/retrieval dates to integers, setting Jan 1 2016 equal to 1
events$d_day <- as.numeric(events$d_date) - as.numeric(as.Date("2015-12-31"))
events$r_day <- as.numeric(events$r_date) - as.numeric(as.Date("2015-12-31"))

# Create location by date matrix (1 indicates date was during a sampling event, 0 otherwise)
  
  # Create matrix with all 0s
  eventlocs <- sort(unique(events$StdLocName))
  event_mat <- matrix(0, nrow = length(eventlocs), ncol = max(events$r_day)) 
  rownames(event_mat) <- eventlocs
  
  # Replace 0s with 1s during each sampling event
  for (i in 1:nrow(event_mat)) {
    temp <- events[events$StdLocName == eventlocs[i], ]
    
    for (j in 1:nrow(temp)) {
      event_mat[i, temp$d_day[j]:temp$r_day[j]] <- 1
    }
  }
  # checks at random locations (only 1s during sampling event and 0s outside of event?):
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
  head(eventvec); head(events[,c(1,13:14)])
  
# Check that all photo dates were during listed sampling events
summary(dat$locdate %in% eventvec) # Yes, all photos during events

#-----------------------------------------------------------------------------------#
# Remove objects that are no longer needed
#-----------------------------------------------------------------------------------# 

rm(list = setdiff(ls(), c("dat", "locs", "events", "event_mat", "species")))
