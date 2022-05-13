########################################################################################################
# SODN -- Camera trap data, 2016-2022
# Data exploration

# ER Zylstra
# 2022-04-28
########################################################################################################

library(dplyr)
library(lubridate)
library(stringr)

# rm(list=ls())

#-----------------------------------------------------------------------------------#
# Import data
#-----------------------------------------------------------------------------------#

# Excluding a few unnecessary columns to avoid formatting issues

# Observations
dat <- read.csv("data/mammals/MAMMALS_ALL_2022-04-18.csv")
  # Change "." to "_" in column names
  colnames(dat) <- str_replace_all(colnames(dat), "[.]", "_")
  
orpi_replace <- read.csv("data/mammals/ORPI_2020_102_54W_adjustedDates.csv")[,c(1:2,4:14,16:19)]
  colnames(orpi_replace)[c(1:13,16:17)] <- colnames(dat)[c(1:13,15:16)]
  
# Camera locations
locs <- read.csv("data/mammals/SODN_Wildlife_Locations_XY_Revised_20220502.csv")[,2:9]

# Deployment schedule
events <- read.csv("data/mammals/SODN_Wildlife_Events_Revised_20220512.csv")[,3:26]

#-----------------------------------------------------------------------------------#
# Replace data from one ORPI camera in 2020 with data that has corrected date 
#-----------------------------------------------------------------------------------#
# Remove datasheet photo
orpi_replace <- filter(orpi_replace, Common_name != 'Datasheet')

# Fix a few species names/codes
orpi_replace <- mutate(orpi_replace, 
                       Species = ifelse(Common_name == "Unknown fox", 
                                        "Caninae sp.",
                                        paste(Genus, Species, sep = " ")),
                       Species_code = ifelse(Common_name == "Unknown fox",
                                             "UNFO", 
                                             Species_code))
orpi_replace <- select(orpi_replace, !Genus)

# Check that number of 2020 observations at ORPI_102_54W in dat is the same as replacement file
nrow(filter(dat, str_detect(ImgPath, "ORPI_102_54W") & FieldSeason == 2020))
nrow(orpi_replace)

# Remove original rows in dat and replace with those from new file with correct dates
orpi_replace$Highlight <- as.character(orpi_replace$Highlight)
dat <- filter(dat, !(str_detect(ImgPath, "ORPI_102_54W") & FieldSeason == 2020))
dat <- bind_rows(dat, orpi_replace)

#-----------------------------------------------------------------------------------#
# Format and organize mammal observation data
#-----------------------------------------------------------------------------------#

# Remove empty or unnecessary columns:
dat <- select(dat,-c(StudyAreaID, UTM_E, UTM_N, UTMZone, FileName, ImgID, ImageNum, Highlight))

# Species
  # Create separate table with "species" information
  species <- count(dat, Species, Common_name, Species_code)
  # 35 different "species" -- includes Unknowns
  species[species$n %in% range(species$n),]
  # Number of observations per "species" ranges from 7 (Bighorn) to >12,000 (WT deer) 
  # Remove Species and Common_name variables from observations frame (have Species_code)
  dat <- select(dat, -c(Species, Common_name))

# Split up ImgPath and append information about park, year, and location into new columns
summary(n.backslashes <- str_count(dat$ImgPath,"\\\\"))  
  # ImgPath always has 7 backslashes (ie, 8 character strings)
n.strings <- mean(n.backslashes) + 1
dat <- cbind(dat,str_split_fixed(dat$ImgPath, "\\\\", n.strings)[,5:7])
names(dat)[(ncol(dat) - 2):ncol(dat)] <- c("Park", "FY_filepath", "Location")
  # checks:
  head(dat)
  count(dat, Park)         #7 parks (no GICL)
  count(dat, FY_filepath)  #7 FYs (16-22)
# Remove ImgPath column
dat <- select(dat, -ImgPath)

# Format date and time
  # Create new date-time column
  dat$datetime <- parse_date_time(dat$ImageDate,c("%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y"))
  # Create new date column
  dat$obsdate <- date(dat$datetime)
  # Create year variable (numeric)
  dat$yr <- year(dat$datetime)
  # Create month variable (numeric)
  dat$mon <- month(dat$datetime)
  # Create day-of-year variable (numeric)
  dat$yday <- yday(dat$datetime)
  # Convert time to a decimal value in [0,24)
  dat$obstime <- hour(dat$datetime) + minute(dat$datetime)/60 + second(dat$datetime)/3600
  # Remove ImageDate column
  dat <- select(dat, -ImageDate)
  
#-----------------------------------------------------------------------------------#
# Attach spatial data to observations
#-----------------------------------------------------------------------------------#

# Extract camera location "names" from observations dataframe
  datlocs <- unique(dat[,c("Park", "Location", "LocationID")])
  datlocs <- arrange(datlocs, Park, Location)
  head(datlocs)
  count(datlocs, Park)
  
# Create a "loc_short" variable in locs and datlocs dataframes with simplified camera location name 
# (so they can be matched easily)
  
  # CAGR
  locs_cagr <- locs[locs$UnitCode == "CAGR", c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")]
  datlocs_cagr <- datlocs[datlocs$Park == "CAGR", ]
  remove <- c("_SOLAR", "CAGR_", "V")
  locs_cagr$loc_short <- str_remove_all(locs_cagr$MarkerName, paste(remove, collapse="|"))
  datlocs_cagr$loc_short <- str_remove_all(datlocs_cagr$Location, paste(remove, collapse="|"))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_cagr$loc_short[!datlocs_cagr$loc_short %in% locs_cagr$loc_short]
  
  # CHIR
  locs_chir <- locs[locs$UnitCode == "CHIR",c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")]  
  datlocs_chir <- datlocs[datlocs$Park == "CHIR",]
  locs_chir$loc_short <- str_remove_all(locs_chir$MarkerName, "V")
  datlocs_chir$loc_short <- datlocs_chir$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_chir$loc_short[!datlocs_chir$loc_short %in% locs_chir$loc_short]
  
  # MOCC
  locs_mocc <- locs[locs$UnitCode == "MOCC",c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")]  
  datlocs_mocc <- datlocs[datlocs$Park == "MOCC",]
  remove <- c("WBC_", "V")
  locs_mocc$loc_short <- str_remove_all(locs_mocc$MarkerName, paste(remove, collapse="|"))
  datlocs_mocc$loc_short <- str_remove_all(datlocs_mocc$Location, paste(remove, collapse="|"))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_mocc$loc_short[!datlocs_mocc$loc_short %in% locs_mocc$loc_short]
  
  # MOWE
  locs_mowe <- locs[locs$UnitCode == "MOWE",c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")] 
  datlocs_mowe <- datlocs[datlocs$Park == "MOWE",]
  locs_mowe$loc_short <- str_remove_all(locs_mowe$MarkerName, "WBC_")
  datlocs_mowe$loc_short <- datlocs_mowe$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_mowe$loc_short[!datlocs_mowe$loc_short %in% locs_mowe$loc_short]
  
  # ORPI
  locs_orpi <- locs[locs$UnitCode == "ORPI",c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")] 
  datlocs_orpi <- datlocs[datlocs$Park == "ORPI",]
  remove <- c("ORPI_", "V")  
  locs_orpi$loc_short <- str_remove_all(locs_orpi$MarkerName, paste(remove, collapse="|"))
  datlocs_orpi$loc_short <- str_remove_all(datlocs_orpi$Location, paste(remove, collapse="|"))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_orpi$loc_short[!datlocs_orpi$loc_short %in% locs_orpi$loc_short]
  
  # SAGW
  locs_sagw <- locs[locs$UnitCode == "SAGW",c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")]  
  datlocs_sagw <- datlocs[datlocs$Park == "SAGW",]
  remove <- c("SAGW ", "W")   
  locs_sagw$loc_short <- str_remove_all(locs_sagw$MarkerName, paste(remove, collapse="|"))
  datlocs_sagw$loc_short <- str_remove_all(datlocs_sagw$Location, paste(remove, collapse="|"))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_sagw$loc_short[!datlocs_sagw$loc_short %in% locs_sagw$loc_short]
  
  # TONT
  locs_tont <- locs[locs$UnitCode == "TONT",c("UnitCode", "StdLocName", "MarkerName", "POINT_X", "POINT_Y")] 
  datlocs_tont <- datlocs[datlocs$Park == "TONT",]
  locs_tont$loc_short <- str_remove_all(locs_tont$MarkerName, "V")
  datlocs_tont$loc_short <- datlocs_tont$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_tont$loc_short[!datlocs_tont$loc_short %in% locs_tont$loc_short]
  
# Combine datlocs dataframes and merge with dat
  datlocs <- rbind(datlocs_cagr, datlocs_chir, datlocs_mocc, datlocs_mowe, 
                   datlocs_orpi, datlocs_sagw, datlocs_tont)
  # Are there duplicates from the list of cameras in the observation file (datlocs)?
  xx <- count(datlocs, Park, Location)
  xx[xx$n > 1, ]
  datlocs[datlocs$Park == "ORPI" & datlocs$Location %in% c("101_007", "102_004"), ]
  # Yes, ORIPI 101_007 is associated with two LocationIDs in the observations file: 117 and 118
  # ORPI 102_004 is associated with two LocationIDs in the observations file: 135 and 136
  locs_orpi[with(locs_orpi, order(loc_short, MarkerName)), ]
  # Both these cameras only listed once in the camera locations file, 
  # so we'll remove duplicates from datlocs to match
  datlocs <- unique(datlocs[ ,c("Park", "Location", "loc_short")])
dat <- left_join(dat, datlocs, by=c("Park", "Location"))
  
# Combine park-specific locs files and join spatial data to observations (dat)
locs2 <- rbind(locs_cagr, locs_chir, locs_mocc, locs_mowe, locs_orpi, locs_sagw, locs_tont)
locs2 <- rename(locs2, Park = UnitCode)
dat <- left_join(dat, locs2[,c("Park", "loc_short", "StdLocName", "POINT_X", "POINT_Y")], 
                 by = c("Park", "loc_short"))
  # check:
  sum(is.na(dat$POINT_X)) #no NAs

#-----------------------------------------------------------------------------------#
# Linking events file to observations
#-----------------------------------------------------------------------------------# 
  
# Notes about sampling "events":
# Occasionally (esp at smaller parks), cameras were immediately re-deployed for continuous sampling
# Occasionally two cameras were deployed at the same location simultaneously
# Sometimes cameras left out for > 1 yr. Not sure how long they actually collected photos.
  
# Only keep necessary columns
events <- select(events, c(StdLocName, ProtocolVersion, DeployDate, RetrievalDate, CameraName, TotalPics))

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
  # SAMPLING EVENTS WITH NO MAMMAL OBSERVATIONS:
    # CHIR in 2016
 
# Calculate length of deployment, in days
events$duration <- as.double(difftime(as.POSIXct(events$r_datetime), as.POSIXct(events$d_datetime), units='days'))
  # Summarize/Visualize 
  summary(events$duration)
  hist(events$duration, breaks = 25)

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
dat$locdate <- paste(dat$StdLocName, dat$o_day, sep="_")

# Create character strings with location and date for day during each sampling event
eventvec <- as.character(vector())
for (i in 1:length(eventlocs)) {
  eventvec <- append(eventvec, paste(eventlocs[i], which(event_mat[i,] == 1), sep="_"))
}
  # check:
  head(eventvec); head(events[,c(1,13:14)])
  
# Check that all photo dates were during listed sampling events
summary(dat$locdate %in% eventvec) # Yes, all photos during events
