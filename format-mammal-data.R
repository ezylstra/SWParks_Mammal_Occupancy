########################################################################################################
# SODN -- Camera trap data, 2016-2022
# Data exploration

# ER Zylstra
# 2022-04-28
########################################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(readr)

# rm(list=ls())

#-----------------------------------------------------------------------------------#
# Import data
#-----------------------------------------------------------------------------------#

# Excluding a few unnecessary columns to avoid formatting issues

# Observations
dat <- read.csv("data/mammals/MAMMALS_ALL_2022-04-18.csv")
  # Change "." to "_" in column names
  colnames(dat) <- str_replace_all(colnames(dat), "[.]", "_")

# Camera locations
locs <- read.csv("data/mammals/SODN_Wildlife_Locations_XY_Revised_20220502.csv")[,2:9]

# Deployment schedule
events <- read.csv("data/mammals/SODN_Wildlife_Events.csv")[,3:26]

#-----------------------------------------------------------------------------------#
# Format and organize mammal observation data
#-----------------------------------------------------------------------------------#

# Remove empty or unnecessary columns:
dat <- select(dat,-c(StudyAreaID, UTM_E, UTM_N, UTMZone, FileName, ImgID, ImageNum, Highlight))

# Species
  # Create separate table with "species" information
  species <- count(dat,Species,Common_name,Species_code)
  # 35 different "species" -- includes Unknowns
  species[species$n %in% range(species$n),]
  # Number of observations per "species" ranges from 7 (Bighorn) to >12,000 (White-tailed deer) 
  # Remove Species and Common.name variables from observations frame (only need Species.code)
  dat <- select(dat,-c(Species,Common_name))

# Split up ImgPath name and append text with information about park, year, and camera location into new columns
summary(n.backslashes <- str_count(dat$ImgPath,"\\\\"))  # ImgPath always has 7 backslashes (ie, 8 character strings)
n.strings <- mean(n.backslashes) + 1
dat <- cbind(dat,str_split_fixed(dat$ImgPath,'\\\\',n.strings)[,5:7])
names(dat)[(ncol(dat)-2):ncol(dat)] <- c('Park','FY_filepath','Location')
  # checks:
  head(dat)
  count(dat,Park)         #7 parks (no GICL)
  count(dat,FY_filepath)  #7 FYs (16-22)
# Remove ImgPath column
dat <- select(dat,-ImgPath)

# Format date and time
  # Create new date-time column
  dat$datetime <- parse_date_time(dat$ImageDate,'%m/%d/%Y %I:%M:%S %p')
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
  dat <- select(dat,-ImageDate)
  
#-----------------------------------------------------------------------------------#
# Attach spatial data to observations
#-----------------------------------------------------------------------------------#

# Extract camera location "names" from observations dataframe
  datlocs <- unique(dat[,c('Park','Location','LocationID')])
  datlocs <- datlocs[with(datlocs,order(Park,Location)),]
  head(datlocs)
  count(datlocs,Park)
  
# Create a "loc_short" variable in locs and datlocs dataframes with simplified camera location name 
# (so they can be matched easily)
  
  # CAGR
  locs_cagr <- locs[locs$UnitCode=='CAGR',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]
  datlocs_cagr <- datlocs[datlocs$Park=='CAGR',]
  remove <- c('_SOLAR','CAGR_','V')
  locs_cagr$loc_short <- str_remove_all(locs_cagr$MarkerName,paste(remove,collapse='|'))
  datlocs_cagr$loc_short <- str_remove_all(datlocs_cagr$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_cagr$loc_short[!datlocs_cagr$loc_short %in% locs_cagr$loc_short]
  
  # CHIR
  locs_chir <- locs[locs$UnitCode=='CHIR',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]  
  datlocs_chir <- datlocs[datlocs$Park=='CHIR',]
  locs_chir$loc_short <- str_remove_all(locs_chir$MarkerName,'V')
  datlocs_chir$loc_short <- datlocs_chir$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_chir$loc_short[!datlocs_chir$loc_short %in% locs_chir$loc_short]
  
  # MOCC
  locs_mocc <- locs[locs$UnitCode=='MOCC',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]  
  datlocs_mocc <- datlocs[datlocs$Park=='MOCC',]
  remove <- c('WBC_','V')
  locs_mocc$loc_short <- str_remove_all(locs_mocc$MarkerName,paste(remove,collapse='|'))
  datlocs_mocc$loc_short <- str_remove_all(datlocs_mocc$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_mocc$loc_short[!datlocs_mocc$loc_short %in% locs_mocc$loc_short]
  
  # MOWE
  locs_mowe <- locs[locs$UnitCode=='MOWE',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')] 
  datlocs_mowe <- datlocs[datlocs$Park=='MOWE',]
  locs_mowe$loc_short <- str_remove_all(locs_mowe$MarkerName,'WBC_')
  datlocs_mowe$loc_short <- datlocs_mowe$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_mowe$loc_short[!datlocs_mowe$loc_short %in% locs_mowe$loc_short]
  
  # ORPI
  locs_orpi <- locs[locs$UnitCode=='ORPI',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')] 
  datlocs_orpi <- datlocs[datlocs$Park=='ORPI',]
  remove <- c('ORPI_','V')  
  locs_orpi$loc_short <- str_remove_all(locs_orpi$MarkerName,paste(remove,collapse='|'))
  datlocs_orpi$loc_short <- str_remove_all(datlocs_orpi$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_orpi$loc_short[!datlocs_orpi$loc_short %in% locs_orpi$loc_short]
  
  # SAGW
  locs_sagw <- locs[locs$UnitCode=='SAGW',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]  
  datlocs_sagw <- datlocs[datlocs$Park=='SAGW',]
  remove <- c('SAGW ','W')   
  locs_sagw$loc_short <- str_remove_all(locs_sagw$MarkerName,paste(remove,collapse='|'))
  datlocs_sagw$loc_short <- str_remove_all(datlocs_sagw$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_sagw$loc_short[!datlocs_sagw$loc_short %in% locs_sagw$loc_short]
  
  # TONT
  locs_tont <- locs[locs$UnitCode=='TONT',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')] 
  datlocs_tont <- datlocs[datlocs$Park=='TONT',]
  locs_tont$loc_short <- str_remove_all(locs_tont$MarkerName,'V')
  datlocs_tont$loc_short <- datlocs_tont$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs_tont$loc_short[!datlocs_tont$loc_short %in% locs_tont$loc_short]
  
# Combine datlocs dataframes and merge with dat
  datlocs <- rbind(datlocs_cagr, datlocs_chir, datlocs_mocc, datlocs_mowe, datlocs_orpi, datlocs_sagw, datlocs_tont)
  # Are there duplicates from the list of cameras in the observation file (datlocs)?
  xx <- count(datlocs,Park,Location)
  xx[xx$n>1,]
  datlocs[datlocs$Park=='ORPI' & datlocs$Location %in% c('101_007','102_004'),]
  # Yes, 101_007 camera at ORPI is associated with two LocationIDs in the observations file: 117 and 118
  # 102_004 camera at ORPI is associated with two LocationIDs in the observations file: 135 and 136
  locs_orpi[with(locs_orpi,order(loc_short,MarkerName)),]
  # Both these cameras only listed once in the camera locations file, so we'll remove duplicates from datlocs to match
  datlocs <- unique(datlocs[,c('Park','Location','loc_short')])
dat <- left_join(dat,datlocs,by=c('Park','Location'))
  
# Combine park-specific locs files and join spatial data to observations (dat)
locs2 <- rbind(locs_cagr, locs_chir, locs_mocc, locs_mowe, locs_orpi, locs_sagw, locs_tont)
locs2 <- rename(locs2,Park=UnitCode)
dat <- left_join(dat,locs2[,c('Park','loc_short','StdLocName','POINT_X','POINT_Y')],by=c('Park','loc_short'))
  # check:
  sum(is.na(dat$POINT_X)) #no NAs

#-----------------------------------------------------------------------------------#
# Linking events file to observations
#-----------------------------------------------------------------------------------#  
  