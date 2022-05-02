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

# Using read_csv instead of read.csv because of formatting issues in the csv files 

# Observations
dat <- as.data.frame(read_csv('data/mammals/MAMMALS_ALL_2022-04-18.csv', trim_ws=TRUE, name_repair=make.names))

# Camera locations
locs <- as.data.frame(read_csv('data/mammals/SODN_Wildlife_Locations_XY_Revised_20220502.csv',trim_ws=TRUE,name_repair=make.names))

# Deployment schedule
events <- as.data.frame(read_csv('data/mammals/SODN_Wildlife_Events.csv',trim_ws=TRUE, name_repair=make.names))[ ,c(1,3:26)]

#-----------------------------------------------------------------------------------#
# Format and organize mammal observation data
#-----------------------------------------------------------------------------------#

# Remove empty or unnecessary columns:
dat <- select(dat,-c(StudyAreaID, UTM_E, UTM_N, UTMZone, FileName, ImgID, ImageNum, Highlight))

# Species
  # Create separate table with "species" information
  species <- count(dat,Species,Common.name,Species.code)
  # 35 different "species" -- includes Unknowns
  species[species$n %in% range(species$n),]
  # Number of observations per "species" ranges from 7 (Bighorn) to >12,000 (White-tailed deer) 
  # Remove Species and Common.name variables from observations frame (only need Species.code)
  dat <- select(dat,-c(Species,Common.name))

# Split up ImgPath name and append text with information about park, year, and camera location into new columns
summary(n.backslashes <- str_count(dat$ImgPath,"\\\\"))  # ImgPath always has 7 backslashes (ie, 8 character strings)
n.strings <- mean(n.backslashes) + 1
dat <- cbind(dat,str_split_fixed(dat$ImgPath,'\\\\',n.strings)[,5:7])
names(dat)[(ncol(dat)-2):ncol(dat)] <- c('Park','FY.filepath','Location')
  # checks:
  head(dat)
  count(dat,Park)         #7 parks (no GICL)
  count(dat,FY.filepath)  #7 FYs (16-22)
# Remove ImgPath column
dat <- select(dat,-ImgPath)

# Format date and time
  # Create new date-time column
  dat$dt <- parse_date_time(dat$ImageDate,'%m/%d/%Y %I:%M:%S %p')
  # Create new date column
  dat$obsdate <- date(dat$dt)
  # Create year variable (numeric)
  dat$yr <- year(dat$dt)
  # Create month variable (numeric)
  dat$mon <- month(dat$dt)
  # Create day-of-year variable (numeric)
  dat$yday <- yday(dat$dt)
  # Convert time to a decimal value in [0,24)
  dat$obstime <- hour(dat$dt) + minute(dat$dt)/60 + second(dat$dt)/3600
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
  
# Create a "loc.short" variable in locs and datlocs dataframes with simplified camera location name 
# (so they can be matched easily)
  
  # CAGR
  locs.cagr <- locs[locs$UnitCode=='CAGR',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]
  datlocs.cagr <- datlocs[datlocs$Park=='CAGR',]
  remove <- c('_SOLAR','CAGR_','V')
  locs.cagr$loc.short <- str_remove_all(locs.cagr$MarkerName,paste(remove,collapse='|'))
  datlocs.cagr$loc.short <- str_remove_all(datlocs.cagr$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.cagr$loc.short[!datlocs.cagr$loc.short %in% locs.cagr$loc.short]
  
  # CHIR
  locs.chir <- locs[locs$UnitCode=='CHIR',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]  
  datlocs.chir <- datlocs[datlocs$Park=='CHIR',]
  locs.chir$loc.short <- str_remove_all(locs.chir$MarkerName,'V')
  datlocs.chir$loc.short <- datlocs.chir$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.chir$loc.short[!datlocs.chir$loc.short %in% locs.chir$loc.short]
  
  # MOCC
  locs.mocc <- locs[locs$UnitCode=='MOCC',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]  
  datlocs.mocc <- datlocs[datlocs$Park=='MOCC',]
  remove <- c('WBC_','V')
  locs.mocc$loc.short <- str_remove_all(locs.mocc$MarkerName,paste(remove,collapse='|'))
  datlocs.mocc$loc.short <- str_remove_all(datlocs.mocc$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.mocc$loc.short[!datlocs.mocc$loc.short %in% locs.mocc$loc.short]
  
  # MOWE
  locs.mowe <- locs[locs$UnitCode=='MOWE',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')] 
  datlocs.mowe <- datlocs[datlocs$Park=='MOWE',]
  locs.mowe$loc.short <- str_remove_all(locs.mowe$MarkerName,'WBC_')
  datlocs.mowe$loc.short <- datlocs.mowe$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.mowe$loc.short[!datlocs.mowe$loc.short %in% locs.mowe$loc.short]
  
  # ORPI
  locs.orpi <- locs[locs$UnitCode=='ORPI',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')] 
  datlocs.orpi <- datlocs[datlocs$Park=='ORPI',]
  remove <- c('ORPI_','V')  
  locs.orpi$loc.short <- str_remove_all(locs.orpi$MarkerName,paste(remove,collapse='|'))
  datlocs.orpi$loc.short <- str_remove_all(datlocs.orpi$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.orpi$loc.short[!datlocs.orpi$loc.short %in% locs.orpi$loc.short]
  
  # SAGW
  locs.sagw <- locs[locs$UnitCode=='SAGW',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')]  
  datlocs.sagw <- datlocs[datlocs$Park=='SAGW',]
  remove <- c('SAGW ','W')   
  locs.sagw$loc.short <- str_remove_all(locs.sagw$MarkerName,paste(remove,collapse='|'))
  datlocs.sagw$loc.short <- str_remove_all(datlocs.sagw$Location,paste(remove,collapse='|'))
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.sagw$loc.short[!datlocs.sagw$loc.short %in% locs.sagw$loc.short]
  
  # TONT
  locs.tont <- locs[locs$UnitCode=='TONT',c('UnitCode','StdLocName','MarkerName','POINT_X','POINT_Y')] 
  datlocs.tont <- datlocs[datlocs$Park=='TONT',]
  locs.tont$loc.short <- str_remove_all(locs.tont$MarkerName,'V')
  datlocs.tont$loc.short <- datlocs.tont$Location
  # check: all cameras in observations file in the camera locations file? (if result is character(0), then yes)
  datlocs.tont$loc.short[!datlocs.tont$loc.short %in% locs.tont$loc.short]
  
# Combine datlocs dataframes and merge with dat
  datlocs <- rbind(datlocs.cagr, datlocs.chir, datlocs.mocc, datlocs.mowe, datlocs.orpi, datlocs.sagw, datlocs.tont)
  # Are there duplicates from the list of cameras in the observation file (datlocs)?
  xx <- count(datlocs,Park,Location)
  xx[xx$n>1,]
  datlocs[datlocs$Park=='ORPI' & datlocs$Location %in% c('101_007','102_004'),]
  # Yes, 101_007 camera at ORPI is associated with two LocationIDs in the observations file: 117 and 118
  # 102_004 camera at ORPI is associated with two LocationIDs in the observations file: 135 and 136
  locs.orpi[with(locs.orpi,order(loc.short,MarkerName)),]
  # Both these cameras only listed once in the camera locations file, so we'll remove duplicates from datlocs to match
  datlocs <- unique(datlocs[,c('Park','Location','loc.short')])
dat <- left_join(dat,datlocs,by=c('Park','Location'))
  
# Combine park-specific locs files and join spatial data to observations (dat)
locs2 <- rbind(locs.cagr, locs.chir, locs.mocc, locs.mowe, locs.orpi, locs.sagw, locs.tont)
locs2 <- rename(locs2,Park=UnitCode)
dat <- left_join(dat,locs2[,c('Park','loc.short','StdLocName','POINT_X','POINT_Y')],by=c('Park','loc.short'))
  # check:
  sum(is.na(dat$POINT_X)) #no NAs

#-----------------------------------------------------------------------------------#
# Linking events file to observations
#-----------------------------------------------------------------------------------#  
  