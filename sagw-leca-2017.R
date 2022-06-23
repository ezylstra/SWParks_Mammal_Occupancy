################################################################################
# Single-season occupancy analysis
# SAGW, LECA (black-tailed jackrabbit), 2017

# ER Zylstra
# 2022-06-15
################################################################################

# Note: Will probably want to create a template script for single-season analyses
# This will serve as a first case study

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(jagsUI)

# rm(list = ls())

# Load photo, location, events, species data 
source("format-mammal-data.R")

  # dat = information about each photo (date, time, species, location)
  # events = information about each camera deployment (dates, location, duration)
  # event_mat = camera location x day matrix with 1/0 indicating whether camera
  #             was deployed or not
  # locs = information about each camera location (park, lat/long, name)
  # species = table with species observed (species code, common name, # of obs)

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

# Identify park, species, and year of interest
park <- "SAGW"
species <- "LECA"
year <- 2017

# Use day number (day 1 = 01 Jan 2016) for column names in event_mat 
colnames(event_mat) <- 1:ncol(event_mat)

#-------------------------------------------------------------------------------#
# Create detection histories for selected park, species, and year
#-------------------------------------------------------------------------------#

# Extract sampling occasion info for selected park and year
occasions <- occasions %>%
  filter(Park == park & yr == year)

# Convert occasion start/end dates to day numbers
occasions$start_day <- 
  as.numeric(date(occasions$start)) - as.numeric(as.Date("2015-12-31"))
occasions$end_day <- 
  as.numeric(date(occasions$end)) - as.numeric(as.Date("2015-12-31"))

# Create a list of days included in sampling occasions
occ_days <- NULL
for (i in 1:nrow(occasions)) {
  occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
}

# Extract photo observations for park, species and year
# Retain a maximum of one observation per day at each location
obs <- dat %>% 
  filter(Park == park & Species_code == species & yr == year) %>%
  select(StdLocName, obsdate, yr, o_day) %>%
  arrange(StdLocName, obsdate) %>%
  distinct

# Extract information about camera locations in selected park
locs_park <- locs %>%
  filter(UnitCode == park) %>%
  select(StdLocName, POINT_X, POINT_Y) %>%
  rename(loc = StdLocName, long = POINT_X, lat = POINT_Y)

# Extract rows from events matrix that correspond to locations in selected park
event_mat <- event_mat[rownames(event_mat) %in% locs_park$loc,]

# Extract columns from events matrix that correspond to sampling occasions
event_mat <- event_mat[,colnames(event_mat) %in% occ_days]

# Convert the events matrix into daily detection histories (ddh)
ddh <- event_mat

  # Change 0s to NA (NA = camera wasn't deployed)
  ddh[ddh == 0] <- NA
  
  # Change 1s to 0 (0 indicates that the species wasn't detected)
  ddh[ddh == 1] <- 0
  
  # Replace 0s with 1s when the species was detected
  for (i in 1:nrow(obs)) {
    ddh[rownames(ddh) == obs$StdLocName[i], 
        colnames(ddh) == as.character(obs$o_day[i])] <- 1
  }
  # checks:
  sum(ddh == 1, na.rm = TRUE)
  sum(obs$o_day %in% occ_days)

# Create a function to aggregate daily detection data during each occasion
  # NA if camera wasn't operational throughout entire occasion (all values = NA)
  # 1 if species was detected one or more times (even if there are NAs)
  # 0 if species was never detected
  paNA <- function(x) {
    if (sum(is.na(x)) == length(x)) {NA} else 
      if (sum(x, na.rm = TRUE) == 0) {0} else {1} 
  }
  
# Create a function to calculate the proportion of a sampling period a camera
# was operational
  propNA <- function(x) {
    (occasions$duration[1] - sum(is.na(x))) / occasions$duration[1]
  }

# Summarize detection data (dh) and effort during each occasion 
dh <- effort <- matrix(NA, 
                       nrow = nrow(ddh), 
                       ncol = ncol(ddh) / occasions$duration[1],
                       dimnames = list(rownames(ddh), NULL))

for (i in 1:ncol(dh)) {
  multiday <- ddh[,colnames(ddh) %in% occasions$start_day[i]:occasions$end_day[i]]
  dh[,i] <- apply(multiday, 1, paNA)
  effort[,i] <- apply(multiday, 1, propNA)
}

#-------------------------------------------------------------------------------#
# Spatial covariates
#-------------------------------------------------------------------------------#
# Will expand this section as more covariates become available

spatial_covs <- locs_park 

# Ensure the order is the same as what's in the detection history matrix
spatial_covs <- spatial_covs[match(rownames(dh), spatial_covs$loc),]

# Identify continuous covariates
covs_cont <- c("long", "lat")

# Scale continuous covariates by mean, SD
for (i in covs_cont) {
  meani <- mean(spatial_covs[,i])
  sdi <- sd(spatial_covs[,i])
  spatial_covs[,paste0(i, "_z")] <- (spatial_covs[,i] - meani) / sdi
}

#-------------------------------------------------------------------------------#
# Package things up for JAGS
#-------------------------------------------------------------------------------#

# Function to create a matrix with information about known latent states, z[i]
# JAGS won't try to estimate z when site is known to be occupied
known_state_occ <- function(dh){
  state <- apply(dh, 1, function(x) ifelse(sum(is.na(x)) == length(x), NA, max(x, na.rm=T)))
  state[state==0] <- NA
  return(state)
}

# Function to create initial values for unknown latent states, z[i]
inits_state_occ <- function(dh){
  state <- apply(dh, 1, function(x) ifelse(sum(is.na(x)) == length(x), NA, max(x, na.rm=T)))
  # Initial value of 1 whenever occupancy state is unknown
  state[state==1] <- 2
  state[is.na(state) | state==0] <- 1
  state[state==2] <- NA
  return(state)
}  

# Bundle data for JAGS
jags_data <- list(y = dh,
                  n_sites = nrow(dh),
                  n_surveys = ncol(dh),
                  lat = spatial_covs$lat_z,
                  effort = effort,
                  z = known_state_occ(dh))

# List of parameters to monitor
params <- c("mean_psi", "beta0", "beta1", 
            "mean_p", "alpha0", "alpha1")

# Initial values
inits <- function(){list(mean_psi = runif(1, 0, 1),
                         beta1 = runif(1, -2, 2),
                         mean_p = runif(1, 0, 1),
                         alpha1 = runif(1, -2, 2),
                         z = inits_state_occ(dh))}

#-------------------------------------------------------------------------------#
# Run model in JAGS
#-------------------------------------------------------------------------------#

nc <- 3      # Number of chains
na <- 2000   # Number of iterations to run in the adaptive phase
nb <- 5000   # Number of iterations to discard (burn-in)
ni <- 20000  # Number of iterations per chain (including burn-in)
nt <- 10     # Thinning rate

out <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "JAGS_SingleSeasonWithCovs.txt",
            n.chains = nc,
            n.adapt = na,
            n.burnin = nb,
            n.iter = ni,
            n.thin = nt,
            parallel = TRUE)

print(out)
plot(out)

# Mean detection probability not estimated well
  # Model parameterization?
  # Priors?
  # Is this the time to learn nimble?
  # Different, better covariates (latitude and effort just placeholders, really)

# For comparison, should run same model in unmarked
