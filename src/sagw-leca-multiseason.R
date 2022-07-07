################################################################################
# Multi-season occupancy analysis
# SAGW, LECA (black-tailed jackrabbit)

# ER Zylstra
# Updated 2022-06-30
################################################################################

# Note: Will probably want to create a template script for this type of analysis
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

# Use day number (day 1 = 01 Jan 2016) for column names in event_mat 
colnames(event_mat) <- 1:ncol(event_mat)

#-------------------------------------------------------------------------------#
# Create detection histories for selected park, species
#-------------------------------------------------------------------------------#

# Extract sampling occasion info for selected park
occasions <- occasions %>%
  filter(Park == park) %>%
  arrange(yr, occasion)

# Add occasion ID
occasions$yr_occ <- paste0(occasions$yr, "_", occasions$occasion)

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

# Extract photo observations for park, species
# Retain a maximum of one observation per day at each location
obs <- dat %>% 
  filter(Park == park & Species_code == species) %>%
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

# In contrast to our single-season model, we'll put detection/survey data in 
# long form. This makes it easier to create universal model scripts that can
# be used with different combinations of temporal and spatial covariates. This 
# may also result in shorter run times. 

# Basically, instead of having site * occasion * season arrays, we'll create 
# long vectors with detection and effort data, along with accompanying variables 
# to indicate site, location, and season associated with the given observations.

# Convert detection data into long form
dh_df <- as.data.frame(dh)
colnames(dh_df) <- sort(occasions$yr_occ)
dh_df$loc <- row.names(dh_df)
dh_long <- dh_df %>%
  pivot_longer(!loc,
               names_to = "occ",
               values_to = "det") %>%
  arrange(loc, occ) %>%
  as.data.frame

# Convert effort data into long form
eff_df <- as.data.frame(effort)
colnames(eff_df) <- sort(occasions$yr_occ)
eff_df$loc <- row.names(eff_df)
eff_long <- eff_df %>%
  pivot_longer(!loc,
               names_to = "occ",
               values_to = "effort") %>%
  arrange(loc, occ) %>%
  as.data.frame

# Merge detection and effort data
surveys <- left_join(dh_long, eff_long)

# Add year and trend columns
surveys$yr <- as.numeric(str_sub(surveys$occ, 1,4))
surveys$trend <- surveys$yr - min(surveys$yr)

# Add site, season, occasion indicators
surveys$season_index <- surveys$trend + 1
surveys$occ_index <- as.numeric(str_sub(surveys$occ,6,nchar(surveys$occ)))
for (i in 1:nrow(surveys)) {
  surveys$site_index[i] <- which(rownames(dh) == surveys$loc[i])
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

# Attach spatial covariates to surveys dataframe so we can use them as 
# covariates for detection
surveys <- left_join(surveys, spatial_covs[,c("loc", "long_z", "lat_z")])

#-------------------------------------------------------------------------------#
# Package things up for JAGS
#-------------------------------------------------------------------------------#

# z (latent occupancy for each site & season) is what we're interested in
# In JAGS, z will be stored in a matrix (n_sites[i] * n_seasons[t])

# In order to create initial values for z, we'll need to summarize detections 
# over occasions at each site in each season. To do this:
  # Create an array with detection data (row = site, col = occ, slice = season)
  # Then use an apply function, summarizing over columns

n_sites <- max(surveys$site_index) # 60
n_seasons <- max(surveys$season_index) # 6
max_n_occasions <- max(surveys$occ_index) # 6
season_occ <- data.frame(season_index = 1:n_seasons)
season_occ$yr <- surveys$yr[match(season_occ$season_index, surveys$season_index)]
for (i in 1:n_seasons) {
  season_occ$n_occasions[i] <- 
    ifelse(sum(surveys$season_index == i) == 0, NA, 
           max(surveys$occ_index[surveys$season_index == i]))
}

# Create an array with detection data (row = site, col = occ, slice = season)
y_array <- array(NA, dim = c(n_sites, max_n_occasions, n_seasons))
for (i in 1:n_seasons) {
  if (is.na(season_occ$n_occasions[i])) {next}
  y_array[,1:season_occ$n_occasions[i],i] <-
    as.matrix(dh_df[,grep(season_occ$yr[i], colnames(dh_df))])
}

# Function to create a matrix with information about known latent states, z[i,t]
# JAGS won't try to estimate z when site is known to be occupied
known_state_occ <- function(y_array){
  state <- apply(y_array, 
                 c(1,3), 
                 function(x) ifelse(sum(is.na(x)) == length(x), 
                                    NA, max(x, na.rm=T)))
  state[state==0] <- NA
  return(state)
}
# check:
known_state_occ(y_array)

# Function to create initial values for unknown latent states, z[i, t]
inits_state_occ <- function(y_array){
  state <- apply(y_array, 
                 c(1,3), 
                 function(x) ifelse(sum(is.na(x)) == length(x), 
                                    NA, 
                                    max(x, na.rm=T)))
  # Initial value of 1 whenever occupancy state is unknown
  state[state==1] <- 2
  state[is.na(state) | state==0] <- 1
  state[state==2] <- NA
  return(state)
}  
# check:
inits_state_occ(y_array)

# Create object with covariates for occupancy in first season
# Here, using latitude and latitude^2
cov_psi <- spatial_covs %>%
  select(lat_z) %>%
  mutate(lat_z2 = lat_z * lat_z) %>%
  as.matrix
n_cov_psi <- ncol(cov_psi)

# Create object with covariates for detection
# Here, using latitude and effort
cov_p <- surveys %>%
  select(lat_z, effort) %>%
  as.matrix
n_cov_p <- ncol(cov_p)

# Bundle data for JAGS
jags_data <- list(y = surveys$det,
                  n_sites = n_sites,
                  n_seasons = n_seasons,
                  n_obs = nrow(surveys),
                  cov_psi = cov_psi,
                  n_cov_psi = ncol(cov_psi),
                  cov_p = cov_p,
                  n_cov_p = ncol(cov_p),
                  lat = spatial_covs$lat_z,
                  long = spatial_covs$long_z,
                  site = surveys$site_index,
                  season = surveys$season_index,
                  z = known_state_occ(y_array))

# List of parameters to monitor
params <- c("mean_psi", "beta_psi0", "beta_psi",
            "mean_p", "beta_p0", "beta_p",
            "beta_eps0", "beta_eps1",
            "beta_gam0", "beta_gam1")

# Initial values
inits <- function(){list(mean_psi = runif(1, 0, 1),
                         beta_psi = runif(n_cov_psi, -2, -2),
                         mean_p = runif(1, 0, 1),
                         beta_p = runif(n_cov_p, -2, 2),
                         beta_eps0 = runif(1, -1, 1),
                         beta_eps1 = runif(1, -2 ,2),
                         beta_gam0 = runif(1, -1, 1),
                         beta_gam1 = runif(1, -2, 2),
                         z = inits_state_occ(y_array))}

#-------------------------------------------------------------------------------#
# Run model in JAGS
#-------------------------------------------------------------------------------#

nc <- 3      # Number of chains
na <- 3000   # Number of iterations to run in the adaptive phase
nb <- 10000  # Number of iterations to discard (burn-in)
ni <- 20000  # Number of iterations per chain (including burn-in)
nt <- 10     # Thinning rate

out <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "JAGS_MultiSeasonWithCovs.txt",
            n.chains = nc,
            n.adapt = na,
            n.burnin = nb,
            n.iter = ni,
            n.thin = nt,
            parallel = TRUE)

print(out)

library(MCMCvis)

# Trace and density plots 
MCMCtrace(out,pdf = FALSE)

#-------------------------------------------------------------------------------#
# For comparison, running the same model in unmarked
#-------------------------------------------------------------------------------#

library(unmarked)
dh_full <- cbind(dh[,1:5], rep(NA, n_sites),
                 dh[,6:9], rep(NA, n_sites), rep(NA, n_sites),
                 matrix(NA, ncol = 6, nrow = n_sites),
                 dh[,10:15],
                 dh[,16:20], rep(NA, n_sites),
                 dh[,21:25], rep(NA, n_sites))
eff_full <- cbind(effort[,1:5], rep(NA, n_sites),
                  effort[,6:9], rep(NA, n_sites), rep(NA, n_sites),
                  matrix(NA, ncol = 6, nrow = n_sites),
                  effort[,10:15],
                  effort[,16:20], rep(NA, n_sites),
                  effort[,21:25], rep(NA, n_sites))

umf <- unmarkedMultFrame(y = dh_full,
                         siteCovs = data.frame(lat = spatial_covs$lat_z,
                                               lat2 = spatial_covs$lat_z^2,
                                               long = spatial_covs$long_z),
                         obsCovs = list(effort = eff_full),
                         numPrimary = max(surveys$season_index))
summary(umf)
summary(m <- colext(~ lat + lat2, ~ long, ~ lat, ~ lat + effort, data = umf)) 
# Estimates are pretty similar, but again, detection estimated poorly
