################################################################################
# Test script for single-season occupancy analysis using spOccupancy

# ER Zylstra
# Updated 2023-01-25
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(spOccupancy)

rm(list = ls())

#------------------------------------------------------------------------------#
# Load all photo, location, events, species data 
#------------------------------------------------------------------------------#

source("src/photo-data/format-mammal-data.R")

# dat = information about each photo (date, time, species, location)
# events = information about each camera deployment (dates, location, duration)
# event_mat = camera location x day matrix with 1/0 indicating whether camera
#             was deployed or not
# locs = information about each camera location (park, lat/long, name)
# species = table with species observed (species code, common name, # of obs)

#------------------------------------------------------------------------------#
# Specify model parameters (objects in all caps)
#------------------------------------------------------------------------------#

# Select park of interest ("CHIR", "ORPI", or "SAGW")
PARK <- "SAGW"

# Select year of interest
YEAR <- 2017

# Figure out which species had a sufficient number of detections
detects <- read.csv("output/species-detections-byparkyr.csv", header = TRUE)
detects <- detects %>%
  filter(Park == PARK & yr == YEAR) %>%
  arrange(desc(propdetect))
# Use detection rate of 5% (n = camera location * sampling occasion) as cutoff
detects %>% 
  filter(propdetect >= 0.05) %>%
  select(c(spp, propdetect))

# Select species of interest from list above
SPECIES <- "LECA"

# Prep detection and covariate data (eventually this will be sourced out)
# source("src/single-season-models/spOccupancy-prep-data.R")

#------------------------------------------------------------------------------#
# Create detection histories for selected park, species, and year
#------------------------------------------------------------------------------#

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

# Extract sampling occasion info for selected park and year
occasions <- occasions %>%
  filter(Park == PARK) %>%
  filter(yr %in% YEAR) %>%
  arrange(yr, occasion)

# Create a list of days included in sampling occasions
occ_days <- NULL
for (i in 1:nrow(occasions)) {
  occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
}

# Extract photo observations for park, species, year
# Retain a maximum of one observation per day at each location
obs <- dat %>% 
  filter(Park == PARK & Species_code == SPECIES & yr == YEAR) %>%
  select(StdLocName, obsdate, yr, o_day) %>%
  arrange(StdLocName, obsdate) %>%
  distinct

# Extract information about camera locations in selected park
events_park <- events %>%
  filter(Park == PARK) %>%
  filter(d_yr == YEAR)
locs_park <- locs %>%
  filter(StdLocName %in% events_park$StdLocName) %>%
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
# sum(ddh == 1, na.rm = TRUE)
# sum(obs$o_day %in% occ_days)

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
                       ncol = nrow(occasions),
                       dimnames = list(rownames(ddh), NULL))

for (i in 1:ncol(dh)) {
  multiday <- ddh[,colnames(ddh) %in% occasions$start_day[i]:occasions$end_day[i]]
  dh[,i] <- apply(multiday, 1, paNA)
  effort[,i] <- apply(multiday, 1, propNA)
}

# Create standardized version of effort variable
effort_mn <- mean(effort)
effort_sd <- sd(effort)
effort_z <- (effort - effort_mn)/effort_sd

#------------------------------------------------------------------------------#
# Extract information about experience of deployment personnel
#------------------------------------------------------------------------------#

# Add a column for each occasion to events_park, indicating whether the 
# occasion was associated with that deployment event or not (1/0) 
occ_matrix <- matrix(NA, nrow = nrow(events_park), ncol = nrow(occasions))
colnames(occ_matrix) <- occasions$yr_occ
for (i in 1:nrow(occ_matrix)) {
  for (j in 1:ncol(occ_matrix)) {
    occ_matrix[i,j] <- 1 * any(occasions$start_day[j]:occasions$end_day[j] %in% 
                                 events_park$d_day[i]:events_park$r_day[i])
  }
}
events_park <- cbind(events_park, occ_matrix)
# Remove rows in the dataframe that aren't associated with any sampling occasion
events_park <- events_park %>%
  mutate(n_occ = rowSums(select(., paste0(YEAR, "_", 1:ncol(occ_matrix))))) %>%
  filter(n_occ > 0) %>%
  arrange(StdLocName, d_date)

# Check if a sampling occasion at a given camera location spanned two 
# deployments (ie, a camera was immediately redeployed during a sampling occ)
# If so, use the deployment experience value from the 2nd deployment
redeploys <- events_park %>%
  group_by(StdLocName) %>%
  summarize(across(starts_with(paste0(YEAR, "_")), sum)) %>%
  as.data.frame()

for (i in 1:nrow(redeploys)) {
  for (j in 2:ncol(redeploys)) {
    if (redeploys[i, j] > 1) {
      ndeploys <- sum(events_park$StdLocName == redeploys$StdLocName[i])
      events_park[events_park$StdLocName == redeploys$StdLocName[i],
                    names(redeploys)[j]] <- c(rep(0, ndeploys - 1), 1)
    }  
  }
}

# Create a matrix with deployment experience values [most experienced person in 
# group] (0 = novice; 1 = experienced; 2 = expert)
deploy_exp <- matrix(NA, nrow = nrow(dh), ncol = ncol(dh))
rownames(deploy_exp) <- rownames(dh) 
colnames(deploy_exp) <- paste0(YEAR, "_", 1:ncol(occ_matrix))
for (i in 1:nrow(deploy_exp)) {
  for (j in 1:ncol(deploy_exp)) {
    if (is.na(dh[i, j])) next
    deploy_exp[i, j] <- 
      events_park$deploy_exp[events_park$StdLocName == rownames(deploy_exp)[i] &
                               events_park[,colnames(deploy_exp)[j]] == 1]
  }
}

#------------------------------------------------------------------------------#
# Create day-of-year variable (using midpoint of each occasion)
#------------------------------------------------------------------------------#

occasions <- occasions %>%
  mutate(start_yday = yday(start),
         end_yday = yday(end),
         mid_yday = round((start_yday + end_yday)/2))

day <- matrix(occasions$mid_yday, 
              nrow = nrow(dh), 
              ncol = ncol(dh), 
              byrow = TRUE)

# Create standardized value of day variable
day_mn <- mean(day)
day_sd <- sd(day)
day_z <- (day - day_mn)/day_sd

#------------------------------------------------------------------------------#
# Spatial covariates (time invariant)
#------------------------------------------------------------------------------#

# Load multi-layer raster with spatial data
spat_raster <- readRDS(paste0("data/covariates/spatial-cov-", PARK, ".rds"))

# Create dataframe with covariate values for each camera location
spatial_covs <- locs_park
# Ensure the order is the same as what's in the detection history matrix
spatial_covs <- spatial_covs[match(rownames(dh), spatial_covs$loc),]

# Extract covariate values for each camera location
spatial_covs <- cbind(spatial_covs, 
                      terra::extract(x = spat_raster, 
                                     y = spatial_covs[,c("long", "lat")],
                                     ID = FALSE))

# Identify continuous covariates that we want to standardize
covs_cont <- names(spatial_covs)
covs_cont <- str_subset(covs_cont, "loc|long|lat|vegclass", negate = TRUE)

# Scale continuous covariates by mean, SD
for (i in covs_cont) {
  meani <- mean(spatial_covs[,i])
  sdi <- sd(spatial_covs[,i])
  spatial_covs[,paste0(i, "_z")] <- (spatial_covs[,i] - meani) / sdi
}

# Create table with pairwise correlations between continuous covariates
correl <- round(cor(spatial_covs[,covs_cont]), 2)
cor_df <- as.data.frame(as.table(correl), stringsAsFactors = FALSE)
cor_df <- cor_df %>%
  filter(Freq != 1) %>%
  mutate(variable1 = pmin(Var1, Var2), .before = Var1) %>%
  mutate(variable2 = pmax(Var1, Var2), .before = Var1) %>%
  select(-c(Var1, Var2)) %>%
  distinct() %>%
  arrange(variable1, variable2) %>%
  rename(corr = Freq)
# View those pairs with high correlations
cor_df %>%
  arrange(desc(corr)) %>%
  filter(abs(corr) > 0.5)

#------------------------------------------------------------------------------#
# Create data object for spOccupancy package
#------------------------------------------------------------------------------#

# First put covariates that could be used in the detection model in a list
# Elements can be n_sites * n_occasions matrices (for survey covariates)
# or vectors of length n_sites (for spatial covariates)
det_covs <- list(day_z = day_z,
                 deploy_exp = deploy_exp,
                 effort_z = effort_z,
                 boundary_z = spatial_covs$boundary_z,
                 east_z = spatial_covs$east_z,
                 elev_z = spatial_covs$elev_z,
                 north_z = spatial_covs$north_z,
                 pois_z = spatial_covs$pois_z,
                 roads_z = spatial_covs$roads_z,
                 slope_z = spatial_covs$slope_z,
                 trail_z = spatial_covs$trail_z)
if (PARK == "CHIR") {
  det_covs <- c(det_covs, 
                list(burn_severity = spatial_covs$burn_severity_2011))
}
if (PARK == "SAGW") {
  det_covs <- c(det_covs, 
                list(vegclass2 = spatial_covs$vegclass2),
                list(vegclass3 = spatial_covs$vegclass3),
                list(wash_z = spatial_covs$wash_z))
}

# Create data object (also a list)
data_list <- list(y = dh,
                  occ.covs = spatial_covs,
                  det.covs = det_covs)

#------------------------------------------------------------------------------#
# Create alternative models/formulas for occupancy and detection
#------------------------------------------------------------------------------#

# Names of spatial covariates
boundary <- "boundary_z"
aspect <- "east_z + north_z"
elevation <- "elev_z"
elevation2 <- "elev_z + I(elev_z ^ 2)"
pois <- "pois_z"
roads <- "roads_z"
slope <- "slope_z"
slope2 <- "slope_z + I(slope_z ^ 2)"
trail <- "trail_z"
if (PARK == "CHIR") {
  burn <- "burn_severity"
}
if (PARK == "SAGW") {
  wash <- "wash_z"
  veg <- "vegclass2 + vegclass3"
}

# Names of all survey covariates (as we'll likely want to include them all)
day <- "day_z + I(day_z^2)"
deploy <- "deploy_exp"
effort <- "effort_z"
survey_full <- paste(c(day, deploy, effort), collapse = "+")

# Create a list of occupancy models (Here, trying each spatial covariate 
# individually, and a few combined with vegclasses)
occ_formula_list <- list(as.formula(paste0("~ ", boundary)),
                         as.formula(paste0("~ ", aspect)),
                         as.formula(paste0("~ ", elevation)),
                         as.formula(paste0("~ ", pois)),
                         as.formula(paste0("~ ", roads)),
                         as.formula(paste0("~ ", slope)),
                         as.formula(paste0("~ ", trail)),
                         as.formula(paste0("~ ", wash)),
                         as.formula(paste0("~ ", veg)),
                         as.formula(paste0("~ ", aspect)),
                         as.formula(paste0("~ ", paste(c(veg, pois), collapse = "+"))),
                         as.formula(paste0("~ ", paste(c(veg, boundary), collapse = "+"))),
                         as.formula(paste0("~ ", paste(c(veg, roads), collapse = "+"))),
                         as.formula(paste0("~ ", paste(c(veg, trail), collapse = "+"))))

# Create a detection model (just trying one for now)
det_formula <- as.formula(paste0("~ ", survey_full))

#------------------------------------------------------------------------------#
# Run models in spOccupancy
#------------------------------------------------------------------------------#

# Set MCMC parameters
n_samples <- 5000
n_burn <- 3000
n_thin <- 1
n_chains <- 3

# Run one model 
# Not specifying priors, but using defaults which are N(0, var = 2.72)
# Not specifying initial values -- by default they come from priors
# Running chains sequentially (n.omp.threads = 1) because vignette states
# this only speeds things up in spatial models
out <- PGOcc(occ.formula = occ_formula_list[[1]],
             det.formula = det_formula, 
             data = data_list, 
             # inits = inits, 
             n.samples = n_samples, 
             # priors = priors, 
             n.omp.threads = 1, 
             verbose = TRUE, 
             n.report = 1000, 
             n.burn = n_burn, 
             n.thin = n_thin, 
             n.chains = n_chains) 
# Took <10 seconds to run
summary(out)

# Running a list of models
out_list <- list()
for (i in 1:length(occ_formula_list)) {
  out <- PGOcc(occ.formula = occ_formula_list[[i]],
               det.formula = det_formula, 
               data = data_list, 
               n.samples = n_samples, 
               n.omp.threads = 1, 
               verbose = TRUE, 
               n.report = 1000, 
               n.burn = n_burn, 
               n.thin = n_thin, 
               n.chains = n_chains) 
  out_list <- c(out_list, list(out))
}
summary(out_list[[1]])


# NEED TO FIX IN COVARIATE OBJECTS ##################################
# There are missing values in data$y with corresponding non-missing values in data$det.covs.
# Removing these site/replicate combinations for fitting the model.
  