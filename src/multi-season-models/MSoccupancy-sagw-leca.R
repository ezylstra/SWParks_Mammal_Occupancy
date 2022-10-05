################################################################################
# Multi-season occupancy analysis
# SAGW, LECA (black-tailed jackrabbit)

# ER Zylstra
# Updated 2022-07-12
################################################################################

# Note: Will probably want to create a template script for this type of analysis
# This will serve as a first case study

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(terra)
library(jagsUI)
library(MCMCvis)

# rm(list = ls())

# Load photo, location, events, species data 
source("src/photo-data/format-mammal-data.R")

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
years <- 2017:2022

#------------------------------------------------------------------------------#
# Create detection histories for selected park, species
#------------------------------------------------------------------------------#

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
# may also result in shorter run times (because we don't have to cycle through
# occasion-site combinations when the camera wasn't operational). 

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

# Remove rows with det = NA (no survey data for that occasion/location)
surveys <- filter(surveys, !is.na(det))

# Scale effort covariate by mean, SD
# (note that this results in standardized values that range from -10 to 0.25)
surveys <- surveys %>%
  mutate(effort_z = (effort - mean(effort))/sd(effort))

# Add year and trend columns
surveys$yr <- as.numeric(str_sub(surveys$occ, 1,4))
surveys$trend <- surveys$yr - min(surveys$yr)

# Add site, season, occasion indicators
surveys$season_index <- surveys$trend + 1
surveys$occ_index <- as.numeric(str_sub(surveys$occ,6,nchar(surveys$occ)))
for (i in 1:nrow(surveys)) {
  surveys$site_index[i] <- which(rownames(dh) == surveys$loc[i])
}

#------------------------------------------------------------------------------#
# Add occasion-specific covariates
#------------------------------------------------------------------------------#

# Day of the year (use midpoint of each occasion)
occasions <- occasions %>%
  mutate(start_yday = yday(start),
         end_yday = yday(end),
         mid_yday = round((start_yday + end_yday)/2),
         mid_yday_z = (mid_yday - mean(mid_yday))/sd(mid_yday))

surveys$day_z <- occasions$mid_yday_z[match(surveys$occ, occasions$yr_occ)]
surveys$day_z2 <- surveys$day_z * surveys$day_z

# Type of camera (new cameras deployed in 2022, all the same before that)
# TODO: check that this is correct and there were no other equipment changes
surveys$camera_new <- ifelse(surveys$yr < 2022, 0, 1)

#------------------------------------------------------------------------------#
# Seasonal (annual) covariates
#------------------------------------------------------------------------------#

# We want to think about variables that will explain transitions between 
# occupancy status from one season to the next (extinction/colonization), so the 
# number of transitions is equal to the number of seasons - 1.

# For now it's just weather covariates, but we can add other types of covariates
# at any time

# Create a dataframe that will contain seasonal covariates
sites <- surveys %>%
  select(loc, site_index) %>%
  distinct %>%
  left_join(., locs_park)

transition_starts <- min(occasions$yr):(max(occasions$yr)-1)

# Load rasters with seasonal (and spatial) weather data
  # Create lists of folder and file names 
  weather_folder <- "data/covariates/weather-derived-rasters/"
  weather_zip <- "data/covariates/weather-derived.zip"
  
  # Unzip weather folder first, if necessary
  if (length(list.files(weather_folder)) == 0) {
    unzip(weather_zip, overwrite = TRUE)
  }
  
# List files in weather folder
weather_files <- list.files(weather_folder, full.names = TRUE)
weather_files <- weather_files[str_detect(weather_files, "2016", negate = TRUE)]
  
# Compile monsoon precipitation data (could do the same for other types of 
# weather data if available)
weather_var <- "monsoon_ppt"
weather_subset <- weather_files[str_detect(weather_files, weather_var)]

# Load each raster and compile into a list
weather_subset_list <- list()
for (i in 1:length(weather_subset)) {
  weather_subset_list[[i]] <- rast(weather_subset[i])
  names(weather_subset_list[[i]]) <- weather_var
}  

# Create a list that will store covariate values for each location in each 
# transition season
sitetrans <- list()
for (i in 1:length(weather_subset_list)) {
  sitetrans[[i]] <- sites
  sitetrans[[i]]$start_yr <- transition_starts[i] 
  sitetrans[[i]]$end_yr <- transition_starts[i] + 1
  sitetrans[[i]]$trans_index <- i
  sitetrans[[i]] <- cbind(sitetrans[[i]],
                          terra::extract(x = weather_subset_list[[i]],
                                         y = sitetrans[[i]][,c("long", "lat")],
                                         ID = FALSE))
}
    
# Combine annual dataframes in list to a single dataframe
sitetrans <- bind_rows(sitetrans)

# Scale continuous covariates by mean, SD
for (i in weather_var) {
  meani <- mean(sitetrans[,i])
  sdi <- sd(sitetrans[,i])
  sitetrans[,paste0(i, "_z")] <- (sitetrans[,i] - meani) / sdi
}
    
# Remove weather_subset_list from workspace, and remove rasters from local repo
rm(weather_subset_list)
invisible(file.remove(list.files(weather_folder, full.names = TRUE)))
    
#------------------------------------------------------------------------------#
# Spatial covariates
#------------------------------------------------------------------------------#

# Load rasters with spatial data
  
  # Create needed lists of folder and file names 
  park_folder <- "data/covariates/rasters-SAGW/"
  park_zip <- "data/covariates/rasters-SAGW.zip"
  raster_filenames <- c("dist_boundary", "dist_pois", "dist_roads", "dist_trail", 
                        "east", "north", "slope", "vegclasses", "dist_wash")
  park_rasters <- c("SAGW_DEM_1as.tif", 
                    paste0(raster_filenames, "_sagw.tif"))
  park_rasters <- paste0(park_folder, park_rasters)

  # Unzip SAGW folder first, if necessary
  if (!all(unlist(lapply(X = park_rasters, FUN = file.exists)))) {
    unzip(park_zip, overwrite = TRUE)
  }

  raster_objects <- ifelse(str_detect(raster_filenames, "dist"), 
                           str_remove(raster_filenames, "dist_"),
                           raster_filenames)
  raster_objects <- c("elev", raster_objects)

  # Load each raster and compile into a list
  raster_list <- list()
  for (i in 1:length(raster_objects)) {
    raster_list[[i]] <- rast(park_rasters[i])
    names(raster_list[[i]]) <- raster_objects[i]
  }

spatial_covs <- locs_park

# Ensure the order is the same as what's in the detection history matrix
spatial_covs <- spatial_covs[match(rownames(dh), spatial_covs$loc),]

# Extract covariate values for each camera location
for (i in 1:length(raster_list)) {
  spatial_covs <- cbind(spatial_covs, 
                        terra::extract(x = raster_list[[i]], 
                                       y = spatial_covs[,c("long", "lat")],
                                       ID = FALSE))
}

# Identify continuous covariates
covs_cont <- c("long", "lat", raster_objects[raster_objects != "vegclasses"])

# Check out pairwise correlations
correl <- round(cor(spatial_covs[,covs_cont]),2)
cor_df <- as.data.frame(as.table(correl))
# Look more carefully at |r| > 0.5
cor_df %>% 
  arrange(desc(abs(Freq))) %>% 
  filter(Freq != 1 & abs(Freq) > 0.5)
  # should probably only use one of: elev, slope (and maybe roads?)
  # should probably only use one of: pois, trail

# Scale continuous covariates by mean, SD
for (i in covs_cont) {
  meani <- mean(spatial_covs[,i])
  sdi <- sd(spatial_covs[,i])
  spatial_covs[,paste0(i, "_z")] <- (spatial_covs[,i] - meani) / sdi
}

# For vegetation classes, make 1 = low gradient desert the reference level
# and create indicators for 2 (low hillslope, foothills) and 3 (med-high 
# gradient, hilly)
spatial_covs <- spatial_covs %>%
  mutate(vegclass2 = ifelse(vegclasses == 2, 1, 0),
         vegclass3 = ifelse(vegclasses == 3, 1, 0)) %>%
  select(-vegclasses)

# Remove raster_list from workspace, and remove rasters from local repo
rm(raster_list)
invisible(file.remove(list.files(park_folder, full.names = TRUE)))

# Attach spatial covariates to surveys dataframe so we can use them as 
# covariates for detection
surveys <- left_join(surveys, 
                     spatial_covs[,c("loc", paste0(covs_cont, "_z"), "vegclass2", "vegclass3")])

# Attach subset of spatial covariates to sitetrans dataframe so we can use them 
# as covariates for extinction/colonization
covs_extcol <- c("elev")
sitetrans <- left_join(sitetrans, 
                       spatial_covs[,c("loc", paste0(covs_extcol, "_z"))])

#------------------------------------------------------------------------------#
# Package things up for JAGS
#------------------------------------------------------------------------------#

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
# known_state_occ(y_array)

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
# inits_state_occ(y_array)

# Identify covariates for each model parameter
  # Initial occupancy (psi)
  covariates_psi <- c("elev", "pois", "vegclass")
  covariates_psi <- ifelse(covariates_psi %in% covs_cont,
                           paste0(covariates_psi, "_z"), covariates_psi)
  # Detection probability (p)
  covariates_p <- c("effort", "day", "camera_new")
  covariates_p <- ifelse(covariates_p %in% c(covs_cont, "effort", "day"),
                         paste0(covariates_p, "_z"), covariates_p)
  # Extinction probability (eps)
  covariates_eps <- c("monsoon_ppt", "elev")
  covariates_eps <- ifelse(covariates_eps %in% c(covs_cont, "monsoon_ppt"),
                           paste0(covariates_eps, "_z"), covariates_eps)
  # Colonization probability (gam)
  covariates_gam <- c("elev")
  covariates_gam <- ifelse(covariates_gam %in% c(covs_cont, "monsoon_ppt"),
                           paste0(covariates_gam, "_z"), covariates_gam)

# Extract matrices of covariate values
# Can create quadratics or interactions here, if necessary
cov_psi <- spatial_covs %>%
  select(contains(covariates_psi)) %>%
  as.matrix
n_cov_psi <- ncol(cov_psi)

cov_p <- surveys %>%
  select(contains(covariates_p)) %>%
  as.matrix
n_cov_p <- ncol(cov_p)

cov_eps <- sitetrans %>%
  select(contains(covariates_eps)) %>%
  as.matrix
n_cov_eps <- ncol(cov_eps)

cov_gam <- sitetrans %>%
  select(contains(covariates_gam)) %>%
  as.matrix
n_cov_gam <- ncol(cov_gam)

# Bundle data for JAGS
jags_data <- list(y = surveys$det,
                  n_sites = n_sites,
                  n_seasons = n_seasons,
                  n_obs = nrow(surveys),
                  n_sitetrans = nrow(sitetrans),
                  cov_psi = cov_psi,
                  n_cov_psi = ncol(cov_psi),
                  cov_p = cov_p,
                  n_cov_p = ncol(cov_p),
                  cov_eps = cov_eps,
                  n_cov_eps = ncol(cov_eps),
                  cov_gam = cov_gam,
                  n_cov_gam = ncol(cov_gam),
                  site = surveys$site_index,
                  season = surveys$season_index,
                  site_ec = sitetrans$site_index,
                  trans_ec = sitetrans$trans_index,
                  z = known_state_occ(y_array))

# List of parameters to monitor
params <- c("mean_psi", "beta_psi0", "beta_psi",
            "mean_p", "beta_p0", "beta_p",
            "mean_eps", "beta_eps0", "beta_eps",
            "mean_gam", "beta_gam0", "beta_gam",
            "PAO")

# Initial values
inits <- function(){list(mean_psi = runif(1, 0, 1),
                         beta_psi = runif(n_cov_psi, -2, -2),
                         mean_p = runif(1, 0, 1),
                         beta_p = runif(n_cov_p, -2, 2),
                         mean_eps = runif(1, 0, 1),
                         beta_eps = runif(n_cov_eps, -2, 2),
                         mean_gam = runif(1, 0, 1),
                         beta_gam = runif(n_cov_gam, -2, 2),
                         z = inits_state_occ(y_array))}

#------------------------------------------------------------------------------#
# Run model in JAGS
#------------------------------------------------------------------------------#

nc <- 3      # Number of chains
na <- 3000   # Number of iterations to run in the adaptive phase
nb <- 10000  # Number of iterations to discard (burn-in)
ni <- 20000  # Number of iterations per chain (including burn-in)
nt <- 10     # Thinning rate

out <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "JAGS/JAGS_MultiSeasonWithCovs.txt",
            n.chains = nc,
            n.adapt = na,
            n.burnin = nb,
            n.iter = ni,
            n.thin = nt,
            parallel = TRUE)

print(out)

# Trace and density plots
# MCMCtrace(out,pdf = FALSE)

# Save model to file in output/models
model_file <- paste0("output/models/",
                     tolower(park), "-",
                     tolower(species), "-",
                     "MS-test2.rds")
saveRDS(object = out,
        file = model_file)

#------------------------------------------------------------------------------#
# Interpreting/plotting effects of covariates (examples)
#------------------------------------------------------------------------------#

# Load JAGS model if we ran it previously:
# model_file <- paste0("output/models/",
#                      tolower(park), "-",
#                      tolower(species), "-",
#                      "MS-test2.rds")
# out <- readRDS(file = model_file)

# Extract posterior samples from jagsUI object and combine into one dataframe
samples <- out$samples
samples <- do.call(rbind, samples)

# Estimate detection probability for a particular level of effort
# (assuming average day number and cameras used prior to 2022)

  # Camera operating for all 7 days
  eff7 <- 1
  eff7_z <- (eff7 - mean(surveys$effort)) / sd(surveys$effort)
    # logit(p[w]) <- beta_p0 + beta_p[1] * lat[w] + beta_p[2] * effort[w]
    # Assuming latitude = mean (so lat_z = 0), we're left with:
    # logit(p[eff7_z]) <- beta_p0 + beta_p[2] * eff7_z
  # Create a matrix of covariate values (including the intercept)
  X_p <- cbind(int = 1, effort = eff7_z)
  # Create a matrix with posterior samples for the parameters we need
  betas_p <- samples[,c("beta_p0", "beta_p[2]")]
  # A little matrix math gives you a vector of predicted values on logit scale
  # Vector length is equal to the number of posterior samples (here, 3000)
  pred_logit_p <- X_p %*% t(betas_p)
  # Backtransform to the probability scale
  pred_prob_p <- exp(pred_logit_p) / (1 + exp(pred_logit_p))
  # Summarize these values to get an estimate with 95% credible interval
  mean(pred_prob_p)
  quantile(pred_prob_p, probs = c(0.025, 0.975))
  
  # Camera operating for 3 of 7 days  
  eff3 <- 3/7
  eff3_z <- (eff3 - mean(surveys$effort)) / sd(surveys$effort)
  X_p <- cbind(int = 1, effort = eff3_z)
  pred_logit_p <- X_p %*% t(betas_p)
  pred_prob_p <- exp(pred_logit_p) / (1 + exp(pred_logit_p))
  mean(pred_prob_p)
  quantile(pred_prob_p, probs = c(0.025, 0.975))

# Estimate how colonization probability would change with a 1-unit (SD) 
# increase in elev
  
  # Here it's important to remember that logit is the log of the odds
  # Odds are the probability event happens / probability event doesn't happen
  # So logit(gamma[i]) = log(gamma[i] / (1-gamma[i])) = log(odds)
  # Odds a site gets colonized at mean elevation = exp(beta_gam0)
  # Odds a site gets colonized with 1-SD increase in elevation = 
    # exp(beta_gam0 + beta_gam1 * 1) = exp(beta_gam0)*exp(beta_gam1)
  # So the odds will change by a FACTOR of exp(beta_gam1):
  # Odds[elev+1SD] = Odds[mean elev] * exp(beta_gam1)
  
  beta_gam <- samples[,c("beta_gam")]
  change <- exp(beta_gam)
  mean(change) # 0.48
  quantile(change, probs = c(0.025, 0.975)) #0.27, 0.76
  # The odds a site is colonized are estimated to be 52% lower for each 
  # 1-SD increase in elevation (95% CI = 24-73%)
  
# Plot effect of distance-to-POI on initial occupancy (assuming vegclass = 1
# and mean elevation)

  # logit(psi[i]) <- beta_psi0 + beta_psi[2] * distance[i]
  
  # Generate a vector of distances that span the range at surveyed locations
  dists <- seq(min(spatial_covs$pois), max(spatial_covs$pois), length = 100)
  # Standardize these values
  dists_z <- (dists - mean(spatial_covs$pois)) / sd(spatial_covs$pois)
  # Create a matrix of covariate values (including the intercept [1])
  X_psi <- cbind(int = 1, dist = dists_z)
  # Create a matrix with posterior samples for the parameters we need
  betas_psi <- samples[,c("beta_psi0", "beta_psi[2]")]
  # A little matrix math gives you a matrix of predicted values on logit scale
  # Matrix dimensions = 100 x 3000 (3000 predicted values for each lat in seq)
  pred_logit_psi <- X_psi %*% t(betas_psi)
  # Backtransform to the probability scale
  pred_prob_psi <- exp(pred_logit_psi) / (1 + exp(pred_logit_psi))
  # Calculate mean and credible interval for each value
  mean_psi <- apply(pred_prob_psi, 1, mean)
  cri_psi <- apply(pred_prob_psi, 1, quantile, probs = c(0.025, 0.975)) 
  
  # Plot predictions
  par(mfrow=c(1,1))
  plot(mean_psi ~ dists, type = "l", bty = "l", ylim = c(0, 1), 
       xlab = "Distance to POI", ylab = "Predicted occupancy", las = 1)
  polygon(c(dists, rev(dists)), c(cri_psi[1,], rev(cri_psi[2,])), 
          col = rgb(0, 0, 0, 0.2), border = NA)
  
#------------------------------------------------------------------------------#
# For comparison, running a simple model in JAGS and unmarked
#------------------------------------------------------------------------------#
  
  # Identify covariates for each model parameter
  # Initial occupancy (psi)
  covariates_psi <- c("elev")
  covariates_psi <- ifelse(covariates_psi %in% covs_cont,
                           paste0(covariates_psi, "_z"), covariates_psi)
  # Detection probability (p)
  covariates_p <- c("effort")
  covariates_p <- ifelse(covariates_p %in% c(covs_cont, "effort", "day"),
                         paste0(covariates_p, "_z"), covariates_p)
  # Extinction probability (eps)
  covariates_eps <- c("elev")
  covariates_eps <- ifelse(covariates_eps %in% c(covs_cont, "monsoon_ppt"),
                           paste0(covariates_eps, "_z"), covariates_eps)
  # Colonization probability (gam)
  covariates_gam <- c("elev")
  covariates_gam <- ifelse(covariates_gam %in% c(covs_cont, "monsoon_ppt"),
                           paste0(covariates_gam, "_z"), covariates_gam)
  
  # Extract matrices of covariate values
  # Can create quadratics or interactions here, if necessary
  cov_psi <- spatial_covs %>%
    select(contains(covariates_psi)) %>%
    as.matrix
  n_cov_psi <- ncol(cov_psi)
  
  cov_p <- surveys %>%
    select(contains(covariates_p)) %>%
    as.matrix
  n_cov_p <- ncol(cov_p)
  
  cov_eps <- sitetrans %>%
    select(contains(covariates_eps)) %>%
    as.matrix
  n_cov_eps <- ncol(cov_eps)
  
  cov_gam <- sitetrans %>%
    select(contains(covariates_gam)) %>%
    as.matrix
  n_cov_gam <- ncol(cov_gam)
  
  # Bundle data for JAGS
  jags_data <- list(y = surveys$det,
                    n_sites = n_sites,
                    n_seasons = n_seasons,
                    n_obs = nrow(surveys),
                    n_sitetrans = nrow(sitetrans),
                    cov_psi = cov_psi,
                    n_cov_psi = ncol(cov_psi),
                    cov_p = cov_p,
                    n_cov_p = ncol(cov_p),
                    cov_eps = cov_eps,
                    n_cov_eps = ncol(cov_eps),
                    cov_gam = cov_gam,
                    n_cov_gam = ncol(cov_gam),
                    site = surveys$site_index,
                    season = surveys$season_index,
                    site_ec = sitetrans$site_index,
                    trans_ec = sitetrans$trans_index,
                    z = known_state_occ(y_array))
  
  # List of parameters to monitor
  params <- c("mean_psi", "beta_psi0", "beta_psi",
              "mean_p", "beta_p0", "beta_p",
              "mean_eps", "beta_eps0", "beta_eps",
              "mean_gam", "beta_gam0", "beta_gam",
              "PAO")
  
  # Initial values
  inits <- function(){list(mean_psi = runif(1, 0, 1),
                           beta_psi = runif(n_cov_psi, -2, -2),
                           mean_p = runif(1, 0, 1),
                           beta_p = runif(n_cov_p, -2, 2),
                           mean_eps = runif(1, 0, 1),
                           beta_eps = runif(n_cov_eps, -2, 2),
                           mean_gam = runif(1, 0, 1),
                           beta_gam = runif(n_cov_gam, -2, 2),
                           z = inits_state_occ(y_array))}
  
  nc <- 3      # Number of chains
  na <- 3000   # Number of iterations to run in the adaptive phase
  nb <- 10000  # Number of iterations to discard (burn-in)
  ni <- 20000  # Number of iterations per chain (including burn-in)
  nt <- 10     # Thinning rate
  
  out.simple <- jags(data = jags_data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "JAGS/JAGS_MultiSeasonWithCovs.txt",
                     n.chains = nc,
                     n.adapt = na,
                     n.burnin = nb,
                     n.iter = ni,
                     n.thin = nt,
                     parallel = TRUE)
  
  print(out.simple)
  
  library(unmarked)
  dh_full <- cbind(dh[,1:5], rep(NA, n_sites),
                   dh[,6:9], rep(NA, n_sites), rep(NA, n_sites),
                   matrix(NA, ncol = 6, nrow = n_sites),
                   dh[,10:15],
                   dh[,16:20], rep(NA, n_sites),
                   dh[,21:25], rep(NA, n_sites))
  
  eff_mean <- mean(c(effort[effort!=0]))
  eff_sd <- sd(c(effort[effort!=0]))
  effort_z <- (effort - eff_mean)/eff_sd
  eff_full <- cbind(effort_z[,1:5], rep(NA, n_sites),
                    effort_z[,6:9], rep(NA, n_sites), rep(NA, n_sites),
                    matrix(NA, ncol = 6, nrow = n_sites),
                    effort_z[,10:15],
                    effort_z[,16:20], rep(NA, n_sites),
                    effort_z[,21:25], rep(NA, n_sites))
  
  umf <- unmarkedMultFrame(y = dh_full,
                           siteCovs = data.frame(elev = spatial_covs$elev_z),
                           obsCovs = list(effort = eff_full),
                           numPrimary = max(surveys$season_index))
  summary(umf)
  summary(m <- colext(~ elev, ~ elev, ~ elev, ~ effort, data = umf))
  # Estimates are pretty similar
  