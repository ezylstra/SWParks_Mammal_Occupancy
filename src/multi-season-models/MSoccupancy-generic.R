################################################################################
# Multi-season occupancy analysis

# ER Zylstra
# Updated 2022-10-07
################################################################################

#------------------------------------------------------------------------------#
# Create detection histories for selected park, species
#------------------------------------------------------------------------------#

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

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

# Type of camera (new cameras deployed in 2022, all the same before that)
# TODO: check that this is correct and there were no other equipment changes
surveys$camera_new <- ifelse(surveys$yr < 2022, 0, 1)

#------------------------------------------------------------------------------#
# Seasonal (annual) covariates
#------------------------------------------------------------------------------#

# We want to think about variables that will explain transitions between 
# occupancy status from one season to the next (extinction/colonization). Note 
# that the number of transitions is equal to the number of seasons - 1.

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

# Compile monsoon precipitation data (could do the same for other types of 
# weather data if available)
weather_var <- "monsoon_ppt"
monsoon_files <- weather_files[str_detect(weather_files, weather_var)]
# Remove monsoon rasters associated with periods outside the years of interest
# (monsoon rainfall in year x could explain transitions between years x and x+1)
monsoon_yrs <- paste0(as.character(transition_starts), collapse = "|")
monsoon_files <- monsoon_files[str_detect(monsoon_files, monsoon_yrs)]

# Load each raster and compile into a list
monsoon_list <- list()
for (i in 1:length(monsoon_files)) {
  monsoon_list[[i]] <- rast(monsoon_files[i])
  names(monsoon_list[[i]]) <- weather_var
}  

# Create a list that will store covariate values for each location in each 
# transition season
sitetrans <- list()
for (i in 1:length(monsoon_list)) {
  sitetrans[[i]] <- sites
  sitetrans[[i]]$start_yr <- transition_starts[i] 
  sitetrans[[i]]$end_yr <- transition_starts[i] + 1
  sitetrans[[i]]$trans_index <- i
  sitetrans[[i]] <- cbind(sitetrans[[i]],
                          terra::extract(x = monsoon_list[[i]],
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
    
# Remove monsoon_list from workspace, and remove rasters from local repo
rm(monsoon_list)
invisible(file.remove(list.files(weather_folder, full.names = TRUE)))
    
#------------------------------------------------------------------------------#
# Spatial covariates
#------------------------------------------------------------------------------#

# Load rasters with spatial data
  
  # Create needed lists of folder and file names 
  park_folder <- paste0("data/covariates/rasters-", park, "/")
  park_zip <- paste0("data/covariates/rasters-", park, ".zip")

  raster_filenames <- c("dist_boundary", "dist_pois", "dist_roads", "dist_trail", 
                        "east", "north", "slope")
  if (park == "SAGW") {
    raster_filenames <- c(raster_filenames, "vegclasses", "dist_wash")
  }
  if (park == "CHIR") {
    raster_filenames <- c(raster_filenames, "burn_severity_2011")
  }

  park_rasters <- c(paste0(park, "_DEM_1as.tif"), 
                    paste0(raster_filenames, "_", tolower(park), ".tif"))
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
if (park == "SAGW") {
  spatial_covs <- spatial_covs %>%
    mutate(vegclass2 = ifelse(vegclasses == 2, 1, 0),
           vegclass3 = ifelse(vegclasses == 3, 1, 0)) %>%
    select(-vegclasses)
}
  
# Remove raster_list from workspace, and remove rasters from local repo
rm(raster_list)
invisible(file.remove(list.files(park_folder, full.names = TRUE)))

# Attach spatial covariates to surveys dataframe so we can use them as 
# covariates for detection
if (park == "SAGW") {
  surveys <- left_join(surveys, 
                       spatial_covs[,c("loc", paste0(covs_cont, "_z"), "vegclass2", "vegclass3")])
} else {
  surveys <- left_join(surveys, 
                       spatial_covs[,c("loc", paste0(covs_cont, "_z"))])  
}


# Attach subset of spatial covariates to sitetrans dataframe so we can use them 
# as covariates for extinction/colonization
covs_extcol <- c("elev")
sitetrans <- left_join(sitetrans, 
                       spatial_covs[,c("loc", paste0(covs_extcol, "_z"))])

#------------------------------------------------------------------------------#
# Extract covariate values for each parameter
#------------------------------------------------------------------------------#

# Initial occupancy (psi)  
  # Create vector of standardized covariate names
  covariates_psi <- ifelse(covariates_psi %in% covs_cont,
                           paste0(covariates_psi, "_z"), covariates_psi)
  if (!is.na(psi_quadratics)) {
    psi_quadratics <- covariates_psi[psi_quadratics]
  } 
  # Extract covariate values
  cov_psi <- spatial_covs %>%
    select(contains(covariates_psi))
  # Add quadratics if needed
  if (!is.na(psi_quadratics)) {
    cov_psi[,paste0(psi_quadratics, "2")] <- cov_psi[,psi_quadratics] ^ 2
  }
  # Put columns in alphabetical order and convert to a matrix
  cov_psi <- cov_psi %>%
    select(order(colnames(.))) %>%
    as.matrix

# Detection probability (p)  
  # Create vector of standardized covariate names
  covariates_p <- ifelse(covariates_p %in% c(covs_cont, "effort", "day"),
                         paste0(covariates_p, "_z"), covariates_p)
  if (!is.na(p_quadratics)) {
    p_quadratics<- covariates_p[p_quadratics]
  } 
  # Extract covariate values
  cov_p <- surveys %>%
    select(contains(covariates_p))
  # Add quadratics if needed
  if (!is.na(p_quadratics)) {
    cov_p[,paste0(p_quadratics, "2")] <- cov_p[,p_quadratics] ^ 2
  }
  # Put columns in alphabetical order and convert to a matrix
  cov_p <- cov_p %>%
    select(order(colnames(.))) %>%
    as.matrix  

# Extinction probability (eps)  
  # Create vector of standardized covariate names
  covariates_eps <- ifelse(covariates_eps %in% c(covs_cont, "monsoon_ppt"),
                           paste0(covariates_eps, "_z"), covariates_eps)  
  if (!is.na(eps_quadratics)) {
    eps_quadratics<- covariates_eps[eps_quadratics]
  }   
  # Extract covariate values  
  cov_eps <- sitetrans %>%
    select(contains(covariates_eps))
  # Add quadratics if needed
  if (!is.na(eps_quadratics)) {
    cov_eps[,paste0(eps_quadratics, "2")] <- cov_eps[,eps_quadratics] ^ 2
  }
  # Put columns in alphabetical order
  cov_eps <- cov_eps %>%
    select(order(colnames(.)))
  # Add interactions if needed
  if (n_eps_interacts > 0) {
    eps_interacts <- matrix(NA, nrow = 1, ncol = 2)
    for (i in 1:n_eps_interacts) {
      eps_interacts[i,] <- get(paste0("eps_int", i))
      eps_interacts[i,1] <- ifelse(eps_interacts[i,1] %in% c(covs_cont, "monsoon_ppt"),
                                   paste0(eps_interacts[i,1], "_z"),
                                   eps_interacts[i,1])
      eps_interacts[i,2] <- ifelse(eps_interacts[i,2] %in% c(covs_cont, "monsoon_ppt"),
                                   paste0(eps_interacts[i,2], "_z"),
                                   eps_interacts[i,2])    
      cov_eps[,paste0(eps_interacts[i,], collapse = "_")] <- 
        cov_eps[,eps_interacts[i,1]] * cov_eps[,eps_interacts[i,2]]
    }
  }
  # Convert to a matrix
  cov_eps <- as.matrix(cov_eps)

# Colonization probability (gam)  
  # Create vector of standardized covariate names
  covariates_gam <- ifelse(covariates_gam %in% c(covs_cont, "monsoon_ppt"),
                           paste0(covariates_gam, "_z"), covariates_gam)  
  if (!is.na(gam_quadratics)) {
    gam_quadratics<- covariates_gam[gam_quadratics]
  }   
  # Extract covariate values  
  cov_gam <- sitetrans %>%
    select(contains(covariates_gam))
  # Add quadratics if needed
  if (!is.na(gam_quadratics)) {
    cov_gam[,paste0(gam_quadratics, "2")] <- cov_gam[,gam_quadratics] ^ 2
  }
  # Put columns in alphabetical order
  cov_gam <- cov_gam %>%
    select(order(colnames(.)))
  # Add interactions if needed
  if (n_gam_interacts > 0) {
    gam_interacts <- matrix(NA, nrow = 1, ncol = 2)
    for (i in 1:n_gam_interacts) {
      gam_interacts[i,] <- get(paste0("gam_int", i))
      gam_interacts[i,1] <- ifelse(gam_interacts[i,1] %in% c(covs_cont, "monsoon_ppt"),
                                   paste0(gam_interacts[i,1], "_z"),
                                   gam_interacts[i,1])
      gam_interacts[i,2] <- ifelse(gam_interacts[i,2] %in% c(covs_cont, "monsoon_ppt"),
                                   paste0(gam_interacts[i,2], "_z"),
                                   gam_interacts[i,2])    
      cov_gam[,paste0(gam_interacts[i,], collapse = "_")] <- 
        cov_gam[,gam_interacts[i,1]] * cov_gam[,gam_interacts[i,2]]
    }
  }
  # Convert to a matrix
  cov_gam <- as.matrix(cov_gam)

# Calculate the number of covariates for each parameter  
n_cov_psi <- ncol(cov_psi)
n_cov_p <- ncol(cov_p)
n_cov_eps <- ncol(cov_eps)
n_cov_gam <- ncol(cov_gam)

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

# print(out)

# Trace and density plots
# MCMCtrace(out,pdf = FALSE)
