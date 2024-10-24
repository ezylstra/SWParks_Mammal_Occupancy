################################################################################
# Multi-season occupancy analysis for species, park, years, and covariates 
# specified in MSoccupancy-wrapper.R

# ER Zylstra
# Updated 2022-11-18
################################################################################

#------------------------------------------------------------------------------#
# Create detection histories for selected park, species, years
#------------------------------------------------------------------------------#

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

# Extract sampling occasion info for selected park and years
occasions <- occasions %>%
  filter(Park == PARK) %>%
  filter(yr %in% YEARS) %>%
  arrange(yr, occasion)

# Create a list of days included in sampling occasions
occ_days <- NULL
for (i in 1:nrow(occasions)) {
  occ_days <- append(occ_days, occasions$start_day[i]:occasions$end_day[i])
}

# Extract photo observations for park, species, years
# Retain a maximum of one observation per day at each location
obs <- dat %>% 
  filter(Park == PARK & Species_code == SPECIES & yr %in% YEARS) %>%
  select(StdLocName, obsdate, yr, o_day) %>%
  arrange(StdLocName, obsdate) %>%
  distinct

# Extract information about camera locations in selected park
events_park <- events %>%
  filter(Park == PARK) %>%
  filter(d_yr %in% YEARS)
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
                       ncol = ncol(ddh) / occasions$duration[1],
                       dimnames = list(rownames(ddh), NULL))

for (i in 1:ncol(dh)) {
  multiday <- ddh[,colnames(ddh) %in% occasions$start_day[i]:occasions$end_day[i]]
  dh[,i] <- apply(multiday, 1, paNA)
  effort[,i] <- apply(multiday, 1, propNA)
}

#------------------------------------------------------------------------------#
# Put detection and effort data into long form
#------------------------------------------------------------------------------#

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
surveys <- left_join(dh_long, eff_long, by = c("loc", "occ"))

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
# Incorporate information about deployment personnel
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
  
# For each row in the survey dataframe, attach the deployment personnel 
# experience value that corresponds with that camera location and occasion
for (i in 1:nrow(surveys)) {
  de <- events_park$deploy_exp[events_park$StdLocName == surveys$loc[i] &
                               events_park[,surveys$occ[i]] == 1]
  # Infrequently, there may be two events associated with a particular camera
  # location and occasion when a camera was redeployed immediately.  In these
  # instances, better to use the 2nd entry for deploy_exp (which reflects who
  # REdeployed the camera)
  de <- de[length(de)]
  surveys$deploy_exp[i] <- de
}

#------------------------------------------------------------------------------#
# Add occasion-specific covariates to detection data
#------------------------------------------------------------------------------#

# Day of the year (use midpoint of each occasion)
occasions <- occasions %>%
  mutate(start_yday = yday(start),
         end_yday = yday(end),
         mid_yday = round((start_yday + end_yday)/2),
         mid_yday_z = (mid_yday - mean(mid_yday))/sd(mid_yday))

surveys$day <- occasions$mid_yday[match(surveys$occ, occasions$yr_occ)]
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

# For now we just have weather covariates, but we can add other types of 
# covariates at any time

# Create a dataframe that will contain seasonal covariates
sites <- surveys %>%
  select(loc, site_index) %>%
  distinct %>%
  left_join(., locs_park, by = "loc")

transition_starts <- YEARS[-length(YEARS)]

# Load rasters with seasonal weather data (that also varies over space)
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
    
# Combine annual dataframes in sitetrans list to a single dataframe
sitetrans <- bind_rows(sitetrans)

# Scale continuous covariates by mean, SD
for (i in weather_var) {
  meani <- mean(sitetrans[,i])
  sdi <- sd(sitetrans[,i])
  sitetrans[,paste0(i, "_z")] <- (sitetrans[,i] - meani) / sdi
}
    
# Remove raster list from workspace, and remove rasters from local repo
rm(monsoon_list)
invisible(file.remove(list.files(weather_folder, full.names = TRUE)))
    
#------------------------------------------------------------------------------#
# Spatial covariates (time invariant)
#------------------------------------------------------------------------------#

# Load rasters with spatial data
  
  # Create needed lists of folder and file names 
  park_folder <- paste0("data/covariates/rasters-", PARK, "/")
  if (PARK == "ORPI") {
    park_zip1 <- "data/covariates/rasters-ORPI-dist.zip"
    park_zip2 <- "data/covariates/rasters-ORPI-topo.zip"
  } else {  
    park_zip <- paste0("data/covariates/rasters-", PARK, ".zip")
  }

  raster_filenames <- c("dist_boundary", "dist_pois", "dist_roads", "dist_trail", 
                        "east", "north", "slope")
  if (PARK == "SAGW") {
    raster_filenames <- c(raster_filenames, "vegclasses", "dist_wash")
  }
  if (PARK == "CHIR") {
    raster_filenames <- c(raster_filenames, "burn_severity_2011")
  }

  park_rasters <- c(paste0(PARK, "_DEM_1as.tif"), 
                    paste0(raster_filenames, "_", tolower(PARK), ".tif"))
  park_rasters <- paste0(park_folder, park_rasters)

  # Unzip park folder(s) first, if necessary
  if (!all(unlist(lapply(X = park_rasters, FUN = file.exists)))) {
    if (PARK == "ORPI") {
      unzip(park_zip1, overwrite = TRUE)
      unzip(park_zip2, overwrite = TRUE)
    } else {
      unzip(park_zip, overwrite = TRUE)
    }
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

# Scale continuous covariates by mean, SD
for (i in covs_cont) {
  meani <- mean(spatial_covs[,i])
  sdi <- sd(spatial_covs[,i])
  spatial_covs[,paste0(i, "_z")] <- (spatial_covs[,i] - meani) / sdi
}

# For vegetation classes, make 1 = low gradient desert the reference level
# and create indicators for 2 (low hillslope, foothills) and 3 (med-high 
# gradient, hilly)
if (PARK == "SAGW") {
  spatial_covs <- spatial_covs %>%
    mutate(vegclass2 = ifelse(vegclasses == 2, 1, 0),
           vegclass3 = ifelse(vegclasses == 3, 1, 0)) %>%
    select(-vegclasses)
}
  
# Remove raster list from workspace, and remove rasters from local repo
rm(raster_list)
invisible(file.remove(list.files(park_folder, full.names = TRUE)))

# Attach spatial covariates to surveys dataframe so we can use them as 
# covariates for detection
if (PARK == "SAGW") {
  surveys <- left_join(surveys, 
                       spatial_covs[,c("loc", paste0(covs_cont, "_z"), "vegclass2", "vegclass3")],
                       by = "loc")
} else {
  surveys <- left_join(surveys, 
                       spatial_covs[,c("loc", paste0(covs_cont, "_z"))],
                       by = "loc")  
}

# Attach subset of spatial covariates to sitetrans dataframe so we can use them 
# as covariates for extinction/colonization (for now, assuming the only spatial
# covariate that might be used for ext/col is elevation)
covs_extcol <- c("elev")
sitetrans <- left_join(sitetrans, 
                       spatial_covs[,c("loc", paste0(covs_extcol, "_z"))], 
                       by = "loc")

#------------------------------------------------------------------------------#
# Extract covariate values for each parameter 
# (covariates selected in MSoccupancy-wrapper.R)
#------------------------------------------------------------------------------#
covs_cont <- c(covs_cont, "effort", "day", "monsoon_ppt", "burn_severity_2011")
# Not including camera_new or deploy_exp in this list because we don't want to 
# standardize (0 indicates camera type used prior to 2022 and inexperienced 
# people deploying cameras)

# Initial occupancy (psi)
  if (!all(is.na(COVARS_PSI))) {
    # Create vector of standardized covariate names
    covars_psi <- ifelse(COVARS_PSI %in% covs_cont,
                         paste0(COVARS_PSI, "_z"), COVARS_PSI)
    # Extract covariate values
    cov_psi <- spatial_covs %>%
      select(contains(covars_psi))
    # Add quadratics if needed
    if (!all(is.na(PSI_QUADS))) {
      psi_quads <- names(cov_psi)[str_detect(names(cov_psi), 
                                             paste(PSI_QUADS, collapse = "|"))]
      cov_psi[,paste0(psi_quads, "2")] <- cov_psi[,psi_quads] ^ 2
    }
    # Put columns in alphabetical order and convert to a matrix
    cov_psi <- cov_psi %>%
      select(order(colnames(.))) %>%
      as.matrix
  }
    
# Detection probability (p)  
  if (!all(is.na(COVARS_P))) { 
    # Create vector of standardized covariate names
    covars_p <- ifelse(COVARS_P %in% covs_cont,
                       paste0(COVARS_P, "_z"), COVARS_P)
    # Extract covariate values
    cov_p <- surveys %>%
      select(contains(covars_p))
    # Add quadratics if needed
    if (!all(is.na(P_QUADS))) {
      p_quads <- names(cov_p)[str_detect(names(cov_p), 
                                         paste(P_QUADS, collapse = "|"))]
      cov_p[,paste0(p_quads, "2")] <- cov_p[,p_quads] ^ 2
    }
    # Put columns in alphabetical order and convert to a matrix
    cov_p <- cov_p %>%
      select(order(colnames(.))) %>%
      as.matrix  
  }

# Extinction probability (eps)  
  if (!all(is.na(COVARS_EPS))) {
    # Create vector of standardized covariate names
    covars_eps <- ifelse(COVARS_EPS %in% covs_cont,
                         paste0(COVARS_EPS, "_z"), COVARS_EPS)  
    # Extract covariate values  
    cov_eps <- sitetrans %>%
      select(contains(covars_eps))
    # Add quadratics if needed
    if (!all(is.na(EPS_QUADS))) {
      eps_quads <- names(cov_eps)[str_detect(names(cov_eps), 
                                             paste(EPS_QUADS, collapse = "|"))]
      cov_eps[,paste0(eps_quads, "2")] <- cov_eps[,eps_quads] ^ 2
    }
    # Put columns in alphabetical order
    cov_eps <- cov_eps %>%
      select(order(colnames(.)))
    # Add interactions if needed
    if (N_EPS_INTERACTS > 0) {
      eps_interacts <- matrix(NA, nrow = 1, ncol = 2)
      for (i in 1:N_EPS_INTERACTS) {
        eps_interacts[i,] <- get(paste0("EPS_INT", i))
        eps_interacts[i,1] <- ifelse(eps_interacts[i,1] %in% covs_cont,
                                     paste0(eps_interacts[i,1], "_z"),
                                     eps_interacts[i,1])
        eps_interacts[i,2] <- ifelse(eps_interacts[i,2] %in% covs_cont,
                                     paste0(eps_interacts[i,2], "_z"),
                                     eps_interacts[i,2])    
        cov_eps[,paste0(eps_interacts[i,], collapse = "_")] <- 
          cov_eps[,eps_interacts[i,1]] * cov_eps[,eps_interacts[i,2]]
      }
    }
    # Convert to a matrix
    cov_eps <- as.matrix(cov_eps)
  }
  
# Colonization probability (gam)  
  if (!all(is.na(COVARS_GAM))) {   
    # Create vector of standardized covariate names
    covars_gam <- ifelse(COVARS_GAM %in% covs_cont,
                         paste0(COVARS_GAM, "_z"), COVARS_GAM)  
    # Extract covariate values  
    cov_gam <- sitetrans %>%
      select(contains(covars_gam))
    # Add quadratics if needed
    if (!all(is.na(GAM_QUADS))) {
      gam_quads <- names(cov_gam)[str_detect(names(cov_gam), 
                                             paste(GAM_QUADS, collapse = "|"))]
      cov_gam[,paste0(gam_quads, "2")] <- cov_gam[,gam_quads] ^ 2
    }
    # Put columns in alphabetical order
    cov_gam <- cov_gam %>%
      select(order(colnames(.)))
    # Add interactions if needed
    if (N_GAM_INTERACTS > 0) {
      gam_interacts <- matrix(NA, nrow = 1, ncol = 2)
      for (i in 1:N_GAM_INTERACTS) {
        gam_interacts[i,] <- get(paste0("GAM_INT", i))
        gam_interacts[i,1] <- ifelse(gam_interacts[i,1] %in% covs_cont,
                                     paste0(gam_interacts[i,1], "_z"),
                                     gam_interacts[i,1])
        gam_interacts[i,2] <- ifelse(gam_interacts[i,2] %in% covs_cont,
                                     paste0(gam_interacts[i,2], "_z"),
                                     gam_interacts[i,2])    
        cov_gam[,paste0(gam_interacts[i,], collapse = "_")] <- 
          cov_gam[,gam_interacts[i,1]] * cov_gam[,gam_interacts[i,2]]
      }
    }
    # Convert to a matrix
    cov_gam <- as.matrix(cov_gam)
  }
  
# Calculate the number of covariates for each parameter  
if (exists("cov_psi")) {n_cov_psi <- ncol(cov_psi)}
if (exists("cov_p")) {n_cov_p <- ncol(cov_p)}
if (exists("cov_eps")) {n_cov_eps <- ncol(cov_eps)}
if (exists("cov_gam")) {n_cov_gam <- ncol(cov_gam)}

#------------------------------------------------------------------------------#
# Package things up for JAGS
#------------------------------------------------------------------------------#

# z (latent occupancy for each site & season) is what we're interested in
# In JAGS, z will be stored in a matrix (n_sites[i] * n_seasons[t])

# In order to create initial values for z, we'll need to summarize detections 
# over occasions at each site in each season. To do this:
  # Create an array with detection data (row = site, col = occ, slice = season)
  # Then use an apply function, summarizing over columns

n_sites <- max(surveys$site_index)
n_seasons <- max(surveys$season_index)
max_n_occasions <- max(surveys$occ_index)
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

# Bundle data for JAGS that will be needed in any model, regardless of covariate 
# structure
jags_data <- list(y = surveys$det,
                  n_sites = n_sites,
                  n_seasons = n_seasons,
                  n_obs = nrow(surveys),
                  n_sitetrans = nrow(sitetrans),
                  site = surveys$site_index,
                  season = surveys$season_index,
                  site_ec = sitetrans$site_index,
                  trans_ec = sitetrans$trans_index,
                  z = known_state_occ(y_array))

# Vector of parameters to monitor (needed for all models)
params <- c("mean_psi", "beta_psi0", "mean_p", "beta_p0",
            "mean_eps", "beta_eps0", "mean_gam", "beta_gam0", "PAO")

# List of initial values (needed for all models)
inits_list <- list(mean_psi = runif(1, 0, 1),
                   mean_p = runif(1, 0, 1),
                   mean_eps = runif(1, 0, 1),
                   mean_gam = runif(1, 0, 1),
                   z = inits_state_occ(y_array))

# Add in data, parameters, inits associated with covariates
if (exists("cov_psi")) {
  jags_data <- c(jags_data,
                 list(cov_psi = cov_psi),
                 list(n_cov_psi = n_cov_psi))
  params <- c(params, "beta_psi")
  inits_list <- c(inits_list,
                  list(beta_psi = runif(n_cov_psi, -2, 2)))
}
if (exists("cov_p")) {
  jags_data <- c(jags_data,
                 list(cov_p = cov_p),
                 list(n_cov_p = n_cov_p))
  params <- c(params, "beta_p")
  inits_list <- c(inits_list,
                  list(beta_p = runif(n_cov_p, -2, 2)))
}
if (exists("cov_eps")) {
  jags_data <- c(jags_data,
                 list(cov_eps = cov_eps),
                 list(n_cov_eps = n_cov_eps))
  params <- c(params, "beta_eps")
  inits_list <- c(inits_list,
                  list(beta_eps = runif(n_cov_eps, -2, 2)))
}
if (exists("cov_gam")) {
  jags_data <- c(jags_data,
                 list(cov_gam = cov_gam),
                 list(n_cov_gam = n_cov_gam))
  params <- c(params, "beta_gam")
  inits_list <- c(inits_list,
                  list(beta_gam = runif(n_cov_gam, -2, 2)))
}

inits <- function(){inits_list}

# Identify correct JAGS model
jags_model <- "JAGS/JAGS_MS"
if (exists("cov_psi")) {jags_model <- paste0(jags_model, "_", "psi")}
if (exists("cov_p")) {jags_model <- paste0(jags_model, "_", "p")}
if (exists("cov_eps")) {jags_model <- paste0(jags_model, "_", "eps")}
if (exists("cov_gam")) {jags_model <- paste0(jags_model, "_", "gam")}
jags_model <- paste0(jags_model, ".txt")
  
#------------------------------------------------------------------------------#
# Run model in JAGS
#------------------------------------------------------------------------------#

nc <- 3      # Number of chains
na <- 5000   # Number of iterations to run in the adaptive phase
nb <- 10000  # Number of iterations to discard (burn-in)
ni <- 30000  # Number of iterations per chain (including burn-in)
nt <- 20     # Thinning rate
# Note: with this number of iterations, the model can take many minutes to run

out <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = jags_model,
            n.chains = nc,
            n.adapt = na,
            n.burnin = nb,
            n.iter = ni,
            n.thin = nt,
            parallel = TRUE)

# print(out)

# Trace and density plots
# MCMCtrace(out,pdf = FALSE)
