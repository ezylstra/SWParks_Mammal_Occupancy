################################################################################
# Prepping detection and covariate data to be used in multi-season occupancy 
# models (run using the spOccupancy package)

# In most instances, this will be called from another script (something like:
# src/multi-season-models/PARK/spOccupancy-PARK-SPECIES_YEARS.R)

# At the moment, new CSVs with monthly visitors and traffic totals need to be 
# added to the repo each year (replacing previous CSVs and updating the name
# of the CSVs in the "Annual covariates" section below). Would be good to 
# eventually automate this.

# ER Zylstra
# Updated 2023-04-07
################################################################################

#------------------------------------------------------------------------------#
# Create detection histories for selected park, species, and year
#------------------------------------------------------------------------------#

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

# Extract sampling occasion info for selected park and year
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
  rename(loc = StdLocName, long = POINT_X, lat = POINT_Y) %>%
  arrange(loc)

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

# Summarize detection data (dh) and effort during each occasion 
dh <- effort <- array(NA, 
                      dim = c(nrow(ddh), max(occasions$occasion), length(YEARS)),
                      dimnames = list(rownames(ddh),
                                      paste0("occ", 1:max(occasions$occasion)),
                                      YEARS))

for (t in 1:length(YEARS)) {
  occyr <- occasions[occasions$yr == YEARS[t],]
  if (nrow(occyr) == 0) next
  allyr <- ddh[,colnames(ddh) %in% min(occyr$start_day):max(occyr$end_day)]
  nocc <- max(occyr$occasion)
  for (i in 1:nocc) {
    multiday <- allyr[,colnames(allyr) %in% occyr$start_day[i]:occyr$end_day[i]]
    dh[,i,t] <- apply(multiday, 1, paNA)
    effort[,i,t] <- apply(multiday, 1, propNA)
  }  
}

# Change values in effort matrix to NA wherever there are NA values in detection
# history matrix (some effort values are 0 in these instances)
effort[which(is.na(dh))] <- NA

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
  mutate(n_occ = rowSums(select(., occasions$yr_occ))) %>%
  filter(n_occ > 0) %>%
  arrange(StdLocName, d_date)

# Check if a sampling occasion at a given camera location spanned two 
# deployments (ie, a camera was immediately redeployed during a sampling occ)
# If so, use the deployment experience value from the 2nd deployment
redeploys <- events_park %>%
  group_by(StdLocName) %>%
  summarize(across(occasions$yr_occ, sum)) %>%
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

# Create an array with deployment experience values [most experienced person in 
# group] (0 = novice; 1 = experienced; 2 = expert)
deploy_exp <- array(NA, dim = dim(dh), dimnames = dimnames(dh))
for (t in 1:length(YEARS)) {
  YR <- YEARS[t]
  occyr <- occasions[occasions$yr == YR,]
  if (nrow(occyr) == 0) next
  nocc <- max(occyr$occasion)
  for (j in 1:nocc) {
    columnname <- paste0(YR, "_", j)
    for (i in 1:dim(dh)[1]) {
      if (is.na(dh[i,j,t])) next
      deploy_exp[i,j,t] <- 
        events_park$deploy_exp[events_park$StdLocName == rownames(deploy_exp)[i] & 
                                 events_park[,columnname] == 1]
    }
  }
}

# Make sure that values in deploy_exp matrix are NA wherever there are NA values 
# in detection history matrix (this should already be true, but just in case)
deploy_exp[which(is.na(dh))] <- NA

#------------------------------------------------------------------------------#
# Create day-of-year variable (using midpoint of each occasion)
#------------------------------------------------------------------------------#

occasions <- occasions %>%
  mutate(start_yday = yday(start),
         end_yday = yday(end),
         mid_yday = round((start_yday + end_yday)/2))

day <- array(NA, dim = dim(dh), dimnames = dimnames(dh))
for (t in 1:length(YEARS)) {
  occyr <- occasions[occasions$yr == YEARS[t],]
  if (nrow(occyr) == 0) next
  impute <- matrix(rep(occyr$mid_yday, nrow(day)), 
                   ncol = nrow(occyr),
                   byrow = TRUE)
  day[,1:ncol(impute),t] <- impute
}
  
# Change values in day matrix to NA wherever there are NA values in detection
# history matrix
day[which(is.na(dh))] <- NA

#------------------------------------------------------------------------------#
# Format arrays with detection data and survey covariate values for spOccupancy
#------------------------------------------------------------------------------#

# Transpose dimensions:
# Currently, all arrays have dims = sites * occasions * years
# For spOccupancy, we want dims = sites * years * occasions
dh <- aperm(dh, c(1, 3, 2))
effort <- aperm(effort, c(1, 3, 2))
deploy_exp <- aperm(deploy_exp, c(1, 3, 2))
day <- aperm(day, c(1, 3, 2))

# Create standardized version of effort variable
effort_mn <- mean(effort, na.rm = TRUE)
effort_sd <- sd(effort, na.rm = TRUE)
effort_z <- (effort - effort_mn)/effort_sd

# Create standardized value of day variable
day_mn <- mean(day, na.rm = TRUE)
day_sd <- sd(day, na.rm = TRUE)
day_z <- (day - day_mn)/day_sd

#------------------------------------------------------------------------------#
# Annual covariates
#------------------------------------------------------------------------------#

# Year (to be used for trend models)
years <- matrix(YEARS, 
                nrow = dim(dh)[1],
                ncol = dim(dh)[2],
                byrow = TRUE)
# Standardize year
years_mn <- mean(years)
years_sd <- sd(years)
years_z <- (years - years_mn)/years_sd

# Indicator for 2022, when different types of cameras were used 
# (will need to revisit this covariate after 2023 season)
camera_2022 <- matrix(rep(c(0, 1, 0), 
                          times = c(sum(YEARS < 2022), 1, sum(YEARS > 2022))),
                      nrow = dim(dh)[1],
                      ncol = dim(dh)[2],
                      byrow = TRUE)

# Monthly visitation data (currently only available for Saguaro, both districts 
# combined)
if (PARK == "SAGW") {
  # Read in data
  monthlyvisits <- read.csv("data/covariates/SAGU_MonthlyVisits_1979-2022.csv")
  # Identify months when surveys occurred
  surveymonths <- unique(c(month(occasions$start), month(occasions$end)))
  # Calculate the total number of visitors during survey months each year
  vis <- monthlyvisits %>%
    pivot_longer(cols = JAN:DEC, 
                 names_to = "month",
                 names_transform = list(month = str_to_title),
                 values_to = "visitors") %>%
    mutate(mon = match(month, month.abb)) %>%
    filter(Year %in% YEARS & mon %in% surveymonths) %>%
    group_by(Year) %>%
    summarize(visitors = sum(visitors)) %>%
    data.frame()
  visits <- matrix(vis$visitors, 
                   nrow = dim(dh)[1],
                   ncol = dim(dh)[2],
                   byrow = TRUE)
  # Standardize
  visits_mn <- mean(visits)
  visits_sd <- sd(visits)
  visits_z <- (visits - visits_mn)/visits_sd
}

# Monthly traffic data (currently only available for SAGW)
if (PARK == "SAGW") {
  # Read in data
  monthlytraffic <- read.csv("data/covariates/SAGW_MonthlyTraffic_1992-2022.csv")
  # Identify months when surveys occurred
  surveymonths <- unique(c(month(occasions$start), month(occasions$end)))
  # Calculate total traffic (averaged across locations) during survey months 
  # each year
  traff <- monthlytraffic %>%
    pivot_longer(cols = JAN:DEC, 
                 names_to = "month",
                 names_transform = list(month = str_to_title),
                 values_to = "traffic") %>%
    mutate(mon = match(month, month.abb)) %>%
    filter(Year %in% YEARS & mon %in% surveymonths) %>%
    group_by(Year, Loc) %>%
    summarize(traffic = sum(traffic), .groups = "keep") %>%
    group_by(Year) %>%
    summarize(traffic = mean(traffic)) %>%
    data.frame()
  traffic <- matrix(traff$traffic, 
                    nrow = dim(dh)[1],
                    ncol = dim(dh)[2],
                    byrow = TRUE)
  # Standardize
  traffic_mn <- mean(traffic)
  traffic_sd <- sd(traffic)
  traffic_z <- (traffic - traffic_mn)/traffic_sd
}
# Note: traffic and visitors were pretty well correlated for 2017-2020, but 
# not in 2021-2022 (overall r = 0.50)

# Only other annual covariates are weather related

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

# Put lat/longs for camera locations in matrix
loc_matrix <- as.matrix(locs_park[,c("long", "lat")])

# Extract and compile monsoon precipitation data
  monsoon_files <- weather_files[str_detect(weather_files, "monsoon_ppt")]
  # Remove monsoon rasters associated with periods outside the years of interest
  # (monsoon rainfall in year x could explain occupancy in year x + 1 since 
  # surveys are done in the first half of the year)
  monsoon_yrs <- paste0(as.character(YEARS - 1), collapse = "|")
  monsoon_files <- monsoon_files[str_detect(monsoon_files, monsoon_yrs)]

  # Load each raster and compile into a list
  monsoon_list <- list()
  for (i in 1:length(monsoon_files)) {
    monsoon_list[[i]] <- rast(monsoon_files[i])
    names(monsoon_list[[i]]) <- "monsoon_ppt"
  }  
  
  # Extract values 
  monsoon_ppt <- matrix(NA, 
                        nrow = dim(dh)[1], 
                        ncol = dim(dh)[2], 
                        dimnames = dimnames(dh)[1])
  # Check that site names in matrix are in same order as locs_park
  # all.equal(rownames(monsoon_ppt), locs_park$loc)
  
  for (i in 1:length(monsoon_list)) {
    monsoon_ppt[,i] <- terra::extract(x = monsoon_list[[i]],
                             y = locs_park[,c("long", "lat")],
                             ID = FALSE)[, "monsoon_ppt"]
  }

  # Standardize values
  monsoon_ppt_mn <- mean(monsoon_ppt)
  monsoon_ppt_sd <- sd(monsoon_ppt)
  monsoon_ppt_z <- (monsoon_ppt - monsoon_ppt_mn)/monsoon_ppt_sd 
  
# Extract and compile 10-month precipitation data (10-months prior to survey
# season in each park) [Don't have this set up for CHIR yet since there are only
# a couple years when sampling was done in May-June]
  if (PARK == "ORPI") {
    ppt10_files <- weather_files[str_detect(weather_files, "ORPI_MayFeb")]
    ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS)]
  }
  if (PARK == "SAGW") {
    ppt10_files <- weather_files[str_detect(weather_files, "SAGW_MarDec")]
    ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS-1)]    
  }

  # Load each raster and compile into a list
  ppt10_list <- list()
  for (i in 1:length(ppt10_files)) {
    ppt10_list[[i]] <- rast(ppt10_files[i])
    names(ppt10_list[[i]]) <- "ppt10"
  }  
  
  # Extract values 
  ppt10 <- matrix(NA, 
                  nrow = dim(dh)[1], 
                  ncol = dim(dh)[2], 
                  dimnames = dimnames(dh)[1])
  for (i in 1:length(ppt10_list)) {
    ppt10[,i] <- terra::extract(x = ppt10_list[[i]],
                                y = locs_park[,c("long", "lat")],
                                ID = FALSE)[, "ppt10"]
  } 
  
  # Standardize values
  ppt10_mn <- mean(ppt10)
  ppt10_sd <- sd(ppt10)
  ppt10_z <- (ppt10 - ppt10_mn)/ppt10_sd 

# Remove rasters from local repo
invisible(file.remove(list.files(weather_folder, full.names = TRUE)))

#------------------------------------------------------------------------------#
# Spatial covariates (time invariant)
#------------------------------------------------------------------------------#

# Load multi-layer raster with spatial data
park_raster <- readRDS(paste0("data/covariates/spatial-cov-", PARK, ".rds"))

# We have two distance-to-boundary layers, one that applies to the entire park
# boundary and one that applies to boundaries that are adjacent to unprotected
# lands (boundaryUP). For now, we'll remove the original boundary layer.
park_raster <- subset(park_raster, "boundary", negate = TRUE)
names(park_raster)[names(park_raster) == "boundaryUP"] <- "boundary"

# Create dataframe with covariate values for each camera location
spatial_covs <- locs_park
# Ensure the order is the same as what's in the detection history matrix
spatial_covs <- spatial_covs[match(rownames(dh), spatial_covs$loc),]

# Extract covariate values for each camera location
spatial_covs <- cbind(spatial_covs, 
                      terra::extract(x = park_raster, 
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

# Add a column with an index for each site (to use as a non-spatial random 
# effect in our models as a simple way to account for the non-independence
# of data that come from the same site over multiple years). 
spatial_covs$site <- 1:nrow(spatial_covs)

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
# cor_df %>%
#   arrange(desc(corr)) %>%
#   filter(abs(corr) > 0.5)

#------------------------------------------------------------------------------#
# Create data object for spOccupancy package
#------------------------------------------------------------------------------#

# First, put covariates that could be used in the occurrence part of the model 
# in a list. Elements can be vectors with length equal to the number of sites 
# (for spatial covariates), or they can be n_sites * n_years matrices (for 
# annual varying covariates that may or may not vary spatially)
occ_covs <- list(boundary_z = spatial_covs$boundary_z, 
                 east_z = spatial_covs$east_z,
                 elev_z = spatial_covs$elev_z,
                 north_z = spatial_covs$north_z,
                 pois_z = spatial_covs$pois_z,
                 roads_z = spatial_covs$roads_z,
                 slope_z = spatial_covs$slope_z,
                 trail_z = spatial_covs$trail_z,
                 roadbound_z = spatial_covs$roadbound_z,
                 trailpoi_z = spatial_covs$trailpoi_z,
                 site = spatial_covs$site,
                 years = years,
                 years_z = years_z,
                 monsoon_ppt_z = monsoon_ppt_z)
if (PARK == "CHIR") {
  occ_covs <- c(occ_covs, 
                list(burn_severity = spatial_covs$burn_severity_2011))
}
if (PARK == "ORPI") {
  occ_covs <- c(occ_covs,
                list(ppt10_z = ppt10_z))
}
if (PARK == "SAGW") {
  occ_covs <- c(occ_covs,
                list(wash_z = spatial_covs$wash_z,
                     vegclass2 = spatial_covs$vegclass2,
                     vegclass3 = spatial_covs$vegclass3,
                     ppt10_z = ppt10_z,
                     visits_z = visits_z,
                     traffic_z = traffic_z))
}

# Then, put covariates that could be used in the detection part of the model in 
# a list. Elements can be vectors with length equal to the number of sites 
# (for spatial covariates), can be n_sites * n_years matrices (for 
# annual varying covariates that may or may not vary spatially), or can be 
# observation-level covariates that vary over sites, years, and occasions.
det_covs <- list(boundary_z = spatial_covs$boundary_z, 
                 east_z = spatial_covs$east_z,
                 elev_z = spatial_covs$elev_z,
                 north_z = spatial_covs$north_z,
                 pois_z = spatial_covs$pois_z,
                 roads_z = spatial_covs$roads_z,
                 slope_z = spatial_covs$slope_z,
                 trail_z = spatial_covs$trail_z,
                 roadbound_z = spatial_covs$roadbound_z,
                 trailpoi_z = spatial_covs$trailpoi_z,
                 site = spatial_covs$site,
                 camera_2022 = camera_2022,
                 years = years,
                 years_z = years_z,
                 day_z = day_z,
                 deploy_exp = deploy_exp,
                 effort_z = effort_z)
if (PARK == "CHIR") {
  det_covs <- c(det_covs, 
                list(burn_severity = spatial_covs$burn_severity_2011))
}
if (PARK == "SAGW") {
  det_covs <- c(det_covs,
                list(wash_z = spatial_covs$wash_z,
                     vegclass2 = spatial_covs$vegclass2,
                     vegclass3 = spatial_covs$vegclass3,
                     visits_z = visits_z,
                     traffic_z = traffic_z))
}

# spOccupancy can't take lat/long, so we'll need to reproject coordinates to 
# WGS 84, Zone 12 = epsg:32612 (should work for all parks)
loc_utms <- terra::project(loc_matrix,
                           from = "epsg:4269",
                           to = "epsg:32612")

# Create data object (also a list)
data_list <- list(y = dh,
                  occ.covs = occ_covs,
                  det.covs = det_covs,
                  coords = loc_utms)
