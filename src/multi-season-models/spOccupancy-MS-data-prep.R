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
# Updated 2025-07-16 by CL McIntyre
################################################################################

#------------------------------------------------------------------------------#
# Create detection histories for selected park, species, and year
#------------------------------------------------------------------------------#

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv(paste0("data/occasions/occasions-", PARK, ".csv"))

# Extract sampling occasion info for selected park and year
occasions <- occasions %>%
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
  filter(Species_code == SPECIES & yr %in% YEARS) %>%
  dplyr::select(loc, obsdate, o_day) %>%
  arrange(loc, obsdate) %>%
  distinct

# Extract information about camera locations in selected park
events_park <- events %>%
  filter(d_yr %in% YEARS)
locs_park <- locs %>%
  filter(loc %in% events_park$loc) %>%
  dplyr::select(loc, longitude, latitude) %>%
  rename(long = longitude, lat = latitude) %>%
  arrange(loc)

# Extract columns from events matrix that correspond to sampling locations
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
  ddh[rownames(ddh) == obs$loc[i], 
      colnames(ddh) == as.character(obs$o_day[i])] <- 1
}
# checks:
# sum(ddh == 1, na.rm = TRUE); sum(obs$o_day %in% occ_days)

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
                                 events_park$start_day[i]:events_park$end_day[i])
  }
}
events_park <- cbind(events_park, occ_matrix)
# Remove rows in the dataframe that aren't associated with any sampling occasion
events_park <- events_park %>%
  mutate(n_occ = rowSums(dplyr::select(., occasions$yr_occ))) %>%
  filter(n_occ > 0) %>%
  arrange(loc, start_day)

# Check if a sampling occasion at a given camera location spanned two 
# deployments (ie, a camera was immediately redeployed during a sampling occ)
# If so, use the deployment experience value from the 2nd deployment
redeploys <- events_park %>%
  group_by(loc) %>%
  summarize(across(occasions$yr_occ, sum)) %>%
  as.data.frame()
for (i in 1:nrow(redeploys)) {
  for (j in 2:ncol(redeploys)) {
    if (redeploys[i, j] > 1) {
      ndeploys <- sum(events_park$loc == redeploys$loc[i])
      events_park[events_park$loc == redeploys$loc[i],
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
        events_park$deploy_exp[events_park$loc == rownames(deploy_exp)[i] & 
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
         mid_yday = round((start_yday + end_yday) / 2))

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

# Indicator for 2022 and beyond, when different types of cameras were used 
# same camera in 2022-2025
camera <- matrix(rep(c(0, 1), 
                          times = c(sum(YEARS < 2022), sum(YEARS >= 2022))),
                      nrow = dim(dh)[1],
                      ncol = dim(dh)[2],
                      byrow = TRUE)

# Indicator for lens type SAGW & ORPI 
# updated to use lens type in events table (rather than binary assumption)
# starting in 2024 these data are now recorded in the field/spatial data, assumed for other years
# accounts for situations where wrong lens was used

# Create an matrix with lens type values 
# (0 = standard; 1 = sensitive; NA = no sampling occasions that year)
lens <- matrix(NA, nrow = dim(dh)[1], ncol = dim(dh)[2], dimnames = c(dimnames(dh)[1],dimnames(dh)[2]), byrow = TRUE)
for (t in 1:length(YEARS)) {
  YR <- YEARS[t]
  for (i in 1:dim(lens)[1]) {
    if (!is.na(table(events_park$loc == rownames(lens)[i]  &  events_park$d_yr == YR)["TRUE"]))
    lens[i,t] <- events_park$lens[events_park$loc == rownames(lens)[i]  &  events_park$d_yr == YR]
  }
}

# Monthly visitation data (now SAGW only rather than all of SAGU) 
# through April 2025 (updated 6/20/2025)
# since monthly visits to the 2 districts have to be populated by hand,  
# only went back to 2016 (instead of 1979 for all of SAGU)
if (PARK == "SAGW") {
  # Read in data
  monthlyvisits <- read.csv("data/covariates/SAGW_MonthlyVisits_2016-2025.csv")
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

# Monthly traffic data (currently only available for SAGW, through April 2025 (updated 6/20/2025)
if (PARK == "SAGW") {
  # Read in data
  monthlytraffic <- read.csv("data/covariates/SAGW_MonthlyTraffic_1992-2025.csv")
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

# Daily water balance variables (deficit, AET) for SAGW (climate zone 4, coordiantes in SAGW)
# downloaded from https://parkfutures.s3.us-west-2.amazonaws.com/maca-tprh-data/index.html
# more info here: https://www.climateanalyzer.biz/water_balance.html
# Daily deficit and actual evapotranspiration (currently only available for SAGW, through December 2024 (updated 6/24/2025)

# Generate cumulative deficit during 10 months prior to sampling. 
# For SAGW, want deficit for Mar-Dec (sampling Jan-Feb)
if (PARK == "SAGW") {
  # Read in data
  dailyWB <- read.csv("data/covariates/SAGU_04_water_balance_historical.csv")
  # Calculate total deficit for Mar-Dec
  deficit <- dailyWB %>%
    mutate(Date = as.Date(Date)) %>%
    mutate(MonthDay = format(Date, format="%m-%d")) %>%
    mutate(Year = year(Date)) %>%
    filter(MonthDay > "02-29") %>%
    filter(Year >= (min(YEARS) - 1) & Year < (max(YEARS))) %>%  # want data for year Mar-Dec for year prior
    group_by(Year) %>%
    summarize(deficit = sum(Deficit.in), .groups = "keep") %>%
    group_by(Year) %>%
    summarize(deficit = mean(deficit)) %>%
    data.frame()   # 2023 is ridiculously low, by factor of 10?
  deficit10 <- matrix(deficit$deficit, 
                    nrow = dim(dh)[1],
                    ncol = dim(dh)[2],
                    byrow = TRUE)
  # Standardize
  deficit10_mn <- mean(deficit10)
  deficit10_sd <- sd(deficit10)
  deficit10_z <- (deficit10 - deficit10_mn)/deficit10_sd
}

# Generate cumulative AET during 10 months prior to sampling. 
# For SAGW, want AET for Mar-Dec (sampling Jan-Feb)
if (PARK == "SAGW") {
  # Read in data
  dailyWB <- read.csv("data/covariates/SAGU_04_water_balance_historical.csv")
  # Calculate total deficit for Mar-Dec
  aet <- dailyWB %>%
    mutate(Date = as.Date(Date)) %>%
    mutate(MonthDay = format(Date, format="%m-%d")) %>%
    mutate(Year = year(Date)) %>%
    filter(MonthDay > "02-29") %>%
    filter(Year >= (min(YEARS) - 1) & Year < (max(YEARS))) %>%  # want data for year Mar-Dec for year prior
    group_by(Year) %>%
    summarize(aet = sum(AET.in), .groups = "keep") %>%
    group_by(Year) %>%
    summarize(aet = mean(aet)) %>%
    data.frame()  # 2023 is ridiculously low, by factor of 10?
  aet10 <- matrix(aet$aet, 
                      nrow = dim(dh)[1],
                      ncol = dim(dh)[2],
                      byrow = TRUE)
  # Standardize
  aet10_mn <- mean(aet10)
  aet10_sd <- sd(aet10)
  aet10_z <- (aet10 - aet10_mn)/aet10_sd
}

# Soil-adjusted vegetation index from Climate Engine (http://climateengine.org)
# Park-wide average of maximum SAVI for April - March (so overlaps sampling)
if (PARK == "SAGW") {
  # Read in data
  SAVI_max <- read.csv("data/covariates/SAGW_SAVImax_Apr-Mar.csv") %>%
    filter(Year >= min(YEARS) & Year <= max(YEARS)) 
  
  savi <- matrix(SAVI_max$SAVImax, 
                  nrow = dim(dh)[1],
                  ncol = dim(dh)[2],
                  byrow = TRUE)
  # Standardize
  savi_mn <- mean(savi)
  savi_sd <- sd(savi)
  savi_z <- (savi - savi_mn)/savi_sd
}


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

# Load shapefile with park boundary
parks <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
park_b <- terra::subset(parks, parks$UNIT_CODE == PARK)
park_b <- as(park_b, "Spatial")

  # Extract and compile monsoon vapor pressure deficit data
  monsoon_vpd_files <- weather_files[str_detect(weather_files, "monsoon_vpd")]
  # Remove monsoon rasters associated with periods outside the years of interest
  # (monsoon vpd in year x could explain occupancy in year x + 1 since 
  # surveys are done in the first half of the year)
  monsoon_vpd_yrs <- paste0(as.character(YEARS - 1), collapse = "|")
  monsoon_vpd_files <- monsoon_vpd_files[str_detect(monsoon_vpd_files, monsoon_vpd_yrs)]

  # Load each raster and compute the mean value across the park in that year
  monsoon_vpd <- rep(NA, length(monsoon_vpd_files))
  for (i in 1:length(monsoon_vpd_files)) {
    monsoon_vpd_raster <- rast(monsoon_vpd_files[i])
    monsoon_vpd[i] <- exact_extract(monsoon_vpd_raster, park_b, "mean")
  }  
  
  monsoon_vpd <- matrix(monsoon_vpd, 
                        nrow = dim(dh)[1],
                        ncol = dim(dh)[2],
                        byrow = TRUE)
  # Standardize
  monsoon_vpd_mn <- mean(monsoon_vpd)
  monsoon_vpd_sd <- sd(monsoon_vpd)
  monsoon_vpd_z <- (monsoon_vpd - monsoon_vpd_mn) / monsoon_vpd_sd 
  
# Extract and compile 10-month vapor pressure deficit data (10-months prior to survey
# season in each park) 
  if (PARK == "ORPI") {
    vpd10_files <- weather_files[str_detect(weather_files, "ORPI_MayFeb_vpd")]
    vpd10_files <- vpd10_files[str_sub(vpd10_files, -8, -5) %in% as.character(YEARS)]
  }
  if (PARK == "SAGW") {
    vpd10_files <- weather_files[str_detect(weather_files, "SAGW_MarDec_vpd")]
    vpd10_files <- vpd10_files[str_sub(vpd10_files, -8, -5) %in% as.character(YEARS - 1)]    
  }
  
  if (PARK == "CHIR") {
    vpd10_files <- weather_files[str_detect(weather_files, "CHIR_JulApr_vpd")]
    vpd10_files <- vpd10_files[str_sub(vpd10_files, -8, -5) %in% as.character(YEARS - 1)]    
  }

 
# Load each raster and compute the mean value across the park in that year
    vpd10 <- rep(NA, length(vpd10_files))
    for (i in 1:length(vpd10_files)) {
      vpd10_raster <- rast(vpd10_files[i])
      vpd10[i] <- exact_extract(vpd10_raster, park_b, "mean")
    }
    
    vpd10 <- matrix(vpd10, 
                    nrow = dim(dh)[1],
                    ncol = dim(dh)[2],
                    byrow = TRUE)
    # Standardize
    vpd10_mn <- mean(vpd10)
    vpd10_sd <- sd(vpd10)
    vpd10_z <- (vpd10 - vpd10_mn) / vpd10_sd 
    
# Extract and compile 6-month vapor pressure deficit data (6-months prior to survey
# season in each park) 
    # if (PARK == "ORPI") {
    #   vpd10_files <- weather_files[str_detect(weather_files, "ORPI_MayFeb_vpd")]
    #   vpd10_files <- vpd10_files[str_sub(vpd10_files, -8, -5) %in% as.character(YEARS)]
    # }
    if (PARK == "SAGW") {
      vpd6_files <- weather_files[str_detect(weather_files, "SAGW_JulDec_vpd")]
      vpd6_files <- vpd6_files[str_sub(vpd6_files, -8, -5) %in% as.character(YEARS - 1)]    
    }
    
    # if (PARK == "CHIR") {
    #   vpd10_files <- weather_files[str_detect(weather_files, "CHIR_JulApr_vpd")]
    #   vpd10_files <- vpd10_files[str_sub(vpd10_files, -8, -5) %in% as.character(YEARS - 1)]    
    # }
    
    
    # Load each raster and compute the mean value across the park in that year
    vpd6 <- rep(NA, length(vpd6_files))
    for (i in 1:length(vpd6_files)) {
      vpd6_raster <- rast(vpd6_files[i])
      vpd6[i] <- exact_extract(vpd6_raster, park_b, "mean")
    }
    
    vpd6 <- matrix(vpd6, 
                    nrow = dim(dh)[1],
                    ncol = dim(dh)[2],
                    byrow = TRUE)
    # Standardize
    vpd6_mn <- mean(vpd6)
    vpd6_sd <- sd(vpd6)
    vpd6_z <- (vpd6 - vpd6_mn) / vpd6_sd 
    

# Extract and compile monsoon precipitation data
    monsoon_files <- weather_files[str_detect(weather_files, "monsoon_ppt")]
    # Remove monsoon rasters associated with periods outside the years of interest
    # (monsoon rainfall in year x could explain occupancy in year x + 1 since 
    # surveys are done in the first half of the year)
    monsoon_yrs <- paste0(as.character(YEARS - 1), collapse = "|")
    monsoon_files <- monsoon_files[str_detect(monsoon_files, monsoon_yrs)]
    
    # Load each raster and compute the mean value across the park in that year
    monsoon_ppt <- rep(NA, length(monsoon_files))
    for (i in 1:length(monsoon_files)) {
      monsoon_raster <- rast(monsoon_files[i])
      monsoon_ppt[i] <- exact_extract(monsoon_raster, park_b, "mean")
    }  
    
    monsoon_ppt <- matrix(monsoon_ppt, 
                          nrow = dim(dh)[1],
                          ncol = dim(dh)[2],
                          byrow = TRUE)
    # Standardize
    monsoon_ppt_mn <- mean(monsoon_ppt)
    monsoon_ppt_sd <- sd(monsoon_ppt)
    monsoon_ppt_z <- (monsoon_ppt - monsoon_ppt_mn) / monsoon_ppt_sd 
    
# Extract and compile 10-month precipitation data (10-months prior to survey
# season in each park) 
    if (PARK == "ORPI") {
      ppt10_files <- weather_files[str_detect(weather_files, "ORPI_MayFeb_ppt")]
      ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS)]
    }
    if (PARK == "SAGW") {
      ppt10_files <- weather_files[str_detect(weather_files, "SAGW_MarDec_ppt")]
      ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS - 1)]    
    }
    
    if (PARK == "CHIR") {
      ppt10_files <- weather_files[str_detect(weather_files, "CHIR_JulApr_ppt")]
      ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS - 1)]    
    }
    
    
    # Load each raster and compute the mean value across the park in that year
    ppt10 <- rep(NA, length(ppt10_files))
    for (i in 1:length(ppt10_files)) {
      ppt10_raster <- rast(ppt10_files[i])
      ppt10[i] <- exact_extract(ppt10_raster, park_b, "mean")
    }
    
    ppt10 <- matrix(ppt10, 
                    nrow = dim(dh)[1],
                    ncol = dim(dh)[2],
                    byrow = TRUE)
    # Standardize
    ppt10_mn <- mean(ppt10)
    ppt10_sd <- sd(ppt10)
    ppt10_z <- (ppt10 - ppt10_mn) / ppt10_sd 
    
  # Extract and compile 6-month precipitation data (6-months prior to survey
    # season in each park) 
    # if (PARK == "ORPI") {
    #   ppt10_files <- weather_files[str_detect(weather_files, "ORPI_MayFeb_ppt")]
    #   ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS)]
    # }
    if (PARK == "SAGW") {
      ppt6_files <- weather_files[str_detect(weather_files, "SAGW_JulDec_ppt")]
      ppt6_files <- ppt6_files[str_sub(ppt6_files, -8, -5) %in% as.character(YEARS - 1)]    
    }
    
    # if (PARK == "CHIR") {
    #   ppt10_files <- weather_files[str_detect(weather_files, "CHIR_JulApr_ppt")]
    #   ppt10_files <- ppt10_files[str_sub(ppt10_files, -8, -5) %in% as.character(YEARS - 1)]    
    # }
    # 
    
    # Load each raster and compute the mean value across the park in that year
    ppt6 <- rep(NA, length(ppt6_files))
    for (i in 1:length(ppt6_files)) {
      ppt6_raster <- rast(ppt6_files[i])
      ppt6[i] <- exact_extract(ppt6_raster, park_b, "mean")
    }
    
    ppt6 <- matrix(ppt6, 
                    nrow = dim(dh)[1],
                    ncol = dim(dh)[2],
                    byrow = TRUE)
    # Standardize
    ppt6_mn <- mean(ppt6)
    ppt6_sd <- sd(ppt6)
    ppt6_z <- (ppt6 - ppt6_mn) / ppt6_sd 
    
    
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
covs_cont <- str_subset(covs_cont, "loc|long|lat|vegclass|rocksizeclass|rockpctclass", negate = TRUE)

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
  dplyr::select(-c(Var1, Var2)) %>%
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
# annual varying covariates that may or may not vary spatially). Including
# unstandardized covariates in case they're needed for backtransforming if the 
# model object is saved and loaded later.
occ_covs <- list(boundary = spatial_covs$boundary, 
                 east = spatial_covs$east,
                 elev = spatial_covs$elev,
                 north = spatial_covs$north,
                 pois = spatial_covs$pois,
                 roads = spatial_covs$roads,
                 slope = spatial_covs$slope,
                 trail = spatial_covs$trail,
                 roadbound = spatial_covs$roadbound,
                 trailpoi = spatial_covs$trailpoi,
                 boundary_z = spatial_covs$boundary_z, 
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
                 monsoon_ppt = monsoon_ppt,
                 monsoon_ppt_z = monsoon_ppt_z,
                 monsoon_vpd = monsoon_vpd,
                 monsoon_vpd_z = monsoon_vpd_z)
if (PARK == "CHIR") {
  occ_covs <- c(occ_covs, 
                list(burn_severity_2011_z = spatial_covs$burn_severity_2011_z))
}
if (PARK == "ORPI") {
  occ_covs <- c(occ_covs,
                list(ppt10 = ppt10,
                     ppt10_z = ppt10_z))
}
if (PARK == "SAGW") {
  occ_covs <- c(occ_covs,
                list(wash = spatial_covs$wash,
                     wash_z = spatial_covs$wash_z,
                     vegclass2 = spatial_covs$vegclass2,
                     vegclass3 = spatial_covs$vegclass3,
                     rocksizeclass2 = spatial_covs$rocksizeclass2,
                     rocksizeclass3 = spatial_covs$rocksizeclass3,
                     rockpctclass2 = spatial_covs$rockpctclass2,
                     rockpctclass3 = spatial_covs$rockpctclass3, 
                     deficit10 = deficit10,
                     deficit10_z = deficit10_z,
                     aet10 = aet10,
                     aet10_z = aet10_z,
                     ppt10 = ppt10,
                     ppt10_z = ppt10_z,
                     ppt6 = ppt6,
                     ppt6_z = ppt6_z,                     
                     vpd10 = vpd10, 
                     vpd10_z = vpd10_z,
                     vpd6 = vpd6, 
                     vpd6_z = vpd6_z,
                     savi = savi,
                     savi_z = savi_z,
                     visits = visits,
                     visits_z = visits_z,
                     traffic = traffic,
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
                 camera = camera,
                 years = years,
                 years_z = years_z,
                 day = day,
                 day_z = day_z,
                 deploy_exp = deploy_exp,
                 effort = effort,
                 effort_z = effort_z)
if (PARK == "CHIR") {
  det_covs <- c(det_covs, 
                list(burn_severity_2011_z = spatial_covs$burn_severity_2011_z))
}
if (PARK == "SAGW") {
  det_covs <- c(det_covs,
                list(wash_z = spatial_covs$wash_z,
                     vegclass2 = spatial_covs$vegclass2,
                     vegclass3 = spatial_covs$vegclass3,
                     visits_z = visits_z,
                     traffic_z = traffic_z))
}
if (PARK != "CHIR" & max(YEARS) > 2022) {
  det_covs <- c(det_covs, 
                list(lens = lens))  
}

# spOccupancy needs camera locations, but can't take lat/long, so we'll need to 
# reproject coordinates to WGS 84, Zone 12 = epsg:32612 (should work for all parks)
loc_matrix <- as.matrix(locs_park[,c("long", "lat")])
loc_utms <- terra::project(loc_matrix,
                           from = "epsg:4269",
                           to = "epsg:32612")

# Create data object (also a list)
data_list <- list(y = dh,
                  occ.covs = occ_covs,
                  det.covs = det_covs,
                  coords = loc_utms)

