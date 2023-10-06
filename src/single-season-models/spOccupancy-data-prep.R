################################################################################
# Prepping detection and covariate data to be used in single-season occupancy 
# models (run using the spOccupancy package)

# In most instances, this will be called from another script (something like:
# src/single-seasons/models/YEAR/spOccupancy-PARK-SPECIES_YEAR.R)

# ER Zylstra
# Updated 2023-10-06
################################################################################

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
  filter(Park == PARK & Species_code %in% SPECIES & yr == YEAR) %>%
  select(loc, obsdate, yr, o_day) %>%
  arrange(loc, obsdate) %>%
  distinct

# Extract information about camera locations in selected park
events_park <- events %>%
  filter(Park == PARK) %>%
  filter(d_yr == YEAR)
locs_park <- locs %>%
  filter(loc %in% events_park$loc) %>%
  select(loc, longitude, latitude) %>%
  rename(long = longitude, lat = latitude)

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
  ddh[rownames(ddh) == obs$loc[i], 
      colnames(ddh) == as.character(obs$o_day[i])] <- 1
}
# checks:
# sum(ddh == 1, na.rm = TRUE)
# sum(obs$o_day %in% occ_days)

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

# Change values in effort matrix to NA wherever there are NA values in detection
# history matrix (effort values should be 0 in these instances)
effort[which(is.na(dh))] <- NA

# Create standardized version of effort variable
effort_mn <- mean(effort, na.rm = TRUE)
effort_sd <- sd(effort, na.rm = TRUE)
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
                                 events_park$start_day[i]:events_park$end_day[i])
  }
}
events_park <- cbind(events_park, occ_matrix)
# Remove rows in the dataframe that aren't associated with any sampling occasion
events_park <- events_park %>%
  mutate(n_occ = rowSums(select(., paste0(YEAR, "_", 1:ncol(occ_matrix))))) %>%
  filter(n_occ > 0) %>%
  arrange(loc, start_day)

# Check if a sampling occasion at a given camera location spanned two 
# deployments (ie, a camera was immediately redeployed during a sampling occ)
# If so, use the deployment experience value from the 2nd deployment
redeploys <- events_park %>%
  group_by(loc) %>%
  summarize(across(starts_with(paste0(YEAR, "_")), sum)) %>%
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

# Create a matrix with deployment experience values [most experienced person in 
# group] (0 = novice; 1 = experienced; 2 = expert)
deploy_exp <- matrix(NA, nrow = nrow(dh), ncol = ncol(dh))
rownames(deploy_exp) <- rownames(dh) 
colnames(deploy_exp) <- paste0(YEAR, "_", 1:ncol(occ_matrix))
for (i in 1:nrow(deploy_exp)) {
  for (j in 1:ncol(deploy_exp)) {
    if (is.na(dh[i, j])) next
    deploy_exp[i, j] <- 
      events_park$deploy_exp[events_park$loc == rownames(deploy_exp)[i] &
                               events_park[,colnames(deploy_exp)[j]] == 1]
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

day <- matrix(occasions$mid_yday, 
              nrow = nrow(dh), 
              ncol = ncol(dh), 
              byrow = TRUE)

# Change values in day matrix to NA wherever there are NA values in detection
# history matrix
day[which(is.na(dh))] <- NA

# Create standardized value of day variable
day_mn <- mean(day, na.rm = TRUE)
day_sd <- sd(day, na.rm = TRUE)
day_z <- (day - day_mn)/day_sd

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
covs_cont <- str_subset(covs_cont, 
                        "loc|long|lat|vegclass|burn_severity_2011", 
                        negate = TRUE)

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
# cor_df %>%
#   arrange(desc(corr)) %>%
#   filter(abs(corr) > 0.5)

#------------------------------------------------------------------------------#
# Create data object for spOccupancy package
#------------------------------------------------------------------------------#

# First put covariates that could be used in the detection model in a list
# Elements can be n_sites * n_occasions matrices (for survey covariates)
# or vectors of length n_sites (for spatial covariates)
det_covs <- list(day = day,
                 day_z = day_z,
                 deploy_exp = deploy_exp,
                 effort = effort, 
                 effort_z = effort_z,
                 boundary_z = spatial_covs$boundary_z, 
                 east_z = spatial_covs$east_z,
                 elev_z = spatial_covs$elev_z,
                 north_z = spatial_covs$north_z,
                 pois_z = spatial_covs$pois_z,
                 roads_z = spatial_covs$roads_z,
                 slope_z = spatial_covs$slope_z,
                 trail_z = spatial_covs$trail_z,
                 roadbound_z = spatial_covs$roadbound_z,
                 trailpoi_z = spatial_covs$trailpoi_z)
if (PARK == "CHIR") {
  det_covs <- c(det_covs, 
                list(burn_severity_2011 = spatial_covs$burn_severity_2011))
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
