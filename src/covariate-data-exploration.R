################################################################################
# SODN -- Camera trap data, 2016-2022
# Initial exploration of covariate data

# ER Zylstra
# 2022-08-25
################################################################################

library(dplyr)
library(terra)

# Load park boundaries
parks <- vect("data/covariates/SWNC_nps_boundary_22020630.shp")

# Load park boundaries + 3-km buffer
parks_b <- vect("data/covariates/SWNC_nps_boundary_3km_22020630.shp")

# Extract polygons for the 3 park units of interest: CHIR, ORPI, and SAGW
parks <- subset(parks, parks$UNIT_CODE %in% c("CHIR", "ORPI", "SAGW"))
parks_b <- subset(parks_b, parks_b$UNIT_CODE %in% c("CHIR", "ORPI", "SAGW"))

# Load 30-m DEMs for each park
dem_chir <- rast("data/covariates/CHIR_DEM_1as.tif")
dem_orpi <- rast("data/covariates/ORPI_DEM_1as.tif")
dem_sagw <- rast("data/covariates/SAGW_DEM_1as.tif")

# Calculate slope (in degrees)
slope_chir <- terrain(dem_chir, v = "slope", unit = "degrees")
slope_orpi <- terrain(dem_orpi, v = "slope", unit = "degrees")
slope_sagw <- terrain(dem_sagw, v = "slope", unit = "degrees")

# Calculate eastness and northness (aspect)
east_chir <- sin(terrain(dem_chir, v = "aspect", unit = "radians"))
east_orpi <- sin(terrain(dem_orpi, v = "aspect", unit = "radians"))
east_sagw <- sin(terrain(dem_sagw, v = "aspect", unit = "radians"))
north_chir <- cos(terrain(dem_chir, v = "aspect", unit = "radians"))
north_orpi <- cos(terrain(dem_orpi, v = "aspect", unit = "radians"))
north_sagw <- cos(terrain(dem_sagw, v = "aspect", unit = "radians"))

# Calculate distance to park boundary
# Note: takes forever, even if cropping raster to park boundary
  # chir_line <- as.lines(subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_bound_chir <- rast(dem_chir)
  # dist_bound_chir <- crop(dist_bound_chir, chir_line)
  # dist_bound_chir <- distance(dist_bound_chir, chir_line)
  # writeRaster(dist_bound_chir, "data/covariates/dist_boundary_chir.tif")
  
  # orpi_line <- as.lines(subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_bound_orpi <- rast(dem_orpi)
  # dist_bound_orpi <- crop(dist_bound_orpi, orpi_line)
  # dist_bound_orpi <- distance(dist_bound_orpi, orpi_line)
  # writeRaster(dist_bound_orpi, "data/covariates/dist_boundary_orpi.tif")
  
  # sagw_line <- as.lines(subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_bound_sagw <- rast(dem_sagw)
  # dist_bound_sagw <- crop(dist_bound_sagw, sagw_line)
  # dist_bound_sagw <- distance(dist_bound_sagw, sagw_line)
  # writeRaster(dist_bound_sagw, "data/covariates/dist_boundary_sagw.tif")

# Load roads files
roads_chir <- vect("data/covariates/roads_chir.shp")
roads_orpi <- vect("data/covariates/roads_orpi.shp")
roads_sagw <- vect("data/covariates/roads_sagw.shp")

# Calculate distance to roads
  # dist_roads_chir <- rast(dem_chir)
  # dist_roads_chir <- crop(dist_roads_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_roads_chir <- distance(dist_roads_chir, roads_chir)
  # writeRaster(dist_roads_chir, "data/covariates/dist_roads_chir.tif")
  
  # dist_roads_sagw <- rast(dem_sagw)
  # dist_roads_sagw <- crop(dist_roads_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_roads_sagw <- distance(dist_roads_sagw, roads_sagw)
  # writeRaster(dist_roads_sagw, "data/covariates/dist_roads_sagw.tif")
  
  # dist_roads_orpi <- rast(dem_orpi)
  # dist_roads_orpi <- crop(dist_roads_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_roads_orpi <- distance(dist_roads_orpi, roads_orpi)
  # writeRaster(dist_roads_orpi, "data/covariates/dist_roads_orpi.tif")

# Load trails
trails <- vect("data/covariates/trails.shp")

# Calculate distance to trails:
  # dist_trail_chir <- rast(dem_chir)
  # dist_trail_chir <- crop(dist_trail_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_trail_chir <- distance(dist_trail_chir, trails)
  # writeRaster(dist_trail_chir, "data/covariates/dist_trail_chir.tif")

  # dist_trail_sagw <- rast(dem_sagw)
  # dist_trail_sagw <- crop(dist_trail_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_trail_sagw <- distance(dist_trail_sagw, trails)
  # writeRaster(dist_trail_sagw, "data/covariates/dist_trail_sagw.tif")

  # dist_trail_orpi <- rast(dem_orpi)
  # dist_trail_orpi <- crop(dist_trail_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_trail_orpi <- distance(dist_trail_orpi, trails)
  # writeRaster(dist_trail_orpi, "data/covariates/dist_trail_orpi.tif")  

# Load all point-of-interest files
# Buildings and lots are polygons, POIs are points
buildings <- vect("data/covariates/buildings.shp")
lots <- vect("data/covariates/parking_lots.shp")
pois <- vect("data/covariates/POIs.shp")

# Convert buildings and lots layers to points
buildings_pt <- as.points(buildings)
lots_pt <- as.points(lots)

# Combine all of them
bl <- union(buildings_pt, lots_pt)
allpois <- union(bl, pois)

# Subset by park
allpois_chir <- subset(allpois, allpois$UNITCODE == "CHIR")
allpois_orpi <- subset(allpois, allpois$UNITCODE == "ORPI")
allpois_sagw <- subset(allpois, allpois$UNITCODE == "SAGU")

# Calculate distance to these features
  dist_bl_chir <- rast(dem_chir)
  dist_bl_chir <- crop(dist_bl_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  start <- Sys.time()
  dist_bl_chir <- distance(dist_bl_chir, allpois_chir)
  Sys.time() - start
  writeRaster(dist_bl_chir, "data/covariates/dist_pois_chir.tif")
  
  dist_bl_sagw <- rast(dem_sagw)
  dist_bl_sagw <- crop(dist_bl_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  start <- Sys.time()
  dist_bl_sagw <- distance(dist_bl_sagw, allpois_sagw)
  Sys.time() - start
  writeRaster(dist_bl_sagw, "data/covariates/dist_pois_sagw.tif")

  dist_bl_orpi <- rast(dem_orpi)
  dist_bl_orpi <- crop(dist_bl_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  start <- Sys.time()
  dist_bl_orpi <- distance(dist_bl_orpi, allpois_orpi)
  Sys.time() - start
  writeRaster(dist_bl_orpi, "data/covariates/dist_pois_orpi.tif")  

# Load precipitation data
# (later I can automate this with apply/loops)
pr2016 <- rast("data/covariates/pr2016.nc")
pr2017 <- rast("data/covariates/pr2017.nc")
pr2018 <- rast("data/covariates/pr2018.nc")
pr2019 <- rast("data/covariates/pr2019.nc")
pr2020 <- rast("data/covariates/pr2020.nc")
pr2021 <- rast("data/covariates/pr2021.nc")
pr2022 <- rast("data/covariates/pr2022.nc")

# check:
# plot(pr2016[[365]]) # Precipitation on 31 Dec 2016
# plot(parks, add = TRUE)
# plot(parks_b, lty = 3, add = TRUE)

# Create annual rasters with cumulative precipitation during monsoon season
for (yr in 2016:2021) {
  monsoon_ppt <- get(paste0("pr",yr))
  startd  <- lubridate::yday(paste0(yr, "-06-15"))
  endd <- lubridate::yday(paste0(yr, "-09-30"))
  monsoon_ppt <- monsoon_ppt[[startd:endd]]
  monsoon_ppt <- app(monsoon_ppt, fun = sum)
  # New rasters name: monsoon_ppt_YEAR
  assign(paste0("monsoon_ppt_",yr), monsoon_ppt)
} 

# TODO: save all these new covariate layers to file?

# TODO: Determine what precipitation metrics we want
  # Cold and warm seasons that NPS has used for other projects?
  # Winter precip? But note that cameras deployed early in the year.
  # 30-year norms for annual precipitation and/or seasonal precipitation?
  # (if so, need to download more annual files)

# Load fire perimeter data
# Note: I already cropped fire data to the 3 parks (because the file was large)
# There were no incidents in ORPI and SAGW (in this layer)
fires <- vect("data/covariates/fire_perimeters_chir.shp")

# There are a lot of duplicated columns, removing things that seem unnecessary
fires <- fires[,1:43]

count(fires_df, UNQE_FIRE_, INCIDENT, FIRE_YEAR)
  # Looks like 27 fires/incidents since 1980
# Look at a couple fires in particular:
horseshoe2 <- subset(fires, fires$INCIDENT == "Horseshoe 2")
madrone <- subset(fires, fires$INCIDENT == "Madrone")

plot(horseshoe2, border = NA, col = "gray")
plot(parks, add = TRUE)
# Horseshoe 2 was the most recent fire, and it burned everything in the park.  

# To calculate the number of fires that occurred in each cell from year X
  # Subset fires layer to have FIRE_YEAR >= X
  # rasterize(fires, dem_chir, sum = TRUE)

last_survey_yr <- 2022
min_year <- last_survey_yr - 15
fire_freq <- subset(fires, fires$FIRE_YEAR > min_year)
fire_freq_raster <- rasterize(fire_freq, dem_chir, sum = TRUE)
