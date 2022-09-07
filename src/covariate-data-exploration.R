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
  # startc <- Sys.time()
  # dist_bound_chir <- distance(dist_bound_chir, chir_line)
  # endc <- Sys.time()
  # endc - startc
  # writeRaster(dist_bound_chir, "data/covariates/dist_boundary_chir.tif")
  
  # orpi_line <- as.lines(subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_bound_orpi <- rast(dem_orpi)
  # dist_bound_orpi <- crop(dist_bound_orpi, orpi_line)
  # starto <- Sys.time()
  # dist_bound_orpi <- distance(dist_bound_orpi, orpi_line)
  # endo <- Sys.time()
  # endo - starto
  # writeRaster(dist_bound_orpi, "data/covariates/dist_boundary_orpi.tif")
  
  # sagw_line <- as.lines(subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_bound_sagw <- rast(dem_sagw)
  # dist_bound_sagw <- crop(dist_bound_sagw, sagw_line)
  # starts <- Sys.time()
  # dist_bound_sagw <- distance(dist_bound_sagw, sagw_line)
  # ends <- Sys.time()
  # ends - starts
  # writeRaster(dist_bound_sagw, "data/covariates/dist_boundary_sagw.tif")

# Load trails
trails <- vect("data/covariates/trails.shp")

# Calculate distance to trails:
  # dist_trail_chir <- rast(dem_chir)
  # dist_trail_chir <- crop(dist_trail_chir, 
  #                         subset(parks, parks$UNIT_CODE == "CHIR"))
  # startc <- Sys.time()
  # dist_trail_chir <- distance(dist_trail_chir, trails)
  # endc <- Sys.time()
  # endc - startc
  # writeRaster(dist_trail_chir, "data/covariates/dist_trail_chir.tif")

  # dist_trail_sagw <- rast(dem_sagw)
  # dist_trail_sagw <- crop(dist_trail_sagw, 
  #                         subset(parks, parks$UNIT_CODE == "SAGW"))
  # starts <- Sys.time()
  # dist_trail_sagw <- distance(dist_trail_sagw, trails)
  # ends <- Sys.time()
  # ends - starts
  # writeRaster(dist_trail_sagw, "data/covariates/dist_trail_sagw.tif")

  # dist_trail_orpi <- rast(dem_orpi)
  # dist_trail_orpi <- crop(dist_trail_orpi, 
  #                         subset(parks, parks$UNIT_CODE == "ORPI"))
  # starto <- Sys.time()
  # dist_trail_orpi <- distance(dist_trail_orpi, trails)
  # endo <- Sys.time()
  # endo - starto
  # writeRaster(dist_trail_orpi, "data/covariates/dist_trail_orpi.tif")  

# Load buildings and parking lots files
buildings <- vect("data/covariates/buildings.shp")
lots <- vect("data/covariates/parking_lots.shp")

# Calculate distance to these features
  # chir_buildings <- subset(buildings, buildings$UNITCODE == "CHIR")
  # chir_lots <- subset(lots, lots$UNITCODE == "CHIR")
  # chir_bl <- union(chir_buildings, chir_lots)
  # chir_bl_lines <- as.lines(chir_bl)
  # dist_buildinglots_chir <- rast(dem_chir)
  # dist_buildinglots_chir <- crop(dist_buildinglots_chir,
  #                                subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_buildinglots_chir <- distance(dist_buildinglots_chir, chir_bl_lines)
  # writeRaster(dist_buildinglots_chir, "data/covariates/dist_buildinglots_chir.tif")
  
  sagw_buildings <- subset(buildings, buildings$UNITCODE == "SAGW")
  sagw_lots <- subset(lots, lots$UNITCODE == "SAGW")
  sagw_bl <- union(sagw_buildings, sagw_lots)
  sagw_bl_lines <- as.lines(sagw_bl)
  dist_buildinglots_sagw <- rast(dem_sagw)
  dist_buildinglots_sagw <- crop(dist_buildinglots_sagw,
                                 subset(parks, parks$UNIT_CODE == "SAGW"))
  starts.bl <- Sys.time()
  dist_buildinglots_sagw <- distance(dist_buildinglots_sagw, sagw_bl_lines)
  ends.bl <- Sys.time()
  ends.bl - starts.bl
  writeRaster(dist_buildinglots_sagw, "data/covariates/dist_buildinglots_sagw.tif")
  rm(dist_buildinglots_sagw)
    
  orpi_buildings <- subset(buildings, buildings$UNITCODE == "ORPI")
  orpi_lots <- subset(lots, lots$UNITCODE == "ORPI")
  orpi_bl <- union(orpi_buildings, orpi_lots)
  orpi_bl_lines <- as.lines(orpi_bl)
  dist_buildinglots_orpi <- rast(dem_orpi)
  dist_buildinglots_orpi <- crop(dist_buildinglots_orpi,
                                 subset(parks, parks$UNIT_CODE == "ORPI"))
  starto.bl <- Sys.time()
  dist_buildinglots_orpi <- distance(dist_buildinglots_orpi, orpi_bl_lines)
  endo.bl <- Sys.time()
  endo.bl - starto.bl
  writeRaster(dist_buildinglots_orpi, "data/covariates/dist_buildinglots_orpi.tif")
  rm(dist_buildinglots_orpi)
  
# Load roads files
roads_chir <- vect("data/covariates/roads_chir.shp")
roads_orpi <- vect("data/covariates/roads_orpi.shp")
roads_sagw <- vect("data/covariates/roads_sagw.shp")

# Calculate distance to roads
  # dist_roads_chir <- rast(dem_chir)
  # dist_roads_chir <- crop(dist_roads_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_roads_chir <- distance(dist_roads_chir, roads_chir)
  # writeRaster(dist_roads_chir, "data/covariates/dist_roads_chir.tif")

  dist_roads_sagw <- rast(dem_sagw)
  dist_roads_sagw <- crop(dist_roads_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  starts.r <- Sys.time()
  dist_roads_sagw <- distance(dist_roads_sagw, roads_sagw)
  ends.r <- Sys.time()
  ends.r - starts.r
  writeRaster(dist_roads_sagw, "data/covariates/dist_roads_sagw.tif")
  rm(dist_roads_sagw)
  
  dist_roads_orpi <- rast(dem_orpi)
  dist_roads_orpi <- crop(dist_roads_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  starto.r <- Sys.time()
  dist_roads_orpi <- distance(dist_roads_orpi, roads_orpi)
  endo.r <- Sys.time()
  endo.r - starto.r
  writeRaster(dist_roads_orpi, "data/covariates/dist_roads_orpi.tif")
  rm(dist_roads_orpi)
  
# Load POI files
pois <- vect("data/covariates/POIs.shp")

# Calculate distance to feature
  # dist_poi_chir <- rast(dem_chir)
  # dist_poi_chir <- crop(dist_poi_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_poi_chir <- distance(dist_poi_chir, pois)
  # writeRaster(dist_poi_chir, "data/covariates/dist_poi_chir.tif")
  
  dist_poi_sagw <- rast(dem_sagw)
  dist_poi_sagw <- crop(dist_poi_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  starts.p <- Sys.time()
  dist_poi_sagw <- distance(dist_poi_sagw, pois)
  ends.p <- Sys.time()
  ends.p - starts.p
  writeRaster(dist_poi_sagw, "data/covariates/dist_poi_sagw.tif")
  rm(dist_poi_sagw)
  
  dist_poi_orpi <- rast(dem_orpi)
  dist_poi_orpi <- crop(dist_poi_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  starto.p <- Sys.time()
  dist_poi_orpi <- distance(dist_poi_orpi, pois)
  endo.p <- Sys.time()
  endo.p - starto.p
  writeRaster(dist_poi_orpi, "data/covariates/dist_poi_orpi.tif")

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
