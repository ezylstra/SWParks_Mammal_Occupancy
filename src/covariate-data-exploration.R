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

# Extract fire data as dataframe
fires_df <- as.data.frame(fires)
dim(fires_df) #202 polygons

# There are a lot of duplicated columns, removing things that seem unnecessary
fires_df <- fires_df[,1:43]

count(fires_df, UNQE_FIRE_, INCIDENT, FIRE_YEAR)
  # Looks like 27 fires/incidents
# Look at a couple fires in particular:
horseshoe2 <- subset(fires, fires$INCIDENT == "Horseshoe 2")
madrone <- subset(fires, fires$INCIDENT == "Madrone")

plot(chir_b, lty = 2)
plot(subset(parks, parks$UNIT_CODE == "CHIR"), add = TRUE)
plot(madrone, col = rgb(170, 210, 240, alpha = 0.4*255, max = 255), 
     border = NA, add = TRUE)
plot(horseshoe2, col = rgb(140, 240, 140, alpha = 0.2*255, max = 255), 
     border = NA, add = TRUE)
