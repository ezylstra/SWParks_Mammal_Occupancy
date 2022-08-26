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
plot(pr2016[[365]])
plot(parks, add = TRUE)
plot(parks_b, lty = 3, add = TRUE)

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

# Calculate distance to park boundary (crop raster to park to reduce run time)
# Note: still takes quite a long time for each park
chir_line <- as.lines(subset(parks, parks$UNIT_CODE == "CHIR"))
dist_bound_chir <- rast(dem_chir)
dist_bound_chir <- crop(dist_bound_chir, chir_line)
startc <- Sys.time()
dist_bound_chir <- distance(dist_bound_chir, chir_line)
endc <- Sys.time()
endc - startc

orpi_line <- as.lines(subset(parks, parks$UNIT_CODE == "ORPI"))
dist_bound_orpi <- rast(dem_orpi)
dist_bound_orpi <- crop(dist_bound_orpi, orpi_line)
starto <- Sys.time()
dist_bound_orpi <- distance(dist_bound_orpi, orpi_line)
endo <- Sys.time()
endo - starto

sagw_line <- as.lines(subset(parks, parks$UNIT_CODE == "SAGW"))
dist_bound_sagw <- rast(dem_sagw)
dist_bound_sagw <- crop(dist_bound_sagw, sagw_line)
starts() <- Sys.time()
dist_bound_sagw <- distance(dist_bound_sagw, sagw_line)
ends <- Sys.time()
ends - starts



# TODO with precip data
  # Convert layer name to date: 
    # see https://tmieno2.github.io/R-as-GIS-for-Economists/gridMET.html
  # Sum values across days of interest
  # Cold and warm seasons that NPS has used for other projects?
  # Monsoon precip (15 Jun - 30 Sep)
  # Winter precip? But note that cameras deployed early in the year.
  # For winter we'll need to combine information from 2 years

