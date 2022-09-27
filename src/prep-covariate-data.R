################################################################################
# SODN -- Camera trap data, 2016-2022
# Prepping covariate data (and saving rasters to zip archives)

# ER Zylstra
# 2022-08-25
################################################################################

library(dplyr)
library(terra)
library(raster)

# Note: Many parts of this script are commented out to avoid accidentally 
# starting functions that take a long time to run and/or to avoid 
# unnecessarily overwriting files 

# Load park boundaries
parks <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")

# Load park boundaries + 3-km buffer
parks_b <- vect("data/covariates/shapefiles/Boundaries_wBuffer_3parks.shp")

# Folders where we'll temporarily store park-specific rasters
chir_folder <- "data/covariates/rasters-CHIR/"
orpi_folder <- "data/covariates/rasters-ORPI/"
sagw_folder <- "data/covariates/rasters-SAGW/"

# Zip folders where we'll store park-specific rasters
chir_zip <- "data/covariates/rasters-CHIR.zip"
orpi_zip <- "data/covariates/rasters-ORPI.zip"
sagw_zip <- "data/covariates/rasters-SAGW.zip"

# Unzip park-specifc folders and load 30-m DEMs for each park
unzip(zipfile = chir_zip)
dem_chir_file <- paste0(chir_folder, "CHIR_DEM_1as.tif")
dem_chir <- rast(dem_chir_file)

unzip(zipfile = orpi_zip)
dem_orpi_file <- paste0(orpi_folder, "ORPI_DEM_1as.tif")
dem_orpi <- rast(dem_orpi_file)

unzip(zipfile = sagw_zip)
dem_sagw_file <- paste0(sagw_folder, "SAGW_DEM_1as.tif")
dem_sagw <- rast(dem_sagw_file)

# dem_chir <- rast("data/covariates/CHIR_DEM_1as.tif")
# dem_orpi <- rast("data/covariates/ORPI_DEM_1as.tif")
# dem_sagw <- rast("data/covariates/SAGW_DEM_1as.tif")

#------------------------------------------------------------------------------#
# Derived topographic variables
#------------------------------------------------------------------------------#

# Calculate slope (in degrees)
slope_chir <- terra::terrain(dem_chir, v = "slope", unit = "degrees")
slope_orpi <- terra::terrain(dem_orpi, v = "slope", unit = "degrees")
slope_sagw <- terra::terrain(dem_sagw, v = "slope", unit = "degrees")

# Calculate eastness and northness (aspect)
east_chir <- sin(terra::terrain(dem_chir, v = "aspect", unit = "radians"))
east_orpi <- sin(terra::terrain(dem_orpi, v = "aspect", unit = "radians"))
east_sagw <- sin(terra::terrain(dem_sagw, v = "aspect", unit = "radians"))
north_chir <- cos(terra::terrain(dem_chir, v = "aspect", unit = "radians"))
north_orpi <- cos(terra::terrain(dem_orpi, v = "aspect", unit = "radians"))
north_sagw <- cos(terra::terrain(dem_sagw, v = "aspect", unit = "radians"))

# Save rasters to park-specific folders
# writeRaster(slope_chir, paste0(chir_folder, "slope_chir.tif"))
# writeRaster(slope_orpi, paste0(orpi_folder, "slope_orpi.tif"))
# writeRaster(slope_sagw, paste0(sagw_folder, "slope_sagw.tif"))
# writeRaster(east_chir, paste0(chir_folder, "east_chir.tif"))
# writeRaster(east_orpi, paste0(orpi_folder, "east_orpi.tif"))
# writeRaster(east_sagw, paste0(sagw_folder, "east_sagw.tif"))
# writeRaster(north_chir, paste0(chir_folder, "north_chir.tif"))
# writeRaster(north_orpi, paste0(orpi_folder, "north_orpi.tif"))
# writeRaster(north_sagw, paste0(sagw_folder, "north_sagw.tif"))

# topo_files <- list.files(path = "data/covariates/topography-rasters",
#                          pattern = ".tif",
#                          full.names = TRUE)
# 
# # Create a zip archive of slope & aspect files (first removing previous archive)
# topo_zipfile <- "data/covariates/topography.zip"
# if (file.exists(topo_zipfile)) {
#   invisible(file.remove(topo_zipfile))
# }
# zip(zipfile = topo_zipfile,
#     files = topo_files)
# # Remove topo files from local repo
# invisible(file.remove(topo_files))

#------------------------------------------------------------------------------#
# Distance to park boundary
#------------------------------------------------------------------------------#

# Note: this takes forever (hours), even if cropping raster to park boundary
  # chir_line <- as.lines(subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_bound_chir <- rast(dem_chir)
  # dist_bound_chir <- terra::crop(dist_bound_chir, chir_line)
  # dist_bound_chir <- terra::distance(dist_bound_chir, chir_line)
  # writeRaster(dist_bound_chir, paste0(chir_folder, "dist_boundary_chir.tif")) 

  # orpi_line <- as.lines(subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_bound_orpi <- rast(dem_orpi)
  # dist_bound_orpi <- terra::crop(dist_bound_orpi, orpi_line)
  # dist_bound_orpi <- terra::distance(dist_bound_orpi, orpi_line)
  # writeRaster(dist_bound_orpi, paste0(orpi_folder, "dist_boundary_orpi.tif"))

  # sagw_line <- as.lines(subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_bound_sagw <- rast(dem_sagw)
  # dist_bound_sagw <- terra::crop(dist_bound_sagw, sagw_line)
  # dist_bound_sagw <- terra::distance(dist_bound_sagw, sagw_line)
  # writeRaster(dist_bound_sagw, paste0(sagw_folder, "dist_boundary_sagw.tif"))

#------------------------------------------------------------------------------#
# Distance to road
#------------------------------------------------------------------------------#

# Load roads files
roads_chir <- vect("data/covariates/shapefiles/roads_chir.shp")
roads_orpi <- vect("data/covariates/shapefiles/roads_orpi.shp")
roads_sagw <- vect("data/covariates/shapefiles/roads_sagw.shp")

# Calculate distance to roads
  # dist_roads_chir <- rast(dem_chir)
  # dist_roads_chir <- terra::crop(dist_roads_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_roads_chir <- terra::distance(dist_roads_chir, roads_chir)
  # writeRaster(dist_roads_chir, paste0(chir_folder, "dist_roads_chir.tif"))
  
  # dist_roads_sagw <- rast(dem_sagw)
  # dist_roads_sagw <- terra::crop(dist_roads_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_roads_sagw <- terra::distance(dist_roads_sagw, roads_sagw)
  # writeRaster(dist_roads_sagw, paste0(sagw_folder, "dist_roads_sagw.tif"))
  
  # dist_roads_orpi <- rast(dem_orpi)
  # dist_roads_orpi <- terra::crop(dist_roads_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_roads_orpi <- terra::distance(dist_roads_orpi, roads_orpi)
  # writeRaster(dist_roads_orpi, paste0(orpi_folder, "dist_roads_orpi.tif"))

#------------------------------------------------------------------------------#
# Distance to trail
#------------------------------------------------------------------------------#

# Load trails
trails <- vect("data/covariates/shapefiles/trails.shp")

# Calculate distance to trails:
  # dist_trail_chir <- rast(dem_chir)
  # dist_trail_chir <- terra::crop(dist_trail_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_trail_chir <- terra::distance(dist_trail_chir, trails)
  # writeRaster(dist_trail_chir, paste0(chir_folder, "dist_trail_chir.tif"))

  # dist_trail_sagw <- rast(dem_sagw)
  # dist_trail_sagw <- terra::crop(dist_trail_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_trail_sagw <- terra::distance(dist_trail_sagw, trails)
  # writeRaster(dist_trail_sagw, paste0(sagw_folder, "dist_trail_sagw.tif"))

  # dist_trail_orpi <- rast(dem_orpi)
  # dist_trail_orpi <- terra::crop(dist_trail_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_trail_orpi <- terra::distance(dist_trail_orpi, trails)
  # writeRaster(dist_trail_orpi, paste0(orpi_folder, "dist_trail_orpi.tif"))

#------------------------------------------------------------------------------#
# Distance to point-of-interest (POIs, buildings, parking lots combined)
#------------------------------------------------------------------------------#

# Load all point-of-interest files
# Buildings and lots are polygons, POIs are points
buildings <- vect("data/covariates/shapefiles/buildings.shp")
lots <- vect("data/covariates/shapefiles/parking_lots.shp")
pois <- vect("data/covariates/shapefiles/POIs.shp")

# Convert buildings and lots layers to points
buildings_pt <- as.points(buildings)
lots_pt <- as.points(lots)

# Combine all of them
bl <- terra::union(buildings_pt, lots_pt)
allpois <- terra::union(bl, pois)

# Subset by park
allpois_chir <- subset(allpois, allpois$UNITCODE == "CHIR")
allpois_orpi <- subset(allpois, allpois$UNITCODE == "ORPI")
allpois_sagw <- subset(allpois, allpois$UNITCODE == "SAGU")

# Calculate distance to these features
  # dist_pois_chir <- rast(dem_chir)
  # dist_pois_chir <- terra::crop(dist_pois_chir, subset(parks, parks$UNIT_CODE == "CHIR"))
  # dist_pois_chir <- terra::distance(dist_pois_chir, allpois_chir)
  # writeRaster(dist_pois_chir, paste0(chir_folder, "dist_pois_chir.tif"))
  
  # dist_pois_sagw <- rast(dem_sagw)
  # dist_pois_sagw <- terra::crop(dist_pois_sagw, subset(parks, parks$UNIT_CODE == "SAGW"))
  # dist_pois_sagw <- terra::distance(dist_pois_sagw, allpois_sagw)
  # writeRaster(dist_pois_sagw, paste0(sagw_folder, "dist_pois_sagw.tif"))

  # dist_pois_orpi <- rast(dem_orpi)
  # dist_pois_orpi <- terra::crop(dist_pois_orpi, subset(parks, parks$UNIT_CODE == "ORPI"))
  # dist_pois_orpi <- terra::distance(dist_pois_orpi, allpois_orpi)
  # writeRaster(dist_pois_orpi, paste0(orpi_folder, "dist_pois_orpi.tif"))

#------------------------------------------------------------------------------#
# Fire perimeter data
#------------------------------------------------------------------------------#

# We have a shapefile available that delineates fire perimeters in southwestern 
# parks (SWNC_FirePerimiters_3km.shp). There were no fires in ORPI or SAGW, so I 
# created a shapefile with fire perimeters in CHIR (fire_perimeters_chir.shp).

# The code used to explore that is below, but I removed the shapefile from 
# the repo because Horseshoe 2 was the only fire that occurred in the last 15
# years and it burned everything in the park. 

# fires <- vect("data/covariates/shapefiles/fire_perimeters_chir.shp")
# 
# # There are a lot of duplicated columns, removing things that seem unnecessary
# fires <- fires[,1:43]
# 
# count(as.data.frame(fires), UNQE_FIRE_, INCIDENT, FIRE_YEAR)
#   # Looks like 27 fires/incidents since 1980
# # Look at a couple fires in particular:
# horseshoe2 <- subset(fires, fires$INCIDENT == "Horseshoe 2")
# madrone <- subset(fires, fires$INCIDENT == "Madrone")
# plot(horseshoe2, border = NA, col = "gray")
# plot(parks, add = TRUE)
# 
# # To calculate the number of fires that occurred in each cell from year X
#   # Subset fires layer to have FIRE_YEAR >= X
#   # terra::rasterize(fires, dem_chir, sum = TRUE)
# 
# last_survey_yr <- 2022
# min_year <- last_survey_yr - 15
# fire_freq <- subset(fires, fires$FIRE_YEAR > min_year)
# fire_freq_raster <- terra::rasterize(fire_freq, dem_chir, sum = TRUE)
# 
# Nothing to save here, as there's no spatial variation in fire characteristics

#------------------------------------------------------------------------------#
# Burn severity data (Horseshoe 2 fire in CHIR only)
#------------------------------------------------------------------------------#

# Read in original raster (.adf file) that's located outside of repo
# adf_file_path <- ".../chir2011mtbs/w001001.adf"
adf_file_path <- "C:/Users/erin/Documents/SODN/Mammals/covariates/chirburn/chir2011mtbs/w001001.adf"
burn <- rast(adf_file_path)

# Look at cell values
freq(burn, digits = 3) # Integer values (0:5)

# Convert to factor
burn <- as.factor(burn)

# Reproject and crop raster
burn <- terra::project(burn, crs(dem_chir))
burn <- terra::crop(burn, dem_chir)

# Resample, so geometry aligns with DEM raster
burn <- terra::resample(burn, dem_chir, method = "near")

# Write to file
# writeRaster(burn, paste0(chir_folder, "burn_severity_chir2011.tif")) 

#------------------------------------------------------------------------------#
# Vegetation classes (SAGW only)
#------------------------------------------------------------------------------#

veg <- vect("C:/Users/erin/Documents/SODN/Mammals/covariates/SAGW_veg_polys.shp")

# Reproject 
veg <- terra::project(veg, crs(dem_chir))

# Aggregate polygons into 4-class vegetation layer (+ developed)
veg4 <- terra::aggregate(veg, by = "Group4Des")
plot(veg4, "Group4Des", border = NULL, main = NULL)

# Get rid of variables we don't need
veg4$MapClass_C <- NULL
veg4$MapClass <- NULL
veg4$Group5 <- NULL
veg4$Group5Des <- NULL
veg4$agg_n <- NULL

# Create short name for veg classes:
vegclasses <- unique(as.data.frame(veg4[, c("Group4", "Group4Des")]))
vegclasses$VegClass <- c("Desert washes",
                         "Developed",
                         "Low gradient desert",
                         "Rocky foothills",
                         "Medium-high gradient")
vegclasses$VegClassNo <- c(5, 4, 1, 2, 3)
veg4$VegClass <- vegclasses$VegClass[match(veg4$Group4, vegclasses$Group4)]
veg4$VegClassNo <- vegclasses$VegClassNo[match(veg4$Group4, vegclasses$Group4)]
veg4$Group4Des <- NULL
# as.data.frame(veg4)

# Rasterize
veg4_raster <- terra::rasterize(veg4, dem_sagw, field = "VegClassNo")
veg4_raster <- terra::crop(veg4_raster, subset(parks, parks$UNIT_CODE == "SAGW"))
# vegclasses[,c(1,3,4)]

# We want to convert cells with value = 5 (washes) to the value of the majority 
# of its neighbors

# Convert cell values of 5 to NA
veg4_raster[veg4_raster == 5] <- NA
# Convert to a RasterLayer (package raster) to use focal function
veg4_rasterl <- raster(veg4_raster)

# Create a function to find the modal value of non-NA neighboring cells (if
# multiple nodes, it'll pick the first)
mode_noNA <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

# Function to replace the focal value with the mode of a 3x3 window if NA
# i = 5 will look at the middle cell in a 3x3 window (i = 13 if 5x5)
# Note that this does extend veg classes slightly outside of the park 
# boundaries, but that shouldn't be a problem.
fill_na <- function(x, i = 5) { 
  if (is.na(x)[i]) {
    return(mode_noNA(x))
  } else {
    return(round(x[i], 0))
  }
}

veg3_rasterl <- raster::focal(veg4_rasterl, w = matrix(1, nrow = 3, ncol = 3), 
                              fun = fill_na, na.rm = FALSE)
# Need to a couple more times because a few of the washes were wide (moving 
# window only contained cells with NAs)
veg3_rasterl <- raster::focal(veg3_rasterl, w = matrix(1, nrow = 3, ncol = 3), 
                              fun = fill_na, na.rm = FALSE)
veg3_rasterl <- raster::focal(veg3_rasterl, w = matrix(1, nrow = 3, ncol = 3), 
                              fun = fill_na, na.rm = FALSE)

# Convert to SpatRaster
veg3 <- rast(veg3_rasterl)
veg3 <- terra::project(veg3, crs(dem_chir))
# Convert to factor
veg3 <- as.factor(veg3)

# See how the camera locations relate to the veg classes:  
locs <- read.csv("data/mammals/SODN_Wildlife_Locations_XY_Revised_20220502.csv")[,2:9]
locs <- locs %>%
  dplyr::filter(UnitCode == "SAGW") %>%
  dplyr::select(c(MarkerName, POINT_X, POINT_Y)) %>%
  dplyr::rename(x = POINT_X, 
                y = POINT_Y) 

locs_sp <- locs %>%
  dplyr::select(-MarkerName) %>%
  as.matrix %>%
  vect(., crs = crs(dem_chir))

plot(veg3)
plot(parks, add = T)
plot(locs_sp, add = T)

# Veg classes at each camera location (using original classes)
camera_veg4 <- cbind(locs, 
                     VegClass = terra::extract(veg4, locs_sp)[, c("VegClass")])
count(camera_veg4, VegClass)
# 29 in (Low gradient desert with high-cover mixed cactus...)
# 13 in (Low hillslope and mountain foothills, rocky, ....)
# 16 in (Medium to high gradient....)
# 1 camera in a desert wash (#34)
plot(locs_sp[34], cex = 2, col = "yellow", add = TRUE)
# Near VC 
# 1 camera NA (#7)
plot(locs_sp[7], cex = 2, col = "yellow", add = TRUE)
# At boundary south of Picture Rocks road
# Original veg file doesn't line up perfectly with park boundary.  On 
# eastern side there are some missing areas, including where this camera is.

# Veg classes at each camera location (after removing desert washes)
camera_veg3 <- cbind(locs, terra::extract(veg3, locs_sp, ID = FALSE))
count(camera_veg3, label)
vegclasses[,c(1,3,4)]
# 31 in class 1 (Low gradient desert)
# 13 in class 2 (Rocky foothills)
# 16 in class 3 (Medium-high gradient)
# 0 in class 4 (Developed)

# Write to file
# writeRaster(veg3, paste0(sagw_folder, "vegclasses_sagw.tif"))
  # This is a 4-class categorical raster with:
  # 1 = Low gradient desert with high-cover mixed cactus and 2-15% tree cover
  # 2 = Low hillslope and mountain foothills, rocky, often north facing, cooler, wetter
  # 3 = Medium-high gradient, contrasting topography (hilly), often Jojoba dominant
  # 4 = Developed

#------------------------------------------------------------------------------#
# Distance to desert wash (based on veg classes, in SAGW only)
#------------------------------------------------------------------------------#  

#Extract desert washes as polygon layer
washes <- subset(veg4, veg4$VegClass == "Desert washes")

# Calculate distance to washes:
dist_wash_sagw <- rast(dem_sagw)

# Rasterize and crop washes raster
wash_raster <- terra::rasterize(washes, dem_sagw, field = "VegClassNo")
# Cells in wash have value = 5, otherwise NA
wash_raster <- terra::crop(wash_raster, subset(parks, parks$UNIT_CODE == "SAGW"))

# Calculate distance to wash
# If x is a SpatRaster and y is missing, terra::distance() computes the 
# distance, for all cells that are NA in SpatRaster x to the nearest cell that 
# is not NA.
dist_wash_sagw <- terra::distance(wash_raster) 

# Write to file
# writeRaster(dist_wash_sagw, paste0(sagw_folder, "dist_wash_sagw.tif"))

#------------------------------------------------------------------------------#
# Weather data
#------------------------------------------------------------------------------#

# Extract precip rasters from zip file
weather_zipfile <- "data/covariates/weather-orig.zip"
zipfiles <- unzip(zipfile = weather_zipfile, list = TRUE)
precip_files <- zipfiles$Name[grepl("pr", zipfiles$Name)]
unzip(zipfile = weather_zipfile,
      files = precip_files)

# Load precipitation data
yrs <- 2016:2021
for (yr in yrs) {
  # Create SpatRasters titled "prYYYY"
  assign(paste0("pr", yr), 
         rast(paste0("data/covariates/weather-orig-rasters/pr", yr, ".nc")))
}

# check:
# plot(pr2016[[365]]) # Precipitation on 31 Dec 2016
# plot(parks, add = TRUE)
# plot(parks_b, lty = 3, add = TRUE)

# Create annual rasters with cumulative precipitation during monsoon season
for (yr in yrs) {
  monsoon_ppt <- get(paste0("pr",yr))
  startd  <- lubridate::yday(paste0(yr, "-06-15"))
  endd <- lubridate::yday(paste0(yr, "-09-30"))
  monsoon_ppt <- monsoon_ppt[[startd:endd]]
  monsoon_ppt <- app(monsoon_ppt, fun = sum)
  # New rasters name: monsoon_ppt_YEAR
  assign(paste0("monsoon_ppt_",yr), monsoon_ppt)
} 

# Save rasters to weather-derived folder
weather_folder <- "data/covariates/weather-derived-rasters/"
monsoon_rasters <- ls()[grep("monsoon_ppt_", ls())]
for (object in monsoon_rasters) {
  writeRaster(get(object), paste0(weather_folder, object, ".tif"))
}

weather_files <- list.files(path = "data/covariates/weather-derived-rasters",
                            pattern = ".tif",
                            full.names = TRUE)

# Create a zip archive of weather rasters 
# Note: we're first removing previous archive!
weather_derived_zipfile <- "data/covariates/weather-derived.zip"
if (file.exists(weather_derived_zipfile)) {
  invisible(file.remove(weather_derived_zipfile))
}
zip(zipfile = weather_derived_zipfile,
    files = weather_files)
# Remove weather rasters from local repo
invisible(file.remove(weather_files))

# TODO: Determine what other precipitation metrics we want
  # Cold and warm seasons that NPS has used for other projects?
  # Winter precip? But note that cameras deployed early in the year.
  # 30-year norms for annual precipitation and/or seasonal precipitation?
  # (if so, need to download more annual files)

#------------------------------------------------------------------------------#
# Put all park-specific rasters in zip files
#------------------------------------------------------------------------------#

chir_files <- list.files(path = chir_folder, full.names = TRUE)
orpi_files <- list.files(path = orpi_folder, full.names = TRUE)
sagw_files <- list.files(path = sagw_folder, full.names = TRUE)

# Create zip archives of files 
# Can first remove previous archives since we unzipped at the top of script

if (file.exists(chir_zip)) {
  invisible(file.remove(chir_zip))
}
zip(zipfile = chir_zip,
    files = chir_files)

if (file.exists(orpi_zip)) {
  invisible(file.remove(orpi_zip))
}
zip(zipfile = orpi_zip,
    files = orpi_files)

if (file.exists(sagw_zip)) {
  invisible(file.remove(sagw_zip))
}
zip(zipfile = sagw_zip,
    files = sagw_files)

# Remove rasters from local repo
invisible(file.remove(chir_files))
invisible(file.remove(orpi_files))
invisible(file.remove(sagw_files))

#------------------------------------------------------------------------------#
# Drainages
#------------------------------------------------------------------------------#

# There are 14,000+ features in the full dataset and eg, 700+ in CHIR, 1900+ in ORPI
# It's not obvious to me that I can select just the bigger drainages based on
# some kind of feature code.  
# For now, it's just not feasible to calculate/use "distance" to drainage 
# based on this shapefile (look at the density of drainages in ORPI for example)

