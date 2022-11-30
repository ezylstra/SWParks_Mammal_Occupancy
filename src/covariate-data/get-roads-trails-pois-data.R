################################################################################
# SODN -- Camera trap data, 2016-2022
# Get vector files with roads, trails, points-of-interest (POIs)

# ER Zylstra
# 2022-09-02
################################################################################

library(dplyr)
library(terra)
library(tigris)
library(httr)
library(geojsonsf)
library(sf)

# Load park boundaries
parks <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")

# Load park boundaries + 3-km buffer
parks_b <- vect("data/covariates/shapefiles/Boundaries_wBuffer_3parks.shp")

#------------------------------------------------------------------------------#
# Roads
#------------------------------------------------------------------------------#

# Load 2021 road layers using tigris package
roads_co <- roads(state = "Arizona", county = "Cochise", year = 2021)
roads_pi <- roads(state = "Arizona", county = "Pima", year = 2021)

# Convert to SpatVectors
roads_co <- vect(roads_co)
roads_pi <- vect(roads_pi)

# Crop roads to buffered park boundaries
roads_chir <- crop(roads_co, subset(parks_b, parks$UNIT_CODE == "CHIR"))
roads_orpi <- crop(roads_pi, subset(parks_b, parks$UNIT_CODE == "ORPI"))
roads_sagw <- crop(roads_pi, subset(parks_b, parks$UNIT_CODE == "SAGW"))

# Feature class codes (MTFCC)
  # S1100 (n = 4): Primary road, limited access highways 
  # S1200 (92): Secondary road, main arteries that are not limited access
  # S1400 (14350): Local neighborhood road, rural road, city street (paved)
  # S1500 (291): 4WD, unpaved dirt trail where 4WD vehicle is required
  
  # There are a limited number of features in the S16XX - S18XX classes, 
  # but they're either well outside of park boundaries or irrelevant for 
  # describing distributions of mammal species.  Will remove these features
  road_features <- c("S1100", "S1200", "S1400", "S1500")
  roads_chir <- subset(roads_chir, roads_chir$MTFCC %in% road_features)
  roads_orpi <- subset(roads_orpi, roads_orpi$MTFCC %in% road_features)
  roads_sagw <- subset(roads_sagw, roads_sagw$MTFCC %in% road_features)

# Write to file
# writeVector(roads_chir, filename = "data/covariates/shapefiles/roads_chir.shp")
# writeVector(roads_orpi, filename = "data/covariates/shapefiles/roads_orpi.shp")
# writeVector(roads_sagw, filename = "data/covariates/shapefiles/roads_sagw.shp")

#------------------------------------------------------------------------------#
# Trails
#------------------------------------------------------------------------------#

# Extracting layer from NPS website
base_url <- httr::parse_url("https://mapservices.nps.gov/arcgis/rest/services/NationalDatasets/NPS_Public_Trails/FeatureServer/0")
base_url$path <- paste(base_url$path, "query", sep = "/")
base_url$query <- list(where = "UNITCODE LIKE 'SAGU' OR UNITCODE LIKE 'CHIR' OR UNITCODE LIKE 'ORPI'", 
                       outFields = "*",
                       f = "geojson")
request <- httr::build_url(base_url)
tmp <- tempfile()
download.file(request, tmp)
trails_sf <- sf::st_read(tmp, drivers = "GeoJSON")
trails <- terra::vect(trails_sf)
# Note: this file includes trails in SAGE, but that's probably fine for now

# Reproject to NAD83
trails <- terra::project(trails, crs(parks))

# Write to file
# writeVector(trails, 
#             filename = "data/covariates/shapefiles/trails.shp", 
#             overwrite = TRUE)

#------------------------------------------------------------------------------#
# Points of interest
#------------------------------------------------------------------------------#

# Extracting layer from NPS website
base_url <- httr::parse_url("https://mapservices.nps.gov/arcgis/rest/services/NationalDatasets/NPS_Public_POIs/FeatureServer/0")
base_url$path <- paste(base_url$path, "query", sep = "/")
base_url$query <- list(where = "UNITCODE LIKE 'SAGU' OR UNITCODE LIKE 'CHIR' OR UNITCODE LIKE 'ORPI'", 
                       outFields = "*",
                       f = "geojson")
request <- httr::build_url(base_url)
tmp <- tempfile()
download.file(request, tmp)
pois_sf <- sf::st_read(tmp, drivers = "GeoJSON")
pois <- terra::vect(pois_sf)

# Reproject to NAD83
pois <- terra::project(pois, crs(parks))

# What are POIs?
count(as.data.frame(pois), POITYPE)
  # Amphitheater (1) # Might overlap with buildings layer
  # Bench (14)
  # Building (1) # Should overlap with buildings layer
  # Campground (9)
  # Campsite (25)
  # Headquarters (1) # Should overlap with buildings layer
  # Interpretive Exhibit (40)
  # Mile Marker (139)
  # Parking lot (68) # Should overlap with parking lots layer
  # Picnic Area (14)
  # Picnic Table (16)
  # Restroom (6)
  # Trailhead (40)
  # Viewpoint (2)
  # Visitor Center (4) # Should overlap with buildings layer
  
  # I think it's okay if there are locations that are duplicated in the POIs, 
  # buildings, parking lot layers since we'll combine them all into 
  # a distance-to-feature layer. However, the Mile Markers seem irrelevant.
  # Will delete them from the POIs layer
  pois <- subset(pois, pois$POITYPE != "Mile Marker")

# Write to file
# writeVector(pois, 
#             filename = "data/covariates/shapefiles/POIs.shp",
#             overwrite = TRUE)

#------------------------------------------------------------------------------#
# Buildings (polygons)
#------------------------------------------------------------------------------#

# Extracting layer from NPS website
base_url <- httr::parse_url("https://mapservices.nps.gov/arcgis/rest/services/NationalDatasets/NPS_Public_Buildings/FeatureServer/2")
base_url$path <- paste(base_url$path, "query", sep = "/")
base_url$query <- list(where = "UNITCODE LIKE 'SAGU' OR UNITCODE LIKE 'CHIR' OR UNITCODE LIKE 'ORPI'", 
                       outFields = "*",
                       f = "geojson")
request <- httr::build_url(base_url)
tmp <- tempfile()
download.file(request, tmp)
builds_sf <- geojsonsf::geojson_sf(tmp)
builds <- terra::vect(builds_sf)

# Reproject to NAD83
builds <- terra::project(builds, crs(parks))

# Write to file
# writeVector(builds, 
#             filename = "data/covariates/shapefiles/buildings.shp",
#             overwrite = TRUE)

#------------------------------------------------------------------------------#
#Parking lots (polygons)
#------------------------------------------------------------------------------#

# Extracting layer from NPS website
base_url <- httr::parse_url("https://mapservices.nps.gov/arcgis/rest/services/NationalDatasets/NPS_Public_ParkingLots/FeatureServer/2")
base_url$path <- paste(base_url$path, "query", sep = "/")
base_url$query <- list(where = "UNITCODE LIKE 'SAGU' OR UNITCODE LIKE 'CHIR' OR UNITCODE LIKE 'ORPI'", 
                       outFields = "*",
                       f = "geojson")
request <- httr::build_url(base_url)
tmp <- tempfile()
download.file(request, tmp)
lots_sf <- geojsonsf::geojson_sf(tmp)
lots <- terra::vect(lots_sf)

# Reproject to NAD83
lots <- terra::project(lots, crs(parks))

# Write to file
# writeVector(lots, 
#             filename = "data/covariates/shapefiles/parking_lots.shp",
#             overwrite = TRUE)

