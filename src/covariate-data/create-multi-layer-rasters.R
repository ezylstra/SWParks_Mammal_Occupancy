################################################################################
# SODN -- Camera trap data, 2016-2022
# Creating multi-layer rasters for each park

# ER Zylstra
# 2023-01-20
################################################################################

library(dplyr)
library(stringr)
library(terra)

# Brief descriptions of each spatial (time-invariant) covariate:
  # boundary = distance to nearest park boundary (m)
  # boundaryUP = distance to nearest (unprotected) park boundary (m; original
    # layer after excluding boundaries adjacent to protected lands)
  # east = eastness (aspect) Ranges from 1 (east-facing) to -1 (west-facing)
  # elev = elevation (m)
  # north = northness (aspect) Ranges from 1 (north-facing) to -1 (south-facing)
  # pois = distance to nearest point of interest (m; incl buildings, lots)
  # roads = distance to nearest road (m)
  # slope = slope (degrees)
  # trail = distance to nearest trail (m)
  # roadbound = distance to nearest road or unprotected boundary (m)
  # trailpoi = distance to nearest trail or point of interest (m)
  # wash [SAGW only] = distance to nearest wash (m) 
  # vegclass [SAGW only] = vegetation class where 1 = low gradient desert with 
    # high-cover mixed cactus and 2-15% tree cover; 2 = Low hillslope and 
    # mountain foothills, rocky, often north facing, cooler, wetter; 
    # 3 = Medium-high gradient, contrasting topography (hilly), often Jojoba 
    # dominant; 4 = developed (no cameras located in this vegclass)
  # burn_severity_2011 [CHIR only] = severity of 2011 burn (integer values, 0:4 
    # with 0 = unburned to 4 = high burn severity)

#------------------------------------------------------------------------------#
# Create multi-layer rasters for each park 
#------------------------------------------------------------------------------#

parks <- c("CHIR", "ORPI", "SAGW")

# Load park boundaries shapefile
boundaries <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")

for (PARK in parks) {
  
  # Extract park boundry
  boundary <- subset(boundaries, boundaries$UNIT_CODE == PARK)

  # Unzip park rasters
  park_folder <- paste0("data/covariates/rasters-", PARK, "/")
  if (PARK == "ORPI") {
    park_zip1 <- paste0("data/covariates/rasters-ORPI-dist.zip") 
    park_zip2 <- paste0("data/covariates/rasters-ORPI-topo.zip") 
    unzip(park_zip1, overwrite = TRUE)
    unzip(park_zip2, overwrite = TRUE)
  } else {
    park_zip <- paste0("data/covariates/rasters-", PARK, ".zip")
    unzip(park_zip)
  }
  
  # Generate list of available rasters
  park_rasters <- list.files(path = park_folder, 
                             pattern = ".tif", 
                             full.names = TRUE)
  # Generate list of (short) names for rasters
  raster_names <- basename(park_rasters) %>%
    str_remove(pattern = ".tif") %>%
    str_remove(pattern = "dist_") %>%
    str_remove(pattern = paste0("_", tolower(PARK)))
  
  # Load rasters into a list
  raster_list <- list()
  for (i in 1:length(park_rasters)) {
    
    raster_list[[i]] <- terra::rast(park_rasters[i])
    names(raster_list[[i]]) <- raster_names[i]

    # Crop raster to park boundary
    raster_list[[i]] <- terra::crop(raster_list[[i]], boundary) 
    
  }
  names(raster_list) <- raster_names

  # For vegclasses, create layers for two dummy variables (classes 2 and 3)
  # [for now, this is just SAGW, but hopefully other parks soon]
  if (PARK == "SAGW") {
    veg_rast <- raster_list[["vegclasses"]]
    veg_rast[veg_rast == 4] <- NA
    vegclass2 <- 1 * (veg_rast == 2)
    names(vegclass2) <- "vegclass2"
    vegclass3 <- 1 * (veg_rast == 3)
    names(vegclass3) <- "vegclass3"
    raster_list <- c(raster_list, vegclass2 = vegclass2, vegclass3 = vegclass3)
  }  

  # For ORPI, boundaryUP layer is slightly smaller than the rest
  if (PARK == "ORPI") {
    raster_list$boundaryUP <- resample(raster_list$boundaryUP,
                                       raster_list$boundary)
  }
  
  # Put layers for each park in same order 
  raster_order <- c("boundary", "boundaryUP", "east", "elev", "north", "pois", 
                    "roads", "roadbound", "slope", "trail", "trailpoi",
                    "burn_severity_2011",
                    "vegclasses", "vegclass2", "vegclass3", "wash")
  raster_order <- raster_order[raster_order %in% names(raster_list)]
  raster_list <- raster_list[raster_order]

  # Create multi-layer raster
  multilayer_raster <- rast(raster_list)

  # Save raster to file
  saveRDS(object = multilayer_raster, 
          file = paste0("data/covariates/spatial-cov-", PARK, ".rds"))
  
  # Remove rasters from folder on local repo
  invisible(file.remove(park_rasters))

} 
