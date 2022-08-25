################################################################################
# SODN -- Camera trap data, 2016-2022
# Get gridMET climate data

# ER Zylstra
# 2022-08-25
################################################################################

library(downloader)
library(terra)

# Load one of the 30-m DEMs to use for CRS
dem_chir <- rast("data/covariates/CHIR_DEM_1as.tif")

# Download gridMET climate data (daily values, at ~4 km resolution across US)
# Annual file is multi-layer raster where each layer represents a single day 
# in given year.

# pr = precipitation (mm)
# tmmn = min near-surface air temperature
# tmmx = max near-surface air temperature

# Select years
yrs <- 2016:2022

# Select type of climate data 
data_types <- c("pr")

# gridMET website with data
url_base <- "http://www.northwestknowledge.net/metdata/data/"

for (yr in yrs) {
  for (data_type in data_types) {
    
    url <- paste0(url_base, data_type, "_", yr, ".nc")
    destfile <- paste0("data/covariates/", data_type, "_", yr, ".nc")
    
    # Download data
    downloader::download(url = url, destfile = destfile, mode = "wb")
    
    # Load raster (eg, object name = pr2016)
    assign(paste0(data_type, yr), 
           rast(paste0("data/covariates/", data_type, "_", yr, ".nc")))
    
    # Crop SpatRaster
    assign(paste0(data_type, yr), 
           crop(get(paste0(data_type, yr)), ext(-113.2, -109.0, 31.7, 32.5)))
    
    # Most layers from NPS are in NAD83 (I think), so best to convert this 
    # raster, which is in lon/lat WGS 84 (difference should be trivial)
    assign(paste0(data_type, yr), 
           terra::project(get(paste0(data_type, yr)), crs(dem_chir)))
    
    # Remove original netCDF file 
    file.remove(destfile)
    
    # Save cropped raster to file
    writeCDF(get(paste0(data_type, yr)), 
             filename = paste0("data/covariates/", data_type, yr, ".nc"))
    
  }
}



