################################################################################
# SODN -- Camera trap data, 2016-2022
# Get gridMET climate data and create covariate rasters

# ER Zylstra
# 2022-08-25
################################################################################

library(downloader)
library(terra)
library(dplyr)
library(stringr)
library(ncdf4)

# Load shapefile with park boundaries
parks <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")

# Download gridMET climate data (daily values, at ~4 km resolution across US)
# Annual file is multi-layer raster where each layer represents a single day 
# in given year.

# pr = precipitation (mm)
# tmmn = min near-surface air temperature
# tmmx = max near-surface air temperature
# vpd = vapor pressure deficit (Pa)

# Select years
yrs <- 2016:2025

# Select type of climate data 
data_types <- c("pr", "vpd")

# gridMET website with data
url_base <- "http://www.northwestknowledge.net/metdata/data/"

# Logical indicating whether we want to re-download gridMET data if done 
# previously
redownload <- FALSE

# Logical indicating whether we want to recreate covariate rasters if they
# already exist
replace <- FALSE

#------------------------------------------------------------------------------#
# Download gridMET data
#------------------------------------------------------------------------------#

# Extract files with original data from zip folder, if one exists:
weather_orig_folder <- "data/covariates/weather-orig-rasters/"
zipfile_orig <- "data/covariates/weather-orig.zip"
if (file.exists(zipfile_orig)) {
  unzip(zipfile = zipfile_orig)
}

for (yr in yrs) {
  for (data_type in data_types) {
    
    url <- paste0(url_base, data_type, "_", yr, ".nc")
    destfile <- paste0(weather_orig_folder, data_type, yr, ".nc")
    
    if(redownload == FALSE & file.exists(destfile)) {next}
    
    # Download data
    downloader::download(url = url, destfile = destfile, mode = "wb")
    
    # Load raster (eg, object name = pr2016)
    assign(paste0(data_type, yr), rast(destfile))
    
    # Crop SpatRaster
    assign(paste0(data_type, yr), 
           crop(get(paste0(data_type, yr)), ext(-113.2, -109.0, 31.7, 32.5)))
    
    # Most layers from NPS are in NAD83, so best to convert this raster, which 
    # is in lon/lat WGS 84 (difference should be trivial)
    assign(paste0(data_type, yr), 
           terra::project(get(paste0(data_type, yr)), crs(parks)))
    
    # Remove original netCDF file 
    invisible(file.remove(destfile))
    
    # Save cropped raster to file
    writeCDF(get(paste0(data_type, yr)), 
             filename = paste0(weather_orig_folder, data_type, yr, ".nc"),
             overwrite = TRUE)
    
  }
}

weather_orig_files <- list.files(path = weather_orig_folder,
                                 pattern = "*.nc",
                                 full.names = TRUE)

# Create a zip archive of all the .nc files (first removing previous archive)
if (file.exists(zipfile_orig)) {
  invisible(file.remove(zipfile_orig))
}
zip(zipfile = zipfile_orig,
    files = weather_orig_files)

#------------------------------------------------------------------------------#
# Derive precipitation covariates
#------------------------------------------------------------------------------#

# Extract files with processed data from zip folder, if one exists:
weather_derived_folder <- "data/covariates/weather-derived-rasters/"
zipfile_derived <- "data/covariates/weather-derived.zip"
if (file.exists(zipfile_derived)) {
  unzip(zipfile = zipfile_derived)
}
weather_derived_files <- list.files(weather_derived_folder,
                                    full.names = TRUE)

# Identify years that we have precipitation data for
pr_files <- str_subset(weather_orig_files, pattern = "pr")
pr_yrs <- basename(pr_files)
pr_yrs <- sort(as.numeric(str_sub(pr_yrs, 3, 6)))

# Is last year of data complete?
# Calculate how many days of data we have and use this to determine the last
# year we can calculate each preciptiation variable
last_yr <- rast(paste0(weather_orig_folder, "/pr", max(pr_yrs), ".nc"))
last_day <- nlyr(last_yr)
rm(last_yr)

# Load original precipitation data
for (yr in pr_yrs) {
  # Create SpatRasters titled "prYYYY"
  assign(paste0("pr", yr), 
         rast(paste0(weather_orig_folder, "/pr", yr, ".nc")))
}

# Create annual rasters with cumulative precipitation during monsoon season
# Raster/filenames will be: monsoon_ppt_YEAR
  
  # Remove last year from list if we don't have data through Sep 30
  yrs <- pr_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-09-30"))) {
    yrs <- yrs[-length(yrs)]
  }
  
  for (yr in yrs) {
    
    new_name <- paste0("monsoon_ppt_", yr)
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    monsoon_ppt <- get(paste0("pr", yr))
    startd  <- lubridate::yday(paste0(yr, "-06-15"))
    endd <- lubridate::yday(paste0(yr, "-09-30"))
    monsoon_ppt <- monsoon_ppt[[startd:endd]]
    monsoon_ppt <- app(monsoon_ppt, fun = sum)
    assign(new_name, monsoon_ppt)
  } 

# Create annual rasters with cumulative precipitation during winter (Oct-Mar)
# Raster/filenames will be: winter_ppt_YEARYEAR

  # Remove last year from list if we don't have data through Mar 30
  yrs <- pr_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-03-30"))) {
    yrs <- yrs[-length(yrs)]
  }
  start_yrs <- yrs[-length(yrs)]
  end_yrs <- start_yrs + 1
  both_yrs <- paste0(start_yrs, end_yrs)

  for (i in 1:length(both_yrs)) {
    
    new_name <- paste0("winter_ppt_", both_yrs[i])
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    # Create rasters with Oct-Dec precipitation
    start1d <- lubridate::yday(paste0(start_yrs[i], "-10-01"))
    end1d <- lubridate::yday(paste0(start_yrs[i], "-12-31"))
    winter1_ppt <- get(paste0("pr", start_yrs[i]))
    winter1_ppt <- winter1_ppt[[start1d:end1d]]
    
    # Create rasters with Jan-Mar precipitation
    start2d <- lubridate::yday(paste0(end_yrs[i], "-01-01"))
    end2d <- lubridate::yday(paste0(end_yrs[i], "-03-31"))
    winter2_ppt <- get(paste0("pr", end_yrs[i]))
    winter2_ppt <- winter2_ppt[[start2d:end2d]]  
    
    # Merge rasters from two calendar years
    winter_ppt <- c(winter1_ppt, winter2_ppt)
    winter_ppt <- app(winter_ppt, fun = sum)
    assign(new_name, winter_ppt) 
  }
  
# Create rasters with cumulative precipitation during 10 months prior to sampling. 

# For SAGW, want precip for Mar-Dec (sampling Jan-Feb)
# Raster/filenames will be: SAGW_MarDec_ppt_YEAR
  
  # Remove last year from list if we don't have data through Dec 31
  yrs <- pr_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-12-31"))) {
    yrs <- yrs[-length(yrs)]
  }

  for (yr in yrs) {
    
    new_name <- paste0("SAGW_MarDec_ppt_", yr)
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    MarDec_ppt <- get(paste0("pr", yr))
    MarDec_ppt <- terra::crop(x = MarDec_ppt, 
                              y = subset(parks, parks$UNIT_CODE == "SAGW"),
                              snap = "out")
    startd  <- lubridate::yday(paste0(yr, "-03-01"))
    endd <- lubridate::yday(paste0(yr, "-12-31"))
    MarDec_ppt <- MarDec_ppt[[startd:endd]]
    MarDec_ppt <- app(MarDec_ppt, fun = sum)
    assign(new_name, MarDec_ppt)
  } 

# For ORPI, want precip for May-Feb (sampling Mar-Apr)
# Raster/filenames will be: ORPI_MayFeb_ppt_YEARYEAR
  
  # Remove last year from list if we don't have data through end of Feb
  yrs <- pr_yrs
  if (last_day < (lubridate::yday(paste0(max(yrs), "-03-01")) - 1)) {
    yrs <- yrs[-length(yrs)]
  } 
  start_yrs <- yrs[-length(yrs)]
  end_yrs <- start_yrs + 1
  both_yrs <- paste0(start_yrs, end_yrs)
  
  for (i in 1:length(both_yrs)) {
    
    new_name <- paste0("ORPI_MayFeb_ppt_", both_yrs[i])
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    # Create rasters with May-Dec precipitation
    start1d <- lubridate::yday(paste0(start_yrs[i], "-05-01"))
    end1d <- lubridate::yday(paste0(start_yrs[i], "-12-31"))
    ppt1 <- get(paste0("pr", start_yrs[i]))
    ppt1 <- terra::crop(x = ppt1, 
                        y = subset(parks, parks$UNIT_CODE == "ORPI"),
                        snap = "out")
    ppt1 <- ppt1[[start1d:end1d]]
    
    # Create rasters with Jan-Feb precipitation
    start2d <- lubridate::yday(paste0(end_yrs[i], "-01-01"))
    end2d <- lubridate::yday(paste0(end_yrs[i], "-03-01")) - 1
    ppt2 <- get(paste0("pr", end_yrs[i]))
    ppt2 <- terra::crop(x = ppt2, 
                        y = subset(parks, parks$UNIT_CODE == "ORPI"),
                        snap = "out")
    ppt2 <- ppt2[[start2d:end2d]]  
    
    # Merge rasters from two calendar years
    MayFeb_ppt <- c(ppt1, ppt2)
    MayFeb_ppt <- app(MayFeb_ppt, fun = sum)
    assign(new_name, MayFeb_ppt)  
  }  

# For CHIR, want precip for Jul-Apr (sampling May-Jun from 2021 on)
# Raster/filenames will be: CHIR_JulApr_ppt_YEARYEAR
  
  # Remove years prior to 2020
  # Remove last year from list if we don't have data through end of Feb
  yrs <- pr_yrs[pr_yrs > 2019]
  if (last_day < (lubridate::yday(paste0(max(yrs), "-04-30")) - 1)) {
    yrs <- yrs[-length(yrs)]
  } 
  start_yrs <- yrs[-length(yrs)]
  end_yrs <- start_yrs + 1
  both_yrs <- paste0(start_yrs, end_yrs)  
  
  for (i in 1:length(both_yrs)) {
    
    new_name <- paste0("CHIR_JulApr_ppt_", both_yrs[i])
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    # Create rasters with Jul-Dec precipitation
    start1d <- lubridate::yday(paste0(start_yrs[i], "-07-01"))
    end1d <- lubridate::yday(paste0(start_yrs[i], "-12-31"))
    ppt1 <- get(paste0("pr", start_yrs[i]))
    ppt1 <- terra::crop(x = ppt1, 
                        y = subset(parks, parks$UNIT_CODE == "CHIR"),
                        snap = "out")
    ppt1 <- ppt1[[start1d:end1d]]
    
    # Create rasters with Jan-Apr precipitation
    start2d <- lubridate::yday(paste0(end_yrs[i], "-01-01"))
    end2d <- lubridate::yday(paste0(end_yrs[i], "-04-30"))
    ppt2 <- get(paste0("pr", end_yrs[i]))
    ppt2 <- terra::crop(x = ppt2, 
                        y = subset(parks, parks$UNIT_CODE == "CHIR"),
                        snap = "out")
    ppt2 <- ppt2[[start2d:end2d]]  

    # Merge rasters from two calendar years
    JulApr_ppt <- c(ppt1, ppt2)
    JulApr_ppt <- app(JulApr_ppt, fun = sum)
    assign(new_name, JulApr_ppt)  
  }    

  
  # Create rasters with cumulative precipitation during 6 months prior to sampling. 
  # SAGW only for now
  
  # For SAGW, want precip for Jul-Dec (sampling Jan-Feb)
  # Raster/filenames will be: SAGW_JulDec_ppt_YEAR
  
  # Remove last year from list if we don't have data through Dec 31
  yrs <- pr_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-12-31"))) {
    yrs <- yrs[-length(yrs)]
  }
  
  for (yr in yrs) {
    
    new_name <- paste0("SAGW_JulDec_ppt_", yr)
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    JulDec_ppt <- get(paste0("pr", yr))
    JulDec_ppt <- terra::crop(x = JulDec_ppt, 
                              y = subset(parks, parks$UNIT_CODE == "SAGW"),
                              snap = "out")
    startd  <- lubridate::yday(paste0(yr, "-07-01"))
    endd <- lubridate::yday(paste0(yr, "-12-31"))
    JulDec_ppt <- JulDec_ppt[[startd:endd]]
    JulDec_ppt <- app(JulDec_ppt, fun = sum)
    assign(new_name, JulDec_ppt)
  } 
  
  
#------------------------------------------------------------------------------#
# Derive vapor pressure deficit covariates
#------------------------------------------------------------------------------#
  
  # Identify years that we have vapor pressure deficit data for
  vpd_files <- str_subset(weather_orig_files, pattern = "vpd")
  vpd_yrs <- basename(vpd_files)
  vpd_yrs <- sort(as.numeric(str_sub(vpd_yrs, 4, 7)))
  
  # Is last year of data complete?
  # Calculate how many days of data we have and use this to determine the last
  # year we can calculate each vapor pressure deficit variable
  last_yr <- rast(paste0(weather_orig_folder, "/vpd", max(vpd_yrs), ".nc"))
  last_day <- nlyr(last_yr)
  rm(last_yr)
  
  # Load original vapor pressure deficit data
  for (yr in vpd_yrs) {
    # Create SpatRasters titled "vpdYYYY"
    assign(paste0("vpd", yr), 
           rast(paste0(weather_orig_folder, "/vpd", yr, ".nc")))
  }
  
  # Create annual rasters with cumulative vapor pressure deficit during monsoon season
  # Raster/filenames will be: monsoon_vpd_YEAR
  
  # Remove last year from list if we don't have data through Sep 30
  yrs <- vpd_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-09-30"))) {
    yrs <- yrs[-length(yrs)]
  }
  
  for (yr in yrs) {
    
    new_name <- paste0("monsoon_vpd_", yr)
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    monsoon_vpd <- get(paste0("vpd", yr))
    startd  <- lubridate::yday(paste0(yr, "-06-15"))
    endd <- lubridate::yday(paste0(yr, "-09-30"))
    monsoon_vpd <- monsoon_vpd[[startd:endd]]
    monsoon_vpd <- app(monsoon_vpd, fun = sum)
    assign(new_name, monsoon_vpd)
  } 
  
  # Create annual rasters with cumulative vapor pressure deficit during winter (Oct-Mar)
  # Raster/filenames will be: winter_vpd_YEARYEAR
  
  # Remove last year from list if we don't have data through Mar 30
  yrs <- vpd_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-03-30"))) {
    yrs <- yrs[-length(yrs)]
  }
  start_yrs <- yrs[-length(yrs)]
  end_yrs <- start_yrs + 1
  both_yrs <- paste0(start_yrs, end_yrs)
  
  for (i in 1:length(both_yrs)) {
    
    new_name <- paste0("winter_vpd_", both_yrs[i])
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    # Create rasters with Oct-Dec vpd
    start1d <- lubridate::yday(paste0(start_yrs[i], "-10-01"))
    end1d <- lubridate::yday(paste0(start_yrs[i], "-12-31"))
    winter1_vpd <- get(paste0("vpd", start_yrs[i]))
    winter1_vpd <- winter1_vpd[[start1d:end1d]]
    
    # Create rasters with Jan-Mar vpd
    start2d <- lubridate::yday(paste0(end_yrs[i], "-01-01"))
    end2d <- lubridate::yday(paste0(end_yrs[i], "-03-31"))
    winter2_vpd <- get(paste0("vpd", end_yrs[i]))
    winter2_vpd <- winter2_vpd[[start2d:end2d]]  
    
    # Merge rasters from two calendar years
    winter_vpd <- c(winter1_vpd, winter2_vpd)
    winter_vpd <- app(winter_vpd, fun = sum)
    assign(new_name, winter_vpd) 
  }
  
  # Create rasters with cumulative vpd during 10 months prior to sampling. 
  
  # For SAGW, want precip for Mar-Dec (sampling Jan-Feb)
  # Raster/filenames will be: SAGW_MarDec_vpd_YEAR
  
  # Remove last year from list if we don't have data through Dec 31
  yrs <- vpd_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-12-31"))) {
    yrs <- yrs[-length(yrs)]
  }
  
  for (yr in yrs) {
    
    new_name <- paste0("SAGW_MarDec_vpd_", yr)
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    MarDec_vpd <- get(paste0("vpd", yr))
    MarDec_vpd <- terra::crop(x = MarDec_vpd, 
                              y = subset(parks, parks$UNIT_CODE == "SAGW"),
                              snap = "out")
    startd  <- lubridate::yday(paste0(yr, "-03-01"))
    endd <- lubridate::yday(paste0(yr, "-12-31"))
    MarDec_vpd <- MarDec_vpd[[startd:endd]]
    MarDec_vpd <- app(MarDec_vpd, fun = sum)
    assign(new_name, MarDec_vpd)
  } 
  
  # For ORPI, want vpd for May-Feb (sampling Mar-Apr)
  # Raster/filenames will be: ORPI_MayFeb_vpd_YEARYEAR
  
  # Remove last year from list if we don't have data through end of Feb
  yrs <- vpd_yrs
  if (last_day < (lubridate::yday(paste0(max(yrs), "-03-01")) - 1)) {
    yrs <- yrs[-length(yrs)]
  } 
  start_yrs <- yrs[-length(yrs)]
  end_yrs <- start_yrs + 1
  both_yrs <- paste0(start_yrs, end_yrs)
  
  for (i in 1:length(both_yrs)) {
    
    new_name <- paste0("ORPI_MayFeb_vpd_", both_yrs[i])
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    # Create rasters with May-Dec vpd
    start1d <- lubridate::yday(paste0(start_yrs[i], "-05-01"))
    end1d <- lubridate::yday(paste0(start_yrs[i], "-12-31"))
    vpd1 <- get(paste0("vpd", start_yrs[i]))
    vpd1 <- terra::crop(x = vpd1, 
                        y = subset(parks, parks$UNIT_CODE == "ORPI"),
                        snap = "out")
    vpd1 <- vpd1[[start1d:end1d]]
    
    # Create rasters with Jan-Feb vpd
    start2d <- lubridate::yday(paste0(end_yrs[i], "-01-01"))
    end2d <- lubridate::yday(paste0(end_yrs[i], "-03-01")) - 1
    vpd2 <- get(paste0("vpd", end_yrs[i]))
    vpd2 <- terra::crop(x = vpd2, 
                        y = subset(parks, parks$UNIT_CODE == "ORPI"),
                        snap = "out")
    vpd2 <- vpd2[[start2d:end2d]]  
    
    # Merge rasters from two calendar years
    MayFeb_vpd <- c(vpd1, vpd2)
    MayFeb_vpd <- app(MayFeb_vpd, fun = sum)
    assign(new_name, MayFeb_vpd)  
  }  
  
  # For CHIR, want vpd for Jul-Apr (sampling May-Jun from 2021 on)
  # Raster/filenames will be: CHIR_JulApr_vpd_YEARYEAR
  
  # Remove years prior to 2020
  # Remove last year from list if we don't have data through end of Feb
  yrs <- vpd_yrs[vpd_yrs > 2019]
  if (last_day < (lubridate::yday(paste0(max(yrs), "-04-30")) - 1)) {
    yrs <- yrs[-length(yrs)]
  } 
  start_yrs <- yrs[-length(yrs)]
  end_yrs <- start_yrs + 1
  both_yrs <- paste0(start_yrs, end_yrs)  
  
  for (i in 1:length(both_yrs)) {
    
    new_name <- paste0("CHIR_JulApr_vpd_", both_yrs[i])
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    # Create rasters with Jul-Dec vpr
    start1d <- lubridate::yday(paste0(start_yrs[i], "-07-01"))
    end1d <- lubridate::yday(paste0(start_yrs[i], "-12-31"))
    vpd1 <- get(paste0("vpd", start_yrs[i]))
    vpd1 <- terra::crop(x = vpd1, 
                        y = subset(parks, parks$UNIT_CODE == "CHIR"),
                        snap = "out")
    vpd1 <- vpd1[[start1d:end1d]]
    
    # Create rasters with Jan-Apr vpr
    start2d <- lubridate::yday(paste0(end_yrs[i], "-01-01"))
    end2d <- lubridate::yday(paste0(end_yrs[i], "-04-30"))
    vpd2 <- get(paste0("vpd", end_yrs[i]))
    vpd2 <- terra::crop(x = vpd2, 
                        y = subset(parks, parks$UNIT_CODE == "CHIR"),
                        snap = "out")
    vpd2 <- vpd2[[start2d:end2d]]  
    
    # Merge rasters from two calendar years
    JulApr_vpd <- c(vpd1, vpd2)
    JulApr_vpd <- app(JulApr_vpd, fun = sum)
    assign(new_name, JulApr_vpd)  
  }   
  
  # Create rasters with cumulative vpd during 6 months prior to sampling.
  # SAGW only for now
  
  # For SAGW, want vpd for Jul-Dec (sampling Jan-Feb)
  # Raster/filenames will be: SAGW_MarDec_vpd_YEAR
  
  # Remove last year from list if we don't have data through Dec 31
  yrs <- vpd_yrs
  if (last_day < lubridate::yday(paste0(max(yrs), "-12-31"))) {
    yrs <- yrs[-length(yrs)]
  }
  
  for (yr in yrs) {
    
    new_name <- paste0("SAGW_JulDec_vpd_", yr)
    new_name_full <- paste0(weather_derived_folder, new_name, ".tif")
    if (replace == FALSE & file.exists(new_name_full)) {next}
    
    JulDec_vpd <- get(paste0("vpd", yr))
    JulDec_vpd <- terra::crop(x = JulDec_vpd, 
                              y = subset(parks, parks$UNIT_CODE == "SAGW"),
                              snap = "out")
    startd  <- lubridate::yday(paste0(yr, "-07-01"))
    endd <- lubridate::yday(paste0(yr, "-12-31"))
    JulDec_vpd <- JulDec_vpd[[startd:endd]]
    JulDec_vpd <- app(JulDec_vpd, fun = sum)
    assign(new_name, JulDec_vpd)
  } 
#------------------------------------------------------------------------------#
# Save rasters to file
#------------------------------------------------------------------------------#  
  
weather_raster_names <- c("monsoon_ppt_", "winter_ppt_", "SAGW_MarDec_ppt_",
                          "ORPI_MayFeb_ppt_", "CHIR_JulApr_ppt_", "SAGW_JulDec_ppt_",
                          "monsoon_vpd_", "winter_vpd_", "SAGW_MarDec_vpd_",
                          "ORPI_MayFeb_vpd_", "CHIR_JulApr_vpd_", "SAGW_JulDec_vpd_")
# Create list of all new rasters
weather_rasters <- ls()[grep(paste(weather_raster_names, collapse = "|"), ls())]

# Save to local folder
for (object in weather_rasters) {
  writeRaster(get(object), paste0(weather_derived_folder, object, ".tif"))
}
  
# Create a zip archive of weather rasters 
# Note: we're first removing previous archive!
weather_derived_files <- list.files(path = weather_derived_folder,
                                    pattern = ".tif",
                                    full.names = TRUE)
if (file.exists(zipfile_derived)) {
  invisible(file.remove(zipfile_derived))
}
zip(zipfile = zipfile_derived,
      files = weather_derived_files)

# Remove all weather rasters from local repo
invisible(file.remove(weather_derived_files))
invisible(file.remove(weather_orig_files))
