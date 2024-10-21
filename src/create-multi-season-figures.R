################################################################################
# Create and save NPS formatted figures for multi-season reports

# ER Zylstra
# Updated 2024-01-11
################################################################################

library(tidyverse)
library(terra)
library(spOccupancy)
library(tidyterra)
library(abind)
require(gridExtra)
require(cowplot)
library(ggspatial)

#------------------------------------------------------------------------------#
# Specify parameters of interest
#------------------------------------------------------------------------------#
# Park, year, and species
PARK <- "SAGW"
YEARS <- 2017:2024
SPECIES <- "CALA"

# Logical indicating whether to create maps with mean occurrence probabilities (with roads/trails)
MAP <- TRUE
# Logical indicating whether to create maps with SD of occurrence probabilities (with roads/trails)
MAP_SD <- FALSE
# Logical indicating whether to create maps with mean occurrence probabilities and raw detections
MAP_DETECT <- TRUE
# Logical indicating whether to create a 4-panel figure with mean and SD
# of occurrence probabilities in first and last year in addition to the single
# panel figures for each parameter. (Only relevant if MAP_SD == TRUE)
FOUR_PANEL <- FALSE
# If creating maps, indicate whether to include lat/long axes labels
LATLONG <- FALSE

# Logical indicating whether to create figures with marginal effects of 
# covariates in the occurrence part of the model
MARG_OCC <- TRUE

# Logical indicating whether to create figures with marginal effects of 
# covariates in the detection part of the model
MARG_DET <- TRUE

# Logical indicating whether to create a figure with naive/estimated occurrence
# over time (including trend, if relevant)
OCC_TIME <- TRUE

# Parameters, for single-panel figures
file_extension1 <- ".png"   # can update to jpg
file_extension2 <- ".pdf"   # keep as pdf
device <- cairo_pdf    # required to embed fonts in pdf
dpi <- 300
width <- 6
height <- 4
units <- "in"

# Parameters, for 4-panel figure (occurrence, SD in first and last year)
width_4 <- 6.5
height_4 <- 8

# At later point, could also create option for a 2-panel figure (1 column, 
# 2 rows) with mean occurrence probabilities in first and last year. 

#------------------------------------------------------------------------------#
# Load/create objects needed for any figure
#------------------------------------------------------------------------------#
# Load functions and raw data
source("src/photo-data/format-mammal-data.R")
source("src/functions.R")

# Load multi-season model object 
modelfile <- paste0("output/multi-season-models/", PARK, "-", min(YEARS), "-",
                    max(YEARS), "-", SPECIES, ".rds")
if (!file.exists(modelfile)) {
  stop("Multi-season model for ", SPECIES, " in ", PARK, ", ", min(YEARS), 
       "-", max(YEARS), 
       ", has not been saved in ouput/multi-season-models folder.")
}
model_list <- readRDS(modelfile)
best <- model_list$model
psi_model <- model_list$psi_model
p_model <- model_list$p_model
data_list <- model_list$data

# Make naive estimate dataframe from data_list
naive <- as.data.frame(data_list$y) %>%
  mutate(loc = row.names(.)) %>% 
  pivot_longer(cols=-loc, names_to = c("year", "occ"), names_sep = 4, values_to = "detect") %>% 
  mutate(occ = str_remove(occ,".occ")) %>% 
  mutate(year=as.numeric(year), occ=as.numeric(occ)) %>%
  mutate(occ.active = ifelse(is.na(detect),0,1)) %>%
  group_by(loc,year) %>%
  mutate(Present = max(detect, na.rm=TRUE), Pct_Present = sum(detect, na.rm=TRUE)/sum(occ.active)) %>%
  mutate(Present = as.character(Present)) %>%
  mutate(Present = ifelse(Present=="1","detected", ifelse(Present=="0","not detected", "no data"))) %>%
  left_join(., locs %>% dplyr::select(loc, longitude, latitude), by = "loc")

naive_spat <- vect(st_as_sf(naive,coords = c("longitude", "latitude"), crs = 4269))
  

# Extract dataframe with covariate values at each camera location
spatial_covs <- data_list$occ.covs[lengths(data_list$occ.covs) == dim(data_list$y)[1]]
spatial_covs <- as.data.frame(spatial_covs)

# Extract names of covariates (with and without "_z" subscripts) from best model
# And for occurrence, extract names of spatial covariates
nonspat_z <- c("years_z", "visits_z", "traffic_z", "monsoon_ppt_z", "ppt10_z")
nonspat <- str_remove(nonspat_z, "_z")
psi_covs_z <- create_cov_list(psi_model)
if (length(psi_covs_z) == 1 & any(psi_covs_z == "1")) {
  psi_covs_z <- character(0)
  psi_covs <- character(0)
  psi_spatcovs_z <- character(0)
  psi_spatcovs <- character(0)
} else {
  psi_covs <- psi_covs_z %>% str_remove_all(pattern = "_z")
  psi_spatcovs_z <- psi_covs_z[!psi_covs_z %in% nonspat_z]
  psi_spatcovs <- psi_covs[!psi_covs %in% nonspat]  
}
p_covs_z <- create_cov_list(p_model)
if (length(p_covs_z) == 1 & any(p_covs_z == "1")) {
  p_covs_z <- character(0)
  p_covs <- character(0)
} else {
  p_covs <- p_covs_z %>% str_remove_all(pattern = "_z")
}

# Create table with parameter estimates
occ_estimates <- parameter_estimates(model = best, 
                                     parameter = "occ",
                                     lower_ci = 0.025,
                                     upper_ci = 0.975)
det_estimates <- parameter_estimates(model = best, 
                                     parameter = "det",
                                     lower_ci = 0.025,
                                     upper_ci = 0.975)
occ_estimates <- occ_estimates %>%
  rename(Covariate = Parameter) %>%
  mutate(Parameter = "Occurrence", .before = "Covariate")
det_estimates <- det_estimates %>%
  rename(Covariate = Parameter) %>%
  mutate(Parameter = "Detection", .before = "Covariate")
estimates <- rbind(occ_estimates, det_estimates)
estimates

# Create basename for output files
base_out <- paste0("output/NPS-figures/multi-season/",
                   PARK, "-", min(YEARS), "-", max(YEARS), "-", SPECIES, "-")

# Create custom theme
# Use Frutiger font - NPS standard - whenever possible. NPS employees can 
# download Fruiter to their NPS computer from 
# https://www.nps.gov/subjects/hfc/nps-typefaces.htm
# (must be on the vpn to do the download)
windowsFonts("Frutiger LT Std 55 Roman" = windowsFont("Frutiger LT Std 55 Roman"))

theme_NPS <- ggplot2::theme_classic() + 
  theme(legend.title = element_text(size = 10, color = "black")) +
  theme(legend.text = element_text(size = 9,color = "black")) +
  theme(axis.title = element_text(size = 10,color = "black")) + 
  theme(axis.text = element_text(size = 9,color = "black")) +
  theme(plot.title = element_text(size = 12, color = "black", hjust = 0.5)) +
  theme(axis.ticks = element_line(color = 'black')) + 
  theme(axis.line = element_line(color = 'black')) +
  theme(plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5)) +
  theme(text = element_text(family = "Frutiger LT Std 55 Roman", face = "plain"))

# Define shape and color for detections
scale_color_detect <- function(...){
  ggplot2:::manual_scale('colour', 
                         values = setNames(c("black","darkgrey","lightgrey"),
                                           c("detected", "no data", "not detected")), 
                         ...)
}

scale_shape_detect <- function(...){
  ggplot2:::manual_scale('shape', 
                         values = setNames(c(19,25,17),
                                           c("detected", "no data", "not detected")), 
                         ...)
}



# Create longer park name for use in plots
park <- ifelse(PARK == "CHIR", "Chiricahua NM",
               ifelse(PARK == "SAGW", "Saguaro NP (TMD)", "Organ Pipe Cactus NM"))

#------------------------------------------------------------------------------#
# Annual occurrence estimates (and trends if "years" is in the model)
#------------------------------------------------------------------------------# 
# This is a plot with annual occurrence estimates that include yearly random 
# effects (black circles, with 95% CIs) and the estimated trend (if "years" is
# in the model). If raw_occ = TRUE, naive occurrence estimates (proportion of 
# sites with a detection) will be included (open circles)

if (OCC_TIME) { 
  
  # Load general information about covariates
  covariates <- read.csv("data/covariates/covariates-MS.csv")
  
  occ_time <- occ_time_plot(model = best, 
                            data_list = data_list,
                            covariate_table = covariates,
                            raw_occ = TRUE,
                            central_meas = mean,
                            lower_ci = 0.025,
                            upper_ci = 0.975)
  
  title <- paste0("Occurrence probability of ",
                  species$Common_name[species$Species_code == SPECIES])
  subtitle <- paste0(park, ", ", min(YEARS), "-", max(YEARS)) 
  
  plot_occtime <- occ_time + 
    theme_NPS +
    theme(legend.title = element_blank()) +
    ggtitle(str_wrap(title, 60), subtitle)
  plotname <- "occ-time"
  
  ggsave(plot_occtime, 
         file = paste0(base_out, plotname, file_extension1),
         dpi = dpi, 
         width = width, 
         height = height, 
         units = units)
  ggsave(plot_occtime, 
         file = paste0(base_out, plotname, file_extension2),
         device = device,
         dpi = dpi, 
         width = width, 
         height = height, 
         units = units)

}

# If yearly random effects were in occurrence model, print here for reference
yrREcols <- grepl("years", colnames(best$beta.star.samples))
apply(best$beta.star.samples[,yrREcols], 2, mean)
# Positive values indicate mean occurrence probabilities were higher than 
# expected based only on fixed effects in the model. Negative values indicate
# values were lower than expected. 

#------------------------------------------------------------------------------#
# Marginal effects of covariates in the occurrence part of the model
# (Creates figures for all continuous covariates)
#------------------------------------------------------------------------------#  
if (MARG_OCC) { 
  
  # Load general information about covariates
  covariates <- read.csv("data/covariates/covariates-MS.csv")
  
  # Identify continuous covariates in occurrence part of the best model
  # Excluding years (trend) since that was covered in section above. 
  psi_continuous <- psi_covs_z[!psi_covs_z %in% c("vegclass2", "vegclass3", "years_z")]
  psi_cont_unique <- unique(psi_continuous)
  psi_n_cont <- length(psi_cont_unique)
  
  # If there are any continuous covariates, create a figure for each:
  if (psi_n_cont > 0) {
    # Loop through each covariate
    for (cov in psi_cont_unique) {
      margplot <- marginal_plot_occ(covariate = cov, 
                                    model = best, 
                                    data_list = data_list,
                                    covariate_table = covariates,
                                    central_meas = mean)
      
      title <- paste0("Occurrence probability of ",
                      species$Common_name[species$Species_code == SPECIES], 
                      " vs. ", tolower(margplot$labels[[1]]))
      subtitle <- paste0(park, ", ", min(YEARS), "-", max(YEARS)) 
      
      plot_marg <- margplot + 
        theme_NPS +
        ggtitle(str_wrap(title, 60), subtitle)
      plotname <- paste0("marg-occ-", str_remove(cov, "_z"))
      
      ggsave(plot_marg, 
             file = paste0(base_out, plotname, file_extension1),
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      ggsave(plot_marg, 
             file = paste0(base_out, plotname, file_extension2),
             device = device,
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      
    } 
  }
}

# If vegetation classes were included as covariates in the model, extract
# occurrence probabilities for each class
if (sum(str_detect(psi_covs, "veg")) > 0) {
  occprobs_veg <- vegclass_estimates(model = best, 
                                     parameter = "occ")
  print(occprobs_veg)
}

# If there are no covariates in the model (ie, a null model), print overall 
# occurrence probability
if (psi_n_cont == 0 & length(psi_covs) == 0) {
  overall_occ <- mean_estimate(model = best, 
                               parameter = "occ",
                               lower_ci = 0.025,
                               upper_ci = 0.975)
  print(overall_occ)
}  

#------------------------------------------------------------------------------#
# Marginal effects of covariates in the detection part of the model
# (Creates figures for all continuous covariates)
#------------------------------------------------------------------------------#  
if (MARG_DET) { 
  
  # Load general information about covariates
  covariates <- read.csv("data/covariates/covariates-MS.csv")
  
  # Identify continuous covariates in detection part of the best model
  p_continuous <- p_covs_z[!p_covs_z %in% c("vegclass2", "vegclass3", 
                                            "camera", "lens")]
  p_cont_unique <- unique(p_continuous)
  p_n_cont <- length(p_cont_unique)
  
  # If there are any continuous covariates, create a figure for each:
  if (p_n_cont > 0) {
    # Loop through each covariate
    for (cov in p_cont_unique) {
      
      # Need a day or effort object with raw values in order to put
      # plot on original scale
      if (cov == "day_z") {
        day <- data_list$det.covs$day
      }
      if (cov == "effort_z") {
        effort <- data_list$det.covs$effort
      }
      
      margplot <- marginal_plot_det(covariate = cov, 
                                    model = best, 
                                    data_list = data_list,
                                    covariate_table = covariates,
                                    central_meas = mean)
      
      title <- paste0("Detection probability of ",
                      species$Common_name[species$Species_code == SPECIES], 
                      " vs. ", tolower(margplot$labels[[1]]))
      subtitle <- paste0(park, ", ", min(YEARS), "-", max(YEARS)) 
      
      plot_marg <- margplot + 
        theme_NPS +
        ggtitle(str_wrap(title, 60), subtitle)
      plotname <- paste0("marg-det-", str_remove(cov, "_z"))
      
      ggsave(plot_marg, 
             file = paste0(base_out, plotname, file_extension1),
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      ggsave(plot_marg, 
             file = paste0(base_out, plotname, file_extension2),
             device = device,
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      
    } 
  }
}

# If camera and/or lens was included as a covariate in the model, extract 
# detection probabilities for each combination of covariate levels
if (sum(str_detect(p_covs, c("camera|lens_2023"))) > 0) {
  detprob_cat <- det_cat_estimates(model = best,
                                   lower_ci = 0.025,
                                   upper_ci = 0.975)
  print(detprob_cat)
}

# If there are no covariates in the model (a null model), print overall 
# detection probability
if (p_n_cont == 0 & length(p_covs) == 0) {
  overall_det <- mean_estimate(model = best, 
                               parameter = "det",
                               lower_ci = 0.025,
                               upper_ci = 0.975)
  print(overall_det)
}  

#------------------------------------------------------------------------------#
# Maps with occurrence probabilities (means, SDs)
#------------------------------------------------------------------------------#
if(MAP | MAP_SD) {
  
  if (length(psi_covs) == 0) {
    stop("No covariates in the model of occurrence. Not creating maps.")
  } 
  
  # Load multi-layer raster with spatial data for park
  park_raster <- readRDS(paste0("data/covariates/spatial-cov-", PARK, ".rds"))
  # We have two distance-to-boundary layers, one that applies to the entire park
  # boundary and one that applies to boundaries that are adjacent to unprotected
  # lands (boundaryUP). For now, we'll remove the original boundary layer.
  park_raster <- subset(park_raster, "boundary", negate = TRUE)
  names(park_raster)[names(park_raster) == "boundaryUP"] <- "boundary"
  
  # Crop and mask by park boundary
  park_boundaries <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
  park_boundary <- subset(park_boundaries, park_boundaries$UNIT_CODE == PARK)
  park_boundary_1km <- buffer(park_boundary, width=1000, singlesided=FALSE)
  
  # Load and clip trails layer to park boundary
  park_trails <- vect("data/covariates/shapefiles/trails.shp")
  # clip to current park
  park_trails <- crop(park_trails, park_boundary)
  
  # Load roads shapefile (within 3km)
  park_roads <- if(PARK=="SAGW") vect("data/covariates/shapefiles/roads_sagw_v2.shp") else vect(paste0("data/covariates/shapefiles/roads_",PARK,"_tigris.shp", sep=""))
  park_roads_1km <- crop(park_roads, park_boundary_1km)
  
  # If there are time-varying covariates (other than year/trend) in the 
  # occurrence part of the model, identify whether we want predictions under 
  # average conditions ("averaged") or under observed conditions in the first 
  # and last year ("observed"). Note that if we're using "averaged" and 
  # years/trend isn't in the model, then predictions from the first and last 
  # year will be very similar (but not identical if we're incorporating random 
  # effects).
  if (any(str_detect(string = psi_covs, 
                     pattern = paste(c("visits", "traffic", "monsoon_ppt", "ppt10"),
                                     collapse = "|")))) {
    ANN_PREDS <- "observed"
  }  
  
  # Generate predicted probabilities (preds_mn and preds_sd rasters)
  source("src/multi-season-models/spOccupancy-MS-predictions.R")
  
  # Create figures
  mn_title <- paste0("Mean occurrence probability of ",
                     species$Common_name[species$Species_code == SPECIES])
  sd_title <- paste0("Standard deviation of occurrence probability of ",
                     species$Common_name[species$Species_code == SPECIES])
  gen_title <- paste0("Occurrence probability of ",
                      species$Common_name[species$Species_code == SPECIES])
  subtitle_fy <- paste0(park, ", ", min(YEARS))
  subtitle_ly <- paste0(park, ", ", max(YEARS))
  
  plot_preds_mn_fy <- ggplot() + 
    geom_spatraster(data = preds_mn_firstyr, mapping = aes(fill = mean_firstyr)) + 
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(fill = '', title = mn_title, subtitle = subtitle_fy) +
    theme_NPS + 
    geom_spatvector(data=park_trails, color="lightgrey", lwd = 0.25, linetype = "longdash") +
    geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
    geom_spatvector(data=park_roads, color="lightgrey", inherit.aes=FALSE, lwd = 0.5) + 
    geom_spatvector(data=park_roads, color="black", inherit.aes=FALSE, lwd = 0.1) + 
    annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
    annotation_scale(location = "br", style="ticks") +
    theme(axis.title = element_blank(),
          axis.line = element_blank())
  plot_preds_sd_fy <- ggplot() + 
    geom_spatraster(data = preds_sd_firstyr, mapping = aes(fill = sd_firstyr)) + 
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(fill = '', title = sd_title, subtitle = subtitle_fy) +
    theme_NPS + 
    geom_spatvector(data=park_trails, color="lightgrey", lwd = 0.25, linetype = "longdash") +
    geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
    geom_spatvector(data=park_roads, color="lightgrey", inherit.aes=FALSE, lwd = 0.5) + 
    geom_spatvector(data=park_roads, color="black", inherit.aes=FALSE, lwd = 0.1) + 
    annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
    annotation_scale(location = "br", style="ticks") +
    theme(axis.title = element_blank(),
          axis.line = element_blank()) 
  plot_preds_mn_ly <- ggplot() + 
    geom_spatraster(data = preds_mn_lastyr, mapping = aes(fill = mean_lastyr)) + 
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(fill = '', title = mn_title, subtitle = subtitle_ly) +
    theme_NPS + 
    geom_spatvector(data=park_trails, color="lightgrey", lwd = 0.25, linetype = "longdash") +
    geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
    geom_spatvector(data=park_roads, color="lightgrey", inherit.aes=FALSE, lwd = 0.5) + 
    geom_spatvector(data=park_roads, color="black", inherit.aes=FALSE, lwd = 0.1) + 
    annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
    annotation_scale(location = "br", style="ticks") +
    theme(axis.title = element_blank(),
          axis.line = element_blank())
  plot_preds_sd_ly <- ggplot() + 
    geom_spatraster(data = preds_sd_lastyr, mapping = aes(fill = sd_lastyr)) + 
    scale_fill_viridis_c(na.value = 'transparent') +
    labs(fill = '', title = sd_title, subtitle = subtitle_ly) +
    theme_NPS + 
    geom_spatvector(data=park_trails, color="lightgrey", lwd = 0.25, linetype = "longdash") +
    geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
    geom_spatvector(data=park_roads, color="lightgrey", inherit.aes=FALSE, lwd = 0.5) + 
    geom_spatvector(data=park_roads, color="black", inherit.aes=FALSE, lwd = 0.1) + 
    annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
    annotation_scale(location = "br", style="ticks") +
    theme(axis.title = element_blank(),
          axis.line = element_blank()) 
  
  if (LATLONG) {
    plot_preds_mn_fy <- plot_preds_mn_fy +
      theme(panel.border = element_rect(color = 'black', fill = NA))
    plot_preds_sd_fy <- plot_preds_sd_fy +
      theme(panel.border = element_rect(color = 'black', fill = NA))
    plot_preds_mn_ly <- plot_preds_mn_ly +
      theme(panel.border = element_rect(color = 'black', fill = NA))
    plot_preds_sd_ly <- plot_preds_sd_ly +
      theme(panel.border = element_rect(color = 'black', fill = NA))
  } else {
    plot_preds_mn_fy <- plot_preds_mn_fy + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_sd_fy <- plot_preds_sd_fy +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_mn_ly <- plot_preds_mn_ly + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_sd_ly <- plot_preds_sd_ly +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
  } 
  
  if (MAP) {
    ggsave(plot_preds_mn_fy, 
           file = paste0(base_out, "map-occ-fy", file_extension1),
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
    ggsave(plot_preds_mn_fy, 
           file = paste0(base_out, "map-occ-fy", file_extension2), 
           device = device,
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
    ggsave(plot_preds_mn_ly, 
           file = paste0(base_out, "map-occ-ly", file_extension1),
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
    ggsave(plot_preds_mn_ly, 
           file = paste0(base_out, "map-occ-ly", file_extension2), 
           device = device,
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
  }
  if (MAP_SD) {
    ggsave(plot_preds_sd_fy, 
           file = paste0(base_out, "map-sd-fy", file_extension1),
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
    ggsave(plot_preds_sd_fy, 
           file = paste0(base_out, "map-sd-fy", file_extension2),
           device = device, 
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
    ggsave(plot_preds_sd_ly, 
           file = paste0(base_out, "map-sd-ly", file_extension1),
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
    ggsave(plot_preds_sd_ly, 
           file = paste0(base_out, "map-sd-ly", file_extension2),
           device = device,
           dpi = dpi, 
           width = width, 
           height = height, 
           units = units)
  }
  
  if (FOUR_PANEL) {
    # For now, assuming we never want lat/long axis labels for 4-panel figure.
    
    if (!MAP) {
      stop("Cannot create 4-panel figure with MAP set to FALSE")
    } 
    if (!MAP_SD) {
      stop("Cannot create 4-panel figure with MAP_SD set to FALSE")
    }
    
    # Find min, max values so color scales are consistent for both years
    minmax_mn <- cbind(minmax(preds_mn_firstyr), minmax(preds_mn_lastyr))
    min_mn <- min(minmax_mn[1, ])
    max_mn <- max(minmax_mn[2, ])
    minmax_sd <- cbind(minmax(preds_sd_firstyr), minmax(preds_sd_lastyr))
    min_sd <- min(minmax_sd[1, ])
    max_sd <- max(minmax_sd[2, ])
    
    col_scale_mn = scale_fill_gradientn(
      colors = hcl.colors(100, palette = "viridis"),
      limits = c(min_mn, max_mn),
      na.value = "transparent",
      guides(fill = ""))
    col_scale_sd = scale_fill_gradientn(
      colors = hcl.colors(100, palette = "viridis"),
      limits = c(min_sd, max_sd),
      na.value = "transparent",
      guides(fill = ""))
    
    plot_preds_mn_fy <- ggplot() + 
      geom_spatraster(data = preds_mn_firstyr, mapping = aes(fill = mean_firstyr)) + 
      col_scale_mn +
      theme_NPS + 
      ggtitle(paste0("Mean, ", min(YEARS))) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 8),
            legend.position = "bottom",
            legend.key.height = unit(0.2, "cm"),
            legend.key.width = unit(1, 'cm'),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_sd_fy <- ggplot() + 
      geom_spatraster(data = preds_sd_firstyr, mapping = aes(fill = sd_firstyr)) + 
      col_scale_sd +
      theme_NPS + 
      ggtitle(paste0("SD, ", min(YEARS))) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 8),
            legend.position = "bottom",
            legend.key.height = unit(0.2, "cm"),
            legend.key.width = unit(1, 'cm'),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_mn_ly <- ggplot() + 
      geom_spatraster(data = preds_mn_lastyr, mapping = aes(fill = mean_lastyr)) + 
      col_scale_mn +
      theme_NPS + 
      ggtitle(paste0("Mean, ", max(YEARS))) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 8),
            legend.position = "bottom",
            legend.key.height = unit(0.2, "cm"),
            legend.key.width = unit(1, 'cm'),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_sd_ly <- ggplot() + 
      geom_spatraster(data = preds_sd_lastyr, mapping = aes(fill = sd_lastyr)) + 
      col_scale_sd +
      theme_NPS + 
      ggtitle(paste0("SD, ", max(YEARS))) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 8),
            legend.position = "bottom",
            legend.key.height = unit(0.2, "cm"),
            legend.key.width = unit(1, 'cm'),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    title <- ggdraw() + 
      draw_label(paste0("Occurrence probability of ",
                        species$Common_name[species$Species_code == SPECIES],
                        " in ", park))
    p <- plot_grid(plot_preds_mn_fy, plot_preds_sd_fy,
                   plot_preds_mn_ly, plot_preds_sd_ly, nrow = 2)
    pp <- plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1.9))
    ggsave(pp,
           file = paste0(base_out, "map-occ-4panel", file_extension1),
           dpi = dpi,
           width = width_4,
           height = height_4,
           units = units)
    ggsave(pp,
           file = paste0(base_out, "map-occ-4panel", file_extension2),
           device = device,
           dpi = dpi,
           width = width_4,
           height = height_4,
           units = units)
  }
  if (MAP_DETECT) {

    if (!MAP) {
      stop("Cannot create map with detections figure with MAP set to FALSE")
    } 
    
    plot_preds_mn_naive_fy <- ggplot() + 
      geom_spatraster(data = preds_mn_firstyr, mapping = aes(fill = mean_firstyr)) + 
      scale_fill_viridis_c(na.value = 'transparent') +
      geom_spatvector(data=naive_spat[naive_spat$year==min(YEARS)], mapping=aes(color=as.factor(Present), shape=as.factor(Present))) +
      scale_color_detect() +
      scale_shape_detect() +
      annotation_scale(location = "br", style="ticks") +
      labs(fill = 'Probability', title = mn_title, subtitle = subtitle_fy, color="Present", shape="Present") +
      theme_NPS + 
      theme(axis.title = element_blank(),
            axis.line = element_blank())
    
    plot_preds_mn_naive_ly <- ggplot() + 
      geom_spatraster(data = preds_mn_lastyr, mapping = aes(fill = mean_lastyr)) + 
      scale_fill_viridis_c(na.value = 'transparent') +
      geom_spatvector(data=naive_spat[naive_spat$year==max(YEARS)], mapping=aes(color=as.factor(Present), shape=as.factor(Present))) +
      scale_color_detect() +
      scale_shape_detect() +
      annotation_scale(location = "br", style="ticks") +
      labs(fill = 'Probability', title = mn_title, subtitle = subtitle_ly, color="Present", shape="Present") +
      theme_NPS + 
      theme(axis.title = element_blank(),
            axis.line = element_blank())
    
    if (LATLONG) {
      plot_preds_mn_naive_fy <- plot_preds_mn_naive_fy +
        theme(panel.border = element_rect(color = 'black', fill = NA))
      plot_preds_mn_naive_ly <- plot_preds_mn_naive_ly +
        theme(panel.border = element_rect(color = 'black', fill = NA))
    } else {
      plot_preds_mn_naive_fy <- plot_preds_mn_naive_fy + 
        theme(axis.text = element_blank(),
              axis.ticks = element_blank())
      plot_preds_mn_naive_ly <- plot_preds_mn_naive_ly +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank())
          } 
    
      ggsave(plot_preds_mn_naive_fy, 
             file = paste0(base_out, "map-occ-detect-fy", file_extension1),
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      ggsave(plot_preds_mn_naive_fy, 
             file = paste0(base_out, "map-occ-detect-fy", file_extension2), 
             device = device,
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      ggsave(plot_preds_mn_naive_ly, 
             file = paste0(base_out, "map-occ-detect-ly", file_extension1),
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
      ggsave(plot_preds_mn_naive_ly, 
             file = paste0(base_out, "map-occ-detect-ly", file_extension2), 
             device = device,
             dpi = dpi, 
             width = width, 
             height = height, 
             units = units)
    }
}



