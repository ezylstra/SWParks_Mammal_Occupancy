################################################################################
# Assess where more/fewer mammal species occur or are detected in a park

# ER Zylstra
# Updated 2023-10-13
################################################################################

library(tidyverse)
library(abind)
library(terra)
library(spOccupancy)
library(tidyterra)
library(RColorBrewer)
library(ggspatial)
library(raster)
library(exactextractr)

#------------------------------------------------------------------------------#
# Load detection data and functions
#------------------------------------------------------------------------------#

# Select park of interest ("CHIR", "ORPI", or "SAGW")
PARK <- "SAGW"

source("src/photo-data/format-mammal-data.R")
source("src/functions.R")

#------------------------------------------------------------------------------#
# Specify parameters of interest
#------------------------------------------------------------------------------#

# Select years
YEARS <- 2017:2025

# Logical indicating whether to include lat/longs on maps
LATLONG <- FALSE

# Create custom NPS theme
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

# Create custom color ramp from Paul Tol's sunset color ramp
# (this way 4 is the same color on all maps)
scale_color_count <- function(...){
  ggplot2:::manual_scale('colour', 
                         values = setNames(c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", 
                                             "#FDB366","#F67E4B", "#DD3D2D", "#A50026"),
                                           c( "0","1","2", "3", "4","5", "6", "7", "8", "9", "10")), 
                         ...)
}

# Create custom symbols
scale_shape_count <- function(...){
  ggplot2:::manual_scale('shape', 
                         values = setNames(c(1, 2, 0, 5, 18, 4, 3, 15, 17, 16, 8, 11),
                                           c("0", "1","2", "3", "4","5", "6", "7", "8", "9", "10")), 
                         ...)
}


# Create longer park name
park <- ifelse(PARK == "CHIR", "Chiricahua NM",
               ifelse(PARK == "SAGW", "Saguaro NP", "Organ Pipe Cactus NM"))

# Figure parameters
file_extension1 <- ".pdf"
file_extension2 <- ".png"
device <- cairo_pdf
dpi <- 300
width <- 6
height <- 4
units <- "in"



#------------------------------------------------------------------------------#
# Summarize detection data
#------------------------------------------------------------------------------#

# Extract camera locations for this park
locs_simple <- locs %>%
  dplyr::select(loc, longitude, latitude) %>%
  rename(lon = longitude,
         lat = latitude)

# Load sampling occasion data and filter by park, years
occasions <- read.csv(paste0("data/occasions/occasions-", PARK, ".csv"))
occasions <- occasions %>%
  filter(yr %in% YEARS)

# Get list of common species (that we ran occupancy models for)
spp_rds <- list.files(path = "output/multi-season-models/", 
                      pattern = paste0(PARK, "-", YEARS[1], "-", 
                                       YEARS[length(YEARS)]),
                      full.names = TRUE)
common_spp <- basename(spp_rds) %>% str_sub(16, 19)

# Filtering out non-natives or unknowns from species list
species <- species %>%
  filter(Nativeness == "Native" & !is.na(Nativeness)) %>%
  dplyr::select(Common_name, Species, Species_code) %>%
  mutate(modeled = 1 * Species_code %in% common_spp,
         rare = ifelse(modeled == 0 & Species_code != "OTVA", 1, 0))
# Labeling all unmodeled species as rare except for rock squirrels, that
# are relatively common but likely had few detections because of their size

# Update common name for PETA (javelina instead of collared peccary)
# and URCI (gray fox instead of common gray fox)
species <- species %>%
  mutate(Common_name = ifelse(Species_code=="PETA","javelina", Common_name)) %>%
  mutate(Common_name = ifelse(Species_code=="URCI","gray fox", Common_name))

# Filtering detection data (but note that we're not filtering by date so all
# detections of rare species are included)
dat_simple <- dat %>%
  filter(Park == PARK & yr %in% YEARS) %>%
  filter(Species_code %in% species$Species_code) %>%
  dplyr::select(Species_code, obsdate, yr, loc)

# Create list of species detected in the park
spp_detect <- dat_simple %>%
  group_by(Species_code) %>%
  summarize(ndetects = length(Species_code),
            nyrs = length(unique(yr)),
            nlocs = length(unique(loc))) %>%
  data.frame()
species <- left_join(species, spp_detect, by = "Species_code") %>%
  filter(!is.na(ndetects))

# Calculate the number of species detected at each camera location
dets <- dat_simple %>%
  left_join(species[, c("Species_code", "modeled", "rare")], 
            by = "Species_code") %>%
  group_by(loc) %>%
  distinct(Species_code, modeled, rare) %>%
  summarize(nspp = length(Species_code),
            nspp_modeled = sum(modeled),
            nspp_rare = sum(rare)) %>%
  right_join(locs_simple, by = "loc") %>%
  mutate(nspp = ifelse(is.na(nspp), 0, nspp)) %>%
  mutate(nspp_modeled = ifelse(is.na(nspp_modeled), 0, nspp_modeled)) %>%
  mutate(nspp_rare = ifelse(is.na(nspp_rare), 0, nspp_rare)) %>%
  mutate(plot = 0) %>% # to get camera location to appear on species maps even if no observations
  data.frame()

# Load park boundary
boundary <- vect("data/covariates/shapefiles/Boundaries_3parks.shp")
boundary <- subset(boundary, boundary$UNIT_CODE == PARK)
# Load DEM
dem <- rast(paste0("data/covariates/DEMs/", PARK, "_DEM_1as.tif"))
contours <- as.contour(dem, maxcells = Inf)
contours <- crop(contours, boundary)
# Make dets into a SpatVector
detsv <- vect(dets, geom = c("lon", "lat"), crs = crs(boundary))

# Buffer park boundary (for road clipping)
park_boundary_1km <- buffer(boundary, width=1000, singlesided=FALSE)

# Load and clip trails layer to park boundary
park_trails <- vect("data/covariates/shapefiles/trails.shp")
# clip to current park
park_trails <- crop(park_trails, boundary)

# Load roads shapefile (within 3km) and clip to within 1km
park_roads_file <- ifelse(PARK=="SAGW", "data/covariates/shapefiles/roads_sagw_v2.shp", ifelse(PARK=="CHIR", "data/covariates/shapefiles/roads_chir_nps_usfs.shp", "data/covariates/shapefiles/roads_orpi_nps.shp"))
park_roads <- vect(park_roads_file)
#park_roads <- if(PARK=="SAGW") vect("data/covariates/shapefiles/roads_sagw_v2.shp") else vect(paste0("data/covariates/shapefiles/roads_",PARK,"_tigris.shp", sep=""))
park_roads_1km <- crop(park_roads, park_boundary_1km)

# Create figure with total number of species observed
mn_title <- "Number of species observed"
subtitle <- paste0(park, ", ", YEARS[1], "-", YEARS[length(YEARS)])
plot_nspp <- ggplot() + 
  #geom_spatvector(data = contours, color = "gray65", fill = NA, linewidth = 0.2) +
  geom_spatvector(data = boundary, color = "darkgreen", fill = "grey", lwd=1) +
  geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
  geom_spatvector(data=park_roads_1km, color="black", inherit.aes=FALSE, lwd = 0.1) + 
  geom_spatvector(data = detsv, aes(color = factor(nspp), shape=factor(nspp)), size = 2) +
  scale_color_count(name = "Species") + 
  scale_shape_count(name = "Species") +
  #scale_color_brewer(palette = "RdYlBu", name = "Species", direction = -1) +
  labs(fill = '', title = mn_title, subtitle = subtitle) +
  theme_NPS + 
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
  annotation_scale(location = "br", style="ticks") +
  theme(axis.title = element_blank(),
        axis.line = element_blank())
if (LATLONG) {
  plot_nspp <- plot_nspp +
    theme(panel.border = element_rect(color = 'black', fill = NA))
} else {
  plot_nspp <- plot_nspp + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}  

ggsave(plot_nspp, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", 
                     YEARS[1], "-", YEARS[length(YEARS)], 
                     "-nspp-detected", file_extension1),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)
ggsave(plot_nspp, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", 
                     YEARS[1], "-", YEARS[length(YEARS)], 
                     "-nspp-detected", file_extension2),
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)

# Create figure with total number of rare (uncommon) species observed
mn_title <- "Number of uncommon species observed"
subtitle <- paste0(park, ", ", YEARS[1], "-", YEARS[length(YEARS)])
footnote <- paste(species$Common_name[species$rare == 1], collapse = ", ")
footnote <- paste0("Species included: ", footnote)
plot_nspp_rare <- ggplot() + 
  #geom_spatvector(data = contours, color = "gray65", fill = NA, linewidth = 0.2) +
  geom_spatvector(data = boundary, color = "darkgreen", fill = "grey", lwd=1) +
  geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
  geom_spatvector(data=park_roads_1km, color="black", inherit.aes=FALSE, lwd = 0.1) + 
  #geom_spatvector(data = detsv[detsv$nspp_rare == 0, ], size = 0.5, color="white") +
  #geom_spatvector(data = detsv[detsv$nspp_rare > 0, ], 
  #                aes(color = factor(nspp_rare)), size = 1) +
  geom_spatvector(data = detsv, aes(color = factor(nspp_rare), shape=factor(nspp_rare)), size = 2) +
  scale_color_count(name = "Species") + 
  scale_shape_count(name = "Species") +
  #scale_color_brewer(palette = "RdYlBu", name = "Species", direction = -1) +
  labs(title = mn_title, subtitle = subtitle, caption = str_wrap(footnote, 80)) +
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
  annotation_scale(location = "br", style="ticks") +
  theme_NPS + 
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7))
if (LATLONG) {
  plot_nspp_rare <- plot_nspp_rare +
    theme(panel.border = element_rect(color = 'black', fill = NA))
} else {
  plot_nspp_rare <- plot_nspp_rare + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}  

ggsave(plot_nspp_rare, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", 
                     YEARS[1], "-", YEARS[length(YEARS)], 
                     "-nspp-uncommmon-detected", file_extension1),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)
ggsave(plot_nspp_rare, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", 
                     YEARS[1], "-", YEARS[length(YEARS)], 
                     "-nspp-rare-detected", file_extension2),
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)

#------------------------------------------------------------------------------#
# Create detection maps for individual common (modeled) species 
#------------------------------------------------------------------------------#
modeled_species <- species %>% filter(modeled==1) %>% dplyr::select(Species_code) %>% pull
# Calculate the number of species detected at each camera location
obs_modeled <- dat_simple %>%
  left_join(species[, c("Species_code", "modeled", "rare")], 
            by = "Species_code") %>%
  group_by(Species_code) %>%
  distinct(Species_code, yr, loc, modeled, rare) %>%
  group_by(Species_code, loc, modeled, rare) %>%
  summarize(yrs = length(yr), .groups = "keep") %>%
  pivot_wider(., names_from = "loc", values_from = "yrs", values_fill = 0) %>%
  pivot_longer(., cols=-c(Species_code, modeled, rare), names_to = "loc", values_to = "yrs") %>%
  filter(modeled==1) %>%
  left_join(locs_simple, by = "loc") %>%
  data.frame()
# Make dets into a SpatVector
obs_modeled_sv <- vect(obs_modeled, geom = c("lon", "lat"), crs = crs(boundary))

# Create figure with total years of species observations
for (i in unique(modeled_species)){
  spp_common <- species$Common_name[species$Species_code==i]
  mn_title <- paste("Number of years with", spp_common, "observations")
  subtitle <- paste0(park, ", ", YEARS[1], "-", YEARS[length(YEARS)])
  plot_nyrs <- ggplot() + 
    #geom_spatvector(data = contours, color = "gray65", fill = NA, linewidth = 0.2) +
    geom_spatvector(data = boundary, color = "darkgreen", fill = "grey", lwd=1) +
    geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
    geom_spatvector(data=park_roads_1km, color="black", inherit.aes=FALSE, lwd = 0.1) + 
    geom_spatvector(data = obs_modeled_sv[obs_modeled_sv$Species_code==i], 
                    aes(color = factor(yrs), shape = factor(yrs)), size = 2) +
    #scale_color_brewer(palette = "RdYlBu", name = "Years", direction = -1) +
    scale_color_count(name = "Years") + 
    scale_shape_count(name = "Years") +
    labs(fill = '', title = mn_title, subtitle = subtitle) +
    theme_NPS + 
    annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
    annotation_scale(location = "br", style="ticks") +
    theme(axis.title = element_blank(),
          axis.line = element_blank())
  if (LATLONG) {
    plot_nyrs <- plot_nyrs +
      theme(panel.border = element_rect(color = 'black', fill = NA))
  } else {
    plot_nyrs <- plot_nyrs + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
    print(plot_nyrs)
  }
  
  ggsave(plot_nyrs, 
         file = paste0("output/NPS-figures/multi-season/", PARK, "-", spp_common, "-",
                       YEARS[1], "-", YEARS[length(YEARS)], 
                       "-nyr-obs", file_extension1),
         device = device, 
         dpi = dpi, 
         width = width, 
         height = height, 
         units = units)
  
  ggsave(plot_nyrs, 
         file = paste0("output/NPS-figures/multi-season/", PARK, "-", spp_common, "-",
                       YEARS[1], "-", YEARS[length(YEARS)], 
                       "-nyr-obs", file_extension2),
         dpi = dpi, 
         width = width, 
         height = height, 
         units = units)
}

#------------------------------------------------------------------------------#
# Create detection maps for individual rare species 
#------------------------------------------------------------------------------#
rare_species <- species %>% filter(rare==1) %>% dplyr::select(Species_code) %>% pull
# Calculate the number of species detected at each camera location
obs_rare <- dat_simple %>%
  left_join(species[, c("Species_code", "modeled", "rare")], 
            by = "Species_code") %>%
  group_by(Species_code) %>%
  distinct(Species_code, yr, loc, modeled, rare) %>%
  group_by(Species_code, loc, modeled, rare) %>%
  summarize(yrs = length(yr), .groups = "keep") %>%
  pivot_wider(., names_from = "loc", values_from = "yrs", values_fill = 0) %>%
  pivot_longer(., cols=-c(Species_code, modeled, rare), names_to = "loc", values_to = "yrs") %>%
  filter(rare==1) %>%
  left_join(locs_simple, by = "loc") %>%
  data.frame()
# Make dets into a SpatVector
obs_rare_sv <- vect(obs_rare, geom = c("lon", "lat"), crs = crs(boundary))

# Create figure with total years of species observations
for (i in unique(rare_species)){
spp_common <- species$Common_name[species$Species_code==i]
mn_title <- paste("Number of years with", spp_common, "observations")
subtitle <- paste0(park, ", ", YEARS[1], "-", YEARS[length(YEARS)])
plot_nyrs <- ggplot() + 
  #geom_spatvector(data = contours, color = "gray65", fill = NA, linewidth = 0.2) +
  geom_spatvector(data = boundary, color = "darkgreen", fill = "grey", lwd=1) +
  geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
  geom_spatvector(data=park_roads_1km, color="black", inherit.aes=FALSE, lwd = 0.1) + 
  geom_spatvector(data = obs_rare_sv[obs_rare_sv$Species_code==i], 
                  aes(color = factor(yrs), shape = factor(yrs)), size = 2) +
  #scale_color_brewer(palette = "RdYlBu", name = "Years", direction = -1) +
  scale_color_count(name = "Years") + 
  scale_shape_count(name = "Years") +
  labs(fill = '', title = mn_title, subtitle = subtitle) +
  theme_NPS + 
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
  annotation_scale(location = "br", style="ticks") +
  theme(axis.title = element_blank(),
        axis.line = element_blank())
if (LATLONG) {
  plot_nyrs <- plot_nyrs +
    theme(panel.border = element_rect(color = 'black', fill = NA))
} else {
  plot_nyrs <- plot_nyrs + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  print(plot_nyrs)
}

ggsave(plot_nyrs, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", spp_common, "-",
                     YEARS[1], "-", YEARS[length(YEARS)], 
                     "-nyr-obs", file_extension1),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)

ggsave(plot_nyrs, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", spp_common, "-",
                     YEARS[1], "-", YEARS[length(YEARS)], 
                     "-nyr-obs", file_extension2),
        dpi = dpi, 
       width = width, 
       height = height, 
       units = units)
}
#------------------------------------------------------------------------------#
# Load occurrence probabilities for common species and summarize
#------------------------------------------------------------------------------#
# Note: creating rasters with predicted occurrence probabilities takes several 
# minutes

for (i in 1:length(spp_rds)) {
  
  SPECIES <- common_spp[i]
  
  model_list <- readRDS(spp_rds[i])
  best <- model_list$model
  best_psi_model <- model_list$psi_model
  best_p_model <- model_list$p_model
  data_list <- model_list$data
  rm(model_list)
  
  if (i == 1) {
    # Use one species to get spatial_covs dataframe, park_raster, etc. 
    source("src/multi-season-models/spOccupancy-MS-data-prep.R")
  }

  # Extract names of covariates (with and without "_z" subscripts) from best model
  # And for occurrence, extract names of spatial covariates
  psi_covs_z <- create_cov_list(best_psi_model)
  if (length(psi_covs_z) == 1 & any(psi_covs_z == "1")) {
    psi_covs_z <- character(0)
    psi_covs <- character(0)
    psi_spatcovs_z <- character(0)
    psi_spatcovs <- character(0)
  } else {
    psi_covs <- psi_covs_z %>% str_remove_all(pattern = "_z")
    psi_spatcovs_z <- psi_covs_z[!psi_covs_z %in% c("years_z", "visits_z", "traffic_z", "monsoon_ppt_z", "ppt10_z", "ppt6_z",
                                                    "monsoon_vpd_z","vpd10_z", "vpd6_z", "aet10_z","deficit10_z", "savi_z")]
    psi_spatcovs <- psi_covs[!psi_covs %in% c("years", "visits", "traffic", "monsoon", "ppt10", "ppt6",
                                              "monsoon_vpd","vpd10", "vpd6", "aet10","deficit10", "savi")]  
  }
  p_covs_z <- create_cov_list(best_p_model)
  if (length(p_covs_z) == 1 & any(p_covs_z == "1")) {
    p_covs_z <- character(0)
    p_covs <- character(0)
  } else {
    p_covs <- p_covs_z %>% str_remove_all(pattern = "_z")
  }

  occ_estimates <- parameter_estimates(model = best, 
                                       parameter = "occ",
                                       lower_ci = 0.025,
                                       upper_ci = 0.975)
  occ_estimates <- occ_estimates %>%
    rename(Covariate = Parameter) %>%
    mutate(Parameter = "Occurrence", .before = "Covariate")
  
  # Note: if there are time-varying covariates (other than year/trend) in the 
  # occurrence part of the model, we'll be estimating occurrence probabilities
  # in the last year under observed conditions (e.g., observed 10-month
  # precipitation)
  if (length(psi_spatcovs) > 0) {
    ANN_PREDS <- "observed"
    source("src/multi-season-models/spOccupancy-MS-predictions.R")
    # Create occrast_SPECIES raster with predicted occurrence probabilities in 
    # last year
    assign(paste0("occrast_", SPECIES), preds_mn_lastyr)
  } else {
    nonspatial <- c("years", "visits", "traffic", "monsoon", "ppt10", "ppt6",
                    "monsoon_vpd","vpd10", "vpd6", "aet10","deficit10", "savi")
    beta_samples <- as.matrix(best$beta.samples)
    # If there are only non-spatial covariates in the model, calculate mean 
    # predicted occurrence probability in the last year.
    if (any(nonspatial %in% psi_covs)) {
      nonspatial_names <- colnames(beta_samples[, -1])
      nonspatial_values <- NA
      for (j in 1:length(nonspatial_names)) {
        nonspatial_values[j] <- get(nonspatial_names[j])[1, length(YEARS)]
      }
      nonspatial_values <- as.matrix(c(1, nonspatial_values))
      pred <- t(nonspatial_values) %*% t(beta_samples)
      pred <- exp(pred) / (1 + exp(pred))
      mean_occ <- mean(pred)
    } else {
      # If there are no covariates in the model, calculate the overall mean
      mean_occ <- mean(exp(beta_samples[,1]) / (1 + exp(beta_samples[,1])))
    }
    assign(paste0("occmean_", SPECIES), mean_occ)
  }  
}

# Combine species rasters (and mean occurrence probabilities for species without
# spatial covariates in the model for occurrence) and sum probabilities to 
# estimate how many of the modeled species are likely present at that location 
# in the last year
occrast_common <- rast(mget(str_subset(ls(), "occrast_")))
occrast_common <- sum(occrast_common)
if (length(str_subset(ls(), "occmean_")) > 0) {
  mean_list <- mget(str_subset(ls(), "occmean_"))
  sum_nonspatial <- sum(unlist(mean_list))
  occrast_common <- occrast_common + sum_nonspatial
}

# Create figure with estimated number of common species
mn_title <- "Estimated number of commmon mammal species"
subtitle <- paste0(park, ", ", YEARS[length(YEARS)])
footnote <- paste(species$Common_name[species$modeled == 1], collapse = ", ")
footnote <- paste0("Species included: ", footnote)
plot_spprich <- ggplot() + 
  geom_spatraster(data = occrast_common, mapping = aes(fill = sum)) + 
  scale_fill_viridis_c(na.value = 'transparent', name = "Species") +
  geom_spatvector(data=park_trails, color="lightgrey", lwd = 0.25, linetype = "longdash") +
  geom_spatvector(data=park_trails, color="black", lwd = 0.1, linetype = "dashed") +
  geom_spatvector(data=park_roads_1km, color="lightgrey", inherit.aes=FALSE, lwd = 0.5) + 
  geom_spatvector(data=park_roads_1km, color="black", inherit.aes=FALSE, lwd = 0.1) + 
  labs(title = mn_title, subtitle = subtitle, 
       caption = str_wrap(footnote, 120)) +
  theme_NPS + 
  annotation_north_arrow(location = "bl", which_north = "true", style = north_arrow_minimal()) +
  annotation_scale(location = "br", style="ticks") +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7))
if (LATLONG) {
  plot_spprich <- plot_spprich +
    theme(panel.border = element_rect(color = 'black', fill = NA))
} else {
  plot_spprich <- plot_spprich + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
}  

ggsave(plot_spprich, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", 
                     YEARS[length(YEARS)], "-spprichness-common",file_extension1),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)

ggsave(plot_spprich, 
       file = paste0("output/NPS-figures/multi-season/", PARK, "-", 
                     YEARS[length(YEARS)], "-spprichness-common",file_extension2),
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)
