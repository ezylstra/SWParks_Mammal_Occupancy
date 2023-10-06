################################################################################
# Assess where more/fewer mammal species occur or are detected in a park

# ER Zylstra
# Updated 2023-07-19
################################################################################

library(tidyverse)
library(abind)
library(terra)
library(spOccupancy)
library(tidyterra)
library(RColorBrewer)

#------------------------------------------------------------------------------#
# Load detection data adn functions
#------------------------------------------------------------------------------#

source("src/photo-data/format-mammal-data.R")
source("src/functions.R")

#------------------------------------------------------------------------------#
# Specify parameters of interest
#------------------------------------------------------------------------------#

# Park, year
PARK <- "SAGW"
YEARS <- 2017:2022

# Logical indicating whether to include lat/longs on maps
LATLONG <- TRUE

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

# Create longer park name
park <- ifelse(PARK == "CHIR", "Chiricahua NM",
               ifelse(PARK == "SAGW", "Saguaro NP", "Organ Pipe Cactus NM"))

# Figure parameters
file_extension <- ".pdf"
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
  filter(Park == PARK) %>%
  select(loc, longitude, latitude) %>%
  rename(lon = longitude,
         lat = latitude)

# Load sampling occasion data and filter by park, years
occasions <- read.csv("data/occasions/occasions-all-parks.csv")
occasions <- occasions %>%
  filter(Park == PARK) %>%
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
  select(Common_name, Species_code) %>%
  mutate(modeled = 1 * Species_code %in% common_spp,
         rare = ifelse(modeled == 0 & Species_code != "OTVA", 1, 0))
# Labeling all unmodeled species as rare except for rock squirrels, that
# are relatively common but likely had few detections because of their size

# Filtering detection data (but note that we're not filtering by date so all
# detections of rare species are included)
dat_simple <- dat %>%
  filter(Park == PARK & yr %in% YEARS) %>%
  filter(Species_code %in% species$Species_code) %>%
  select(Species_code, obsdate, yr, loc)

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
  left_join(locs_simple, by = "loc") %>%
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

# Create figure with total number of species detected
mn_title <- "Number of species detected"
subtitle <- paste0(park, ", ", YEARS[1], "-", YEARS[length(YEARS)])
plot_nspp <- ggplot() + 
  geom_spatvector(data = contours, color = "gray65", fill = NA, linewidth = 0.2) +
  geom_spatvector(data = boundary, color = "black", fill = NA) +
  geom_spatvector(data = detsv, aes(color = factor(nspp)), size = 2) +
  scale_color_brewer(palette = "RdYlBu", name = "", direction = -1) +
  labs(fill = '', title = mn_title, subtitle = subtitle) +
  theme_NPS + 
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
                     "-nspp-detected.pdf"),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)

# Create figure with total number of rare species detected
mn_title <- "Number of rare species detected"
subtitle <- paste0(park, ", ", YEARS[1], "-", YEARS[length(YEARS)])
footnote <- paste(species$Common_name[species$rare == 1], collapse = ", ")
footnote <- paste0("Species included: ", footnote)
plot_nspp_rare <- ggplot() + 
  geom_spatvector(data = contours, color = "gray65", fill = NA, linewidth = 0.2) +
  geom_spatvector(data = boundary, color = "black", fill = NA) +
  geom_spatvector(data = detsv[detsv$nspp_rare > 0, ], 
                  aes(color = factor(nspp_rare)), size = 2) +
  scale_color_brewer(palette = "RdYlBu", name = "", direction = -1) +
  labs(title = mn_title, subtitle = subtitle, caption = str_wrap(footnote, 80)) +
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
                     "-nspp-rare-detected.pdf"),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)

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
    psi_spatcovs_z <- psi_covs_z[!psi_covs_z %in% c("years_z", "visits_z", "traffic_z")]
    psi_spatcovs <- psi_covs[!psi_covs %in% c("years", "visits", "traffic")]  
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
    nonspatial <- c("years", "visits", "traffic")
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
  scale_fill_viridis_c(na.value = 'transparent', name = "") +
  labs(title = mn_title, subtitle = subtitle, 
       caption = str_wrap(footnote, 80)) +
  theme_NPS + 
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
                     YEARS[length(YEARS)], "-spprichness-common.pdf"),
       device = device, 
       dpi = dpi, 
       width = width, 
       height = height, 
       units = units)
