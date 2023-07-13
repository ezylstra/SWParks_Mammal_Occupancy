################################################################################
# Script to make/save NPS formatted figures for single-season reports

# ER Zylstra
# Updated 2023-07-13
################################################################################

library(tidyverse)
library(terra)
library(spOccupancy)
library(tidyterra)

# Specify park, year, species of interest
PARK <- "SAGW"
YEAR <- 2022
SPECIES <- "URCI"

# Logical indicating whether to create a map with mean occurrence probabilities 
MAP <- TRUE
# Logical indicating whether to create a map with SD of occurrence probabilities
MAP_SD <- FALSE
  # If creating maps, indicate whether to include lat/long axes labels
  LATLONG <- FALSE

# Logical indicating whether to create figures with marginal covariate effects
# for the occurrence part of the model
MARG_OCC <- TRUE

# Logical indicating whether to create figures with marginal covariate effects
# for the detection part of the model
MARG_DET <- TRUE

# Load model object 
modelfile <- paste0("output/single-season-models/", PARK, "-", YEAR, "-",
                    SPECIES, ".rds")
if (!file.exists(modelfile)) {
  stop("Single-season model for ", SPECIES, " in ", PARK, " in ", YEAR, 
          " has not been saved in ouput/single-season-models folder.")
}
model_list <- readRDS(modelfile)
best <- model_list$model
psi_model <- model_list$psi_model
p_model <- model_list$p_model
data_list <- model_list$data

# Create basename for output files
base_out <- paste0("output/NPS-figures/single-season/",
                   PARK, "-", YEAR, "-", SPECIES, "-")

# Load functions and raw data
source("src/functions.R")
source("src/photo-data/format-mammal-data.R")

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

# Create longer park name for use in plots
park <- ifelse(PARK == "CHIR", "Chiricahua NM",
               ifelse(PARK == "SAGW", "Saguaro NP", "Organ Pipe Cactus NM"))

#------------------------------------------------------------------------------#
# Maps with occurrence probabilities (mean, SD)
#------------------------------------------------------------------------------#
if (MAP) {
  
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
  
  # Extract dataframe with covariate values at each camera location
  spatial_covs <- data_list$occ.covs
  
  # Extract names of covariates from best model
  psi_covs_z <- create_cov_list(psi_model)
  if (length(psi_covs_z) == 1 & any(psi_covs_z == "1")) {
    psi_covs_z <- character(0)
    psi_covs <- character(0)
  } else {
    psi_covs <- psi_covs_z %>% str_remove_all(pattern = "_z")
  }
  
  # Generate predicted probabilities (preds_mn raster)
  source("src/single-season-models/spOccupancy-predictions.R")
  
  # Create figures
  mn_title <- paste0("Mean occurrence probability of ",
                     species$Common_name[species$Species_code == SPECIES])
  sd_title <- paste0("Standard deviation of occurrence probability of ",
                     species$Common_name[species$Species_code == SPECIES])
  subtitle <- paste0(park, ", ", YEAR)
  
  plot_preds_mn <- ggplot() + 
    geom_spatraster(data = preds_mn, mapping = aes(fill = mean)) + 
    scale_fill_viridis_c(na.value = 'transparent') +
    # geom_spatvector(data = park_boundary, fill = NA, color = "black", size = 0.5) + 
    labs(fill = '', title = mn_title, subtitle = subtitle) +
    theme_NPS + 
    theme(axis.title = element_blank(),
          axis.line = element_blank())
  plot_preds_sd <- ggplot() + 
    geom_spatraster(data = preds_sd, mapping = aes(fill = sd)) + 
    scale_fill_viridis_c(na.value = 'transparent') +
    # geom_spatvector(data = park_boundary, fill = NA, color = "black", size = 0.5) + 
    labs(fill = '', title = sd_title, subtitle = subtitle) +
    theme_NPS +
    theme(axis.title = element_blank(),
          axis.line = element_blank())
  if (LATLONG) {
    plot_preds_mn <- plot_preds_mn +
      theme(panel.border = element_rect(color = 'black', fill = NA))
    plot_preds_sd <- plot_preds_sd +
      theme(panel.border = element_rect(color = 'black', fill = NA))
  } else {
    plot_preds_mn <- plot_preds_mn + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
    plot_preds_sd <- plot_preds_sd + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
  }  

  ggsave(plot_preds_mn, file = paste0(base_out, "occmap-mn.pdf"),
         device = cairo_pdf, dpi = 300, width = 6, height = 4, units = "in")
  if (MAP_SD) {
    ggsave(plot_preds_sd, file = paste0(base_out, "occmap-sd.pdf"),
           device = cairo_pdf, dpi = 300, width = 6, height = 4, units = "in")
  }
}  
  
#------------------------------------------------------------------------------#
# Marginal covariate effects for the occurrence part of the model
# (Creates figures for all covariates)
#------------------------------------------------------------------------------#  
  
  # Load covariates.csv?

#------------------------------------------------------------------------------#
# Marginal covariate effects for the detection part of the model
# (Creates figures for all covariates)
#------------------------------------------------------------------------------#  

# Load covariates.csv?
