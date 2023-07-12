################################################################################
# Script to make/save NPS formatted figures for single-season reports

# ER Zylstra
# Updated 2023-07-12
################################################################################

library(ggplot2)
library(tidyterra)

# Specify park, year, species of interest and load model object 
PARK <- "SAGW"
YEAR <- 2022
SPECIES <- "CALA"

modelfile <- paste0("output/single-season-models/", PARK, "-", YEAR, "-",
                    SPECIES, ".rds")
model <- readRDS(modelfile)

# Create custom theme
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

# For occurrence probability maps
  # Need to load rasters too

  # map <- ggplot()...
  # fig_NPS <- map + theme_NPS
  fig_name <- "occmap"
  fig_file <- paste0("output/NPS-figures/", PARK, "-", YEAR, "-", SPECIES, "-",
                     fig_name, "-4NPS.pdf")
  ggsave(fig_NPS, 
         file = fig_file, 
         device = cairo_pdf, 
         dpi = 300, 
         width = 6, 
         height = 4, 
         units="in")

# For marginal plots
  # Need covariates.csv

