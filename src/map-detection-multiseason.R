################################################################################
# Create maps to visualize camera trap data

# C McIntyre
# 2022-07-22
################################################################################



# relies on the dh_long and locs_park from the sagw-leca-multiseason code
# learning how to make maps in R and branches in GitHub

library(sf) # for working with simple features
library(tmap) # for mapping spatial data
library(leaflet) # Another package for mapping spatial data
library(dplyr)

# re-import observation data to get scientific and common name of species back
# Observations
species_name <- read.csv("data/mammals/MAMMALS_ALL_2022-04-18.csv") %>%
  dplyr::select(Species, Common.name, Species.code) %>%
  distinct() %>%
  rename(Scientific_name = Species, Common_name = Common.name, Species = Species.code)


# Compute presence/absence (naive) and percent of weeks/occ present within each year
naive <- dh_df %>%
  pivot_longer(!loc,
               names_to = c("year","week"),
               names_pattern = "(.*)_(.)",
               values_to = "det") %>%
  mutate(Species = species) %>%
  group_by(loc,year, Species) %>%
  filter(!is.na(det)) %>% #remove NAs from calculations rather than setting them to 0
  summarize(Pct_Present = sum(det)/n(), Present = ifelse(Pct_Present>0,1,0), .groups="keep") %>%
  mutate(Detection = ifelse(Present==1,"Detected","Not Detected")) %>%
  arrange(loc, year, Species) %>%
  left_join(species_name, by = "Species")

# Proportion of cameras with detection - naive estimate, not accounting for 
naive %>% ungroup() %>% group_by(year, Species) %>% summarize(Pct_Detected = sum(Present)/n(), .groups="keep")


# make into simple feature
naive_sf <- st_as_sf(naive %>% 
                       #tibble::rownames_to_column(., "loc") %>%
                       left_join(., locs_park, by = "loc"), 
                     coords=c("long","lat"), crs=4326) # not 100% that WGS84 is correct

naive_present_sf <- naive_sf %>% filter(Present == 1)
naive_absent_sf <- naive_sf %>% filter(Present == 0 | is.na(Present))



# extract species scientific name and common name
sci_name <- as.character(species_name %>% 
                           filter(Species == species) %>% 
                           dplyr::select(Scientific_name))
common_name <- as.character(species_name %>% 
                           filter(Species == species) %>% 
                           dplyr::select(Common_name))

# code below based on NPS IMD Intro to R Training: https://katemmiller.github.io/IMD_R_Training_Intro/

# load NPS boundaries (split into sub-units)
  # NPS boundaries from 20220106
park_boundaries1 <- st_read('./data/locations/SWNCwildlife_nps_boundary_20220106.shp')

# check that projections match
st_crs(naive_sf) == st_crs(park_boundaries1) # FALSE; park_boundaries crs = 4269, naive_sf crs = 4326

# project park_boundaries to 4326 (WGS84)
park_boundaries <- st_transform(park_boundaries1, crs=4326)

# quick plot of park boundaries, by unit_code [column 2]
plot(park_boundaries[2])

# filter park_boundaries to select park unit of interest
unit_boundary <- park_boundaries %>% filter(UNIT_CODE == park)

# quick plot of unit boundary, by unit_code [column 2]
plot(unit_boundary[2])

# extract unit name 
unit_name <- unit_boundary$UNIT_NAME




# Map collage with presence/absence by year for a given species
  # not interactive but potentially informative
# To Do
  # want to define my own colors for detected & not detected


# make bounding box for map to see entire unit boundary
yearly_map_bb <- st_bbox(unit_boundary) 

yearly_map <-
  # unit boundary
  tm_shape(unit_boundary) +
  tm_borders('black',lwd=1) +
  
  # detections, bboc = sets the bounding box equal to the unit_boundary extent
  tm_shape(naive_sf, bbox=yearly_map_bb) +
  tm_dots("Detection", size=1, palette=c("Detected"='blue', "Not Detected" = 'grey')) + 
  
  # facet by year, free.coords = FALSE forces each facet to use the bounding box set above
  tm_facets(by="year", free.coords = FALSE) + 
  
  # Other map features
  #tm_compass(size = 2, type = 'arrow', text.size = 1, position = c('left', 'bottom')) +
  #tm_scale_bar(text.size = 1.25, position = c('center', 'bottom')) + 
  tm_layout(inner.margins = c(0.2, 0.02, 0.02, 0.02), # make room for legend
            outer.margins = 0,
            legend.just = 'right',
            legend.position = c('left', 'top'),
            main.title = paste("Naive estimates of", common_name, "at", unit_name, sep=" "), main.title.size = 1)

yearly_map



