################################################################################
# Create maps to visualize camera trap data

# C McIntyre
# 2022-07-22
################################################################################



# relies on the dh_long and locs_park from the sagw-leca-multiseason code
# learning how to make maps in R and branches im GitHub

library(sf) # for working with simple features
library(tmap) # for mapping spatial data
library(leaflet) # Another package for mapping spatial data
library(dplyr)

# Compute presence/absence (naive occupancy) and percent of weeks/occ present within each year
naive <- dh_df %>%
  pivot_longer(!loc,
               names_to = c("year","week"),
               names_pattern = "(.*)_(.)",
               values_to = "det") %>%
  mutate(Species = species) %>%
  group_by(loc,year, Species) %>%
  filter(!is.na(det)) %>% #remove NAs from calculations rather than setting them to 0
  summarize(Pct_Present = sum(det)/n(), Present = ifelse(Pct_Present>0,1,0), .groups="keep")
  arrange(loc, year, Species) 

# make into 
naive_sf <- st_as_sf(naive %>% 
                       #tibble::rownames_to_column(., "loc") %>%
                       left_join(., locs_park, by = "loc"), 
                     coords=c("long","lat"), crs=4326)

naive_present_sf <- naive_sf %>% filter(Present == 1)
naive_absent_sf <- naive_sf %>% filter(Present == 0 | is.na(Present))

# code below based on NPS IMD Intro to R Training: https://katemmiller.github.io/IMD_R_Training_Intro/

# Map collage with presence/absence by year for a given species
  # probably won't be interactive but could be informative




# Try to make a map

# Load park tiles

NPSbasic = 'https://atlas-stg.geoplatform.gov/styles/v1/atlas-user/ck58pyquo009v01p99xebegr9/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1IjoiYXRsYXMtdXNlciIsImEiOiJjazFmdGx2bjQwMDAwMG5wZmYwbmJwbmE2In0.lWXK2UexpXuyVitesLdwUg'

NPSimagery = 'https://atlas-stg.geoplatform.gov/styles/v1/atlas-user/ck72fwp2642dv07o7tbqinvz4/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1IjoiYXRsYXMtdXNlciIsImEiOiJjazFmdGx2bjQwMDAwMG5wZmYwbmJwbmE2In0.lWXK2UexpXuyVitesLdwUg'

NPSslate = 'https://atlas-stg.geoplatform.gov/styles/v1/atlas-user/ck5cpvc2e0avf01p9zaw4co8o/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1IjoiYXRsYXMtdXNlciIsImEiOiJjazFmdGx2bjQwMDAwMG5wZmYwbmJwbmE2In0.lWXK2UexpXuyVitesLdwUg'

NPSlight = 'https://atlas-stg.geoplatform.gov/styles/v1/atlas-user/ck5cpia2u0auf01p9vbugvcpv/tiles/256/{z}/{x}/{y}@2x?access_token=pk.eyJ1IjoiYXRsYXMtdXNlciIsImEiOiJjazFmdGx2bjQwMDAwMG5wZmYwbmJwbmE2In0.lWXK2UexpXuyVitesLdwUg'

naive_map_lf <-
  leaflet() %>% 
  setView(lng = -111.166, lat = 32.298, zoom = 13) %>% 
  # parktiles
  addTiles(group = 'Map',
           urlTemplate = NPSbasic) %>%
  addTiles(group = 'Imagery',
           urlTemplate = NPSimagery) %>%
  addTiles(group = 'Light',
           urlTemplate = NPSlight) %>%
  addTiles(group = 'Slate',
           urlTemplate = NPSslate) %>% 
  addLayersControl(map = .,
                   baseGroups = c('Map', 'Imagery', 'Light', 'Slate'),
                   options = layersControlOptions(collapsed = T)) %>% 
  # points
  # addCircleMarkers(data = naive_present_sf[year=2017],
  #                  lng = st_coordinates(naive_present_sf)[,1],
  #                  lat = st_coordinates(naive_present_sf)[,2],
  #                  fillColor = '#56B4E9',
  #                  radius = 4,
  #                  stroke = FALSE, # turn off outline
  #                  fillOpacity = 1) %>% 
  addCircleMarkers(data = naive_sf_absent[year=2017],
                   lng = st_coordinates(naive_sf_absent)[,1],
                   lat = st_coordinates(naive_sf_absent)[,2],
                   fillColor = '#99999',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>%
  # scale bar and settings
  addScaleBar(position = 'bottomright') %>% 
  scaleBarOptions(maxWidth = 10, metric = TRUE) 

naive_map_lf

# make map by percent weeks with detection 
naive_sf_0 <- naive_sf %>% filter(Pct_Present ==0)
naive_sf_20 <- naive_sf %>% filter(Pct_Present ==0.2)
naive_sf_40 <- naive_sf %>% filter(Pct_Present ==0.4)
naive_sf_60 <- naive_sf %>% filter(Pct_Present ==0.6)
naive_sf_80 <- naive_sf %>% filter(Pct_Present ==0.8)
naive_sf_100 <- naive_sf %>% filter(Pct_Present ==1)

naive_map_pct <-
  leaflet() %>% 
  setView(lng = -111.166, lat = 32.298, zoom = 12) %>% 
  # parktiles
  addTiles(group = 'Map',
           urlTemplate = NPSbasic) %>%
  addTiles(group = 'Imagery',
           urlTemplate = NPSimagery) %>%
  addTiles(group = 'Light',
           urlTemplate = NPSlight) %>%
  addTiles(group = 'Slate',
           urlTemplate = NPSslate) %>% 
  addLayersControl(map = .,
                   baseGroups = c('Map', 'Imagery', 'Light', 'Slate'),
                   options = layersControlOptions(collapsed = T)) %>% 
  # points
  addCircleMarkers(data = naive_sf_0,
                   lng = st_coordinates(naive_sf_0)[,1],
                   lat = st_coordinates(naive_sf_0)[,2],
                   fillColor = '#FFFFFF',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>% 
  addCircleMarkers(data = naive_sf_20,
                   lng = st_coordinates(naive_sf_20)[,1],
                   lat = st_coordinates(naive_sf_20)[,2],
                   fillColor = '#d5d2db',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>%
  addCircleMarkers(data = naive_sf_40,
                   lng = st_coordinates(naive_sf_40)[,1],
                   lat = st_coordinates(naive_sf_40)[,2],
                   fillColor = '#aba5b7',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>%
  addCircleMarkers(data = naive_sf_60,
                   lng = st_coordinates(naive_sf_60)[,1],
                   lat = st_coordinates(naive_sf_60)[,2],
                   fillColor = '#827893',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>%
  addCircleMarkers(data = naive_sf_80,
                   lng = st_coordinates(naive_sf_80)[,1],
                   lat = st_coordinates(naive_sf_80)[,2],
                   fillColor = '#584b6f',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>%
  addCircleMarkers(data = naive_sf_100,
                   lng = st_coordinates(naive_sf_100)[,1],
                   lat = st_coordinates(naive_sf_100)[,2],
                   fillColor = '#2e1e4b',
                   radius = 4,
                   stroke = FALSE, # turn off outline
                   fillOpacity = 1) %>%
  # scale bar and settings
  addScaleBar(position = 'bottomright') %>% 
  scaleBarOptions(maxWidth = 10, metric = TRUE) 

naive_map_pct 
