# SWParks_Mammal_Occupancy
Using camera trap data to assess occupancy of mammals in Southwestern National Parks

**Currently under development**

## Dependencies

JAGS 4.3.1 (At least at the moment, we're using JAGS to run Bayesian occupancy models. This version of JAGS is needed if running R 4.2.x)

The project also uses the following R packages:

+ coda (optional)
+ dplyr
+ ggplot2
+ gridExtra
+ jagsUI
+ leaflet (mapping spatial data)
+ lubridate
+ MCMCvis (optional)
+ sf
+ stringr
+ tidyr
+ tmap (mapping spatial data)
+ unmarked (optional)

## General workflow

1. Organize and format raw photo observation, events, and location data
2. Create detection histories for selected park(s), species(s), and year(s)
3. Run single-season (static) or multi-season (dynamic) models to estimate occupancy, detection probability or other parameters (e.g., colonization, extinction, trends in occupancy)

## Scripts

1. Organize and format raw data
   1. **format-mammal-data.R**: Import, organize, and clean up photo observation data, events data, and location data (for all sites, years, species). 
   2. **delineate-sampling-occasions.R**: Create a table summarizing sampling occasions (start and end dates, duration) at each park throughout the duration of the study. (this script calls format-mammal-data.R)
   3. **data-exploration.R**: Summarize (and create plots to visualize) when and where cameras were deployed (this script calls format-mammal-data.R)
2. Create detection histories and run occupancy model(s)
   1. **sagw-leca-2017.R**: Example of a single-season occupancy analysis, using data for black-tailed jackrabbits at Saguaro National Park in 2017.
   2. **JAGS_SingleSeasonWithCovs.txt**: Example of a JAGS model file. Used to run the single-season model for jackrabbits in Saguaro.
   3. **map-detection-data.R**: Create map with jackrabbit detection data (using NPS imagery)
   4. **sagw-leca-multiseason.R**: Example of a multi-season (dynamic) occupancy anaysis, using data for blck-tailed jackrabbits at Saguaro National Park.
   5. **JAGS_MultiSeasonWithCovs.txt**: Example of a JAGS model file. Used to run the multi-season, dynamic model for jackrabbits in Saguaro.

## Directory structure

+ data:
   + locations: information about survey/camera locations
   + mammals: photo and sampling event data (because these datasets potentially include information about sensitive species, some of the raw data are not publicly available)
   + occasions: information about sampling occasions in each park and year (created in delineate-sampling-occasions.R)
   + rasters: rasters with covariate data
+ JAGS: JAGS model files (.txt files)
+ output: (note: some files may not be under version control, but directory structure is)
   + models: output from models run in JAGS (jagsUI objects)
+ src: R scripts to process/format data and to run occupancy models in JAGS and interpret/visualize the results.

