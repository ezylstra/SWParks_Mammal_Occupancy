# SWParks_Mammal_Occupancy
Using camera trap data to assess occupancy of mammals in Southwestern National Parks

**Currently under development**

## Dependencies

+ downloader
+ dplyr
+ geojsonsf
+ ggplot2
+ gridExtra
+ httr
+ leaflet (mapping spatial data)
+ lubridate
+ raster
+ sf
+ spOccupancy
+ stringr
+ terra
+ tidyr
+ tidyterra
+ tigris (for downloading TIGER road files)
+ tmap (mapping spatial data)
+ unmarked (optional)

We're not currently using JAGS 4.3.1 to run single-season occupancy models, so 
we don't need to load the jagsUI, coda, or MCMCvis packages to run these 
analyses. Scripts that use JAGS (to run single-season or multi-season models) 
can be found in JAGS-labeled folders, but are not currently part of workflows. 

## General workflow

The first step is to organize and format raw photo observation, events, and 
location data. However, the scripts should only need to be modified once per 
year after all the photos are processed.  

**To run a single-season model (NEED TO ADD THIS)**

**To run a multi-season (dynamic) model** and estimate occupancy, detection 
probability, and dynamic parameters (colonization, extinction), you'll need to
run scripts in the src/multi-season-models folder. These files will call
scripts in other folders (e.g., src/photo-data, JAGS) as needed. 

1. Use [MSoccupancy-wrapper.R](src/multi-season-models/MSoccupancy-wrapper.R) to
   specify and run a multi-season occupancy model.
   + In the first part of the script (~ the first 100 lines), specify the 
     park, species, years, and covariates of interest (specifiable parameters in 
     all caps).  
   + Run the model in JAGS by calling 
     [MSoccupancy-generic.R](src/multi-season-models/MSoccupancy-generic.R) 
     using the source() function (line 103). 
   + The last few lines of the script will save the model output along with 
     other associated objects to file (in the output/models folder).

2. Use [MSoccupancy-results.R](src/multi-season-models/MSoccupancy-results.R) to
   process results from a multi-season model run previously. This script:
   + Loads a workspace with the jagsUI object and other dataframes
   + Creates a table that summarizes covariate effects (with names)
   + Creates a table with the estimated proportion of surveyed sites occupied 
     each year
   + Estimates (and plots) a logit-linear trend in occupancy for surveyed 
     locations
   + Provides example code to calculate and plot marginal covariate effects
   + (in development) Creates maps with the predicted probability of occupancy 
     in the first and last year

## Directory structure

+ data:
   + covariates: data (typically in rasters) that can be used as covariates in 
   occupancy models
   + locations: information about survey/camera locations
   + mammals: photo and sampling event data (_because these datasets potentially 
   include information about sensitive species, some of the raw data are not 
   publicly available_)
   + occasions: information about sampling occasions in each park and year 
   (created in  delineate-sampling-occasions.R)
+ src: 
  + photo-data: R scripts used to organize and format photo and location data.
  + covariate-data: R scripts used to obtain and format covariate data.
  + single-season-models: R scripts used to run process/format data, run 
  single-season occupancy models, and interpret/visualize results.
    + 2022: scripts used to run occupancy models for select species in various
    parks in 2022. 
    + JAGS: scripts used to run the models in JAGS (current not part of the 
    workflow).
  + functions.R: an R script that defines functions to be used in any script.    
  + map-detection-data: R scripts used to map detection data.
  + multi-season-models: R scripts used to run multi-season models in JAGS 
  (_need to work on this_)
+ output: Figures and tables summarizing survey effort and species detection
data (_some files may not be under version control, but directory 
structure is_).
   + models: R workspaces or objects with output from models run in JAGS
+ JAGS: JAGS model files (.txt files).

## Scripts

Scripts used to organize and format photo and location data -- found in the 
src/photo-data folder (should only need to be modified once per year)
   
1. [format-mammal-data.R](src/photo-data/format-mammal-data.R): Import, 
   organize, and clean up photo observation data, events data, and location data 
   (for all sites, years, species). 
   
2. [delineate-sampling-occasions.R](src/photo-data/delineate-sampling-occasions.R): 
   Create a table summarizing sampling occasions (start and end dates, duration) 
   at each park throughout the duration of the study. (this script calls 
   format-mammal-data.R)

3. [summarize-deployments-photos.R](src/photo-data/summarize-deployements-photos.R): 
   Summarize (and create plots to visualize) when and where cameras were 
   deployed. Create tables with the number of observations of each species, each 
   year (this script calls format-mammal-data.R)

Scripts used to obtain and format covariate data -- found in the 
src/covariate-data folder (should only need to be modified if new covariate
data are available or if running multi-season models)

1. [get-roads-trails-pois-data.R](src/covariate-data/get-roads-trails-pois-data.R): 
   Download feature data from NPS or other federal websites.

2. [prep-spatial-covariate-data.R](src/covariate-data/prep-covaraite-data.R): 
   Process spatial (time-invariant) data to generate rasters for each park.
   
3. [create-multi-layer-rasters.R](src/covariate-data/create-multi-layer-rasters.R):
   Create a multi-layer raster for each park that contains all the spatial
   (time-invariant) covariate data at the same resolution and extent.

4. [get-and-prep-climate-data.R](src/covariate-data/get-climate-data.R): 
   Download gridMET climate data and create rasters.

Scripts used to create detection histories and run single-season occupancy models
(**working on this**)
 
Scripts used to create detection histories and run multi-season occupancy models
(**need to review**)

1. [MSoccupancy-wrapper.R](src/multi-season-models/MSoccupancy-wrapper.R): Used
   to specify and run a multi-season occupancy model. Specifiable parameters
   (i.e., park, species, year, covariates) are denoted in all caps. There are
   options to create quadratic terms and/or interactions among covariates. The
   script runs the model in JAGS by calling 
   [MSoccupancy-generic.R](src/multi-season-models/MSoccupancy-generic.R) and 
   saves the jagsUI object with other relevant dataframes in an R workspace 
   saved to the output/models folders.

2. [MSoccupancy-generic.R](src/multi-season-models/MSoccupancy-generic.R): A 
   script called by MSoccupancy-wrapper.R that reads in photo and covariate 
   data, formats the data, and runs the model in JAGS.

3. [MSoccupancy-results.R](src/multi-season-models/MSoccupancy-results.R): Used
   to summarize, interpret, and visualize results from a multi-season occupancy
   model. 

4. [MSoccupancy-sagw-leca.R](src/multi-season-models/MSoccupancy-sagw-leca.R): 
   Example of a multi-season (dynamic) occupancy analysis, using data for 
   black-tailed jackrabbits in Saguaro National Park. _This file is now obsolete,
   since it is now more efficient to specify a model for jackrabbits in SAGW 
   using the MSoccupancy-wrapper.R script._

Scripts used to map detection data
(**need to review**)

1. [map-detection-data.R](src/map-detection-data/map-detection-data.R): Explore 
   how to create maps visualizing where species were detected in a park and 
   given year.

2. [map-detection-multiseason.R](src/map-detection-data/map-detection-multiseason.R): 
   Explore how to create maps visualizing where species were detected in a park 
   in all years.
