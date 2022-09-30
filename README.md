# SWParks_Mammal_Occupancy
Using camera trap data to assess occupancy of mammals in Southwestern National Parks

**Currently under development**

## Dependencies

JAGS 4.3.1 (At least at the moment, we're using JAGS to run Bayesian occupancy 
models. This version of JAGS is needed if running R 4.2.x)

The project also uses the following R packages:

+ coda (optional)
+ downloader
+ dplyr
+ ggplot2
+ gridExtra
+ jagsUI
+ leaflet (mapping spatial data)
+ lubridate
+ MCMCvis (optional)
+ sf
+ stringr
+ terra
+ tidyr
+ tigris (for downloading TIGER road files)
+ tmap (mapping spatial data)
+ unmarked (optional)

## General workflow

The first step is to organize and format raw photo observation, events, and 
location data. However, the scripts should only need to be modified once per 
year after all the photos are processed.  

**To run a multi-season (dynamic) model** to estimate occupancy, detection 
probability or other parameters (e.g., colonization, extinction), you will only
need to run scripts in the src/multi-season-models folder. These files will call
scripts in other folders (e.g., src/photo-data, JAGS) as needed. Eventually, 
there will be a wrapper where you select the species, park, and covariates you
would like to include in a multi-season model.  For now:

1. Run [MSoccupancy-sagw-leca.R](src/multi-season-models/MSoccupancy-sagw-leca.R) 
   or a similarly named script to run an analysis for a particular park 
   (e.g., SAGW) and species (e.g., LECA = black-tailed jackrabbit). In this 
   script, you can specify which covariates should be included in the model. 
   Model output will be saved in the output/models folder.

2. Run [estimate-trends.R](src/multi-season-models/estimate-trends.R) to 
   estimate logit-linear trends in the estimated proportion of sites occupied 
   over time, based on the dynamic, multi-season model (this will load the saved 
   model object from the output/models folder).

Similarly, to run a single-season model, you will only need to run a script in 
the src/single-season-models folder.

1. Run [SSoccupancy-sagw-leca-2017.R](src/single-season-models/SSoccupancy-sagw-leca-2017.R) 
   or a similarly named script to run an analysis for a particular park 
   (e.g., SAGW), species (e.g., LECA = black-tailed jackrabbit), and year (e.g., 
   2017). In this script, you can specify which covariates should be included in 
   the model.

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
+ JAGS: JAGS model files (.txt files)
+ output: (note: some files may not be under version control, but directory 
structure is)
   + models: output from models run in JAGS (jagsUI objects)
+ src: R scripts to process/format data and to run occupancy models in JAGS and 
interpret/visualize the results.

## Scripts

Scripts used to create detection histories and run occupancy models
   
1. [MSoccupancy-sagw-leca.R](src/multi-season-models/MSoccupancy-sagw-leca.R): 
Example of a multi-season (dynamic) occupancy analysis, using data for 
black-tailed jackrabbits in Saguaro National Park.
   
2. [JAGS_MultiSeasonWithCovs.txt](JAGS/JAGS_MultiSeasonWithCovs.txt): Example of 
a JAGS model file. Used to run the multi-season, dynamic model for jackrabbits 
in Saguaro. Don't need to run this file independently - it will be called from 
MSoccupancy scripts.

3. [estimate-trends.R](src/multi-season-models/estimate-trends.R): Script to 
estimate logit-linear trends in the estimated proportion of sites occupied over 
time, based on a dynamic, multi-season model.

4. [SSoccupancy-sagw-leca-2017.R](src/single-season-models/SSoccupancy-sagw-leca-2017.R): 
Example of a single-season occupancy analysis, using data for black-tailed 
jackrabbits at Saguaro National Park in 2017.
   
5. [JAGS_SingleSeasonWithCovs.txt](JAGS/JAGS_SingleSeasonWithCovs.txt): Example 
of a JAGS model file. Used to run the single-season model for jackrabbits in 
Saguaro.  

Scripts used to organize and format photo and location data (should only need to 
be modified once per year)
   
1. [format-mammal-data.R](src/photo-data/format-mammal-data.R): Import, 
organize, and clean up photo observation data, events data, and location data 
(for all sites, years, species). 
   
2. [delineate-sampling-occasions.R](src/photo-data/delineate-sampling-occasions.R): 
Create a table summarizing sampling occasions (start and end dates, duration) at 
each park throughout the duration of the study. (this script calls 
format-mammal-data.R)

3. [summarize-deployments-photos.R](src/photo-data/summarize-deployements-photos.R): 
Summarize (and create plots to visualize) when and where cameras were deployed. 
Create tables with the number of observations of each species, each year (this 
script calls format-mammal-data.R)

Scripts used to obtain and format covariate data 

1. [get-climate-data.R](src/covariate-data/get-climate-data.R): Download gridMET 
climate data

2. [get-roads-trails-pois-data.R](src/covariate-data/get-roads-trails-pois-data.R): 
Download feature data from NPS or other federal websites.

3. [prep-covariate-data.R](src/covariate-data/prep-covaraite-data.R): Process 
data to generate rasters from which we can extract covariate values for camera 
or other locations.

Scripts used to map detection data

1. [map-detection-data.R](src/map-detection-data/map-detection-data.R): Explore 
how to create maps visualizing where species were detected in a park and given 
year.

2. [map-detection-multiseason.R](src/map-detection-data/map-detection-multiseason.R): 
Explore how to create maps visualizing where species were detected in a park in 
all years.
   
