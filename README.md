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
+ stringr (v 1.5.0 or higher)
+ terra
+ tidyr
+ tidyterra
+ tigris (for downloading TIGER road files)
+ tmap (mapping spatial data)
+ unmarked (optional)

We are not currently using JAGS 4.3.1 to run single-season occupancy models, so 
we don't need to load the jagsUI, coda, or MCMCvis packages to run these 
analyses. Scripts that use JAGS (to run single-season or multi-season models) 
can be found in JAGS-labeled folders, but are not part of current workflows. 

## General workflow

The first step is to organize and format raw photo observation, events, and 
location data. However, the scripts should only need to be modified once per 
year after all the photos are processed.  

**To run a single-season model** and estimate occurrence and detection
probabilities, you'll need to access files in the src/single-season-models
folder. These files will call scripts in other folders (e.g., src/functions.R, 
src/photo-data/format-mammal-data.R) as needed. _Note: we are now using the 
spOccupancy package to run single-season models._

1. Open [spOccupancy-TEMPLATE.R](src/single-season-models/spOccupancy-TEMPLATE.R) 
   and run through line 68, where you'll specify the park (PARK), year (YEAR), 
   and species (SPECIES) of interest. 

2. Save this script as: src/single-season-models/YEAR/spOccupancy-PARK-SPECIES-YEAR.R.

3. Continue running this script from line 68 through the section labeled "Run 
   models and compare fit". In this part of the script, you'll specify 
   covariates for a set of candidate models, run the models using the 
   spOccupancy package, and create a table comparing model fit.
   
4. Run the rest of the script, beginning with the section "Look at results and 
   predictions from best model". Here, you'll select a model for inference,
   assess convergence and goodness-of-fit, produce tables with parameter 
   estimates, generate a figure with predicted occurrence probabilities, and 
   calculate marginal covariate effects.

**To run a multi-season (dynamic) model** and estimate annual occurrence and 
detection probabilities, you'll need to access files in the 
src/multi-season-models folder. These files will call scripts in other folders 
(e.g., src/photo-dataformat-mammal-data.R) as needed. These are not explicit 
dynamic models (that estimate colonization and extinction parameters), but
instead model annual occurrence probabilities (that can vary as a function of 
spatial or non-spatial annual covariates). _Note: we are now using the 
spOccupancy package to run multi-season models._

1. Open [spOccupancy-MS-TEMPLATE.R](src/multi-season-models/spOccupancy-MS-TEMPLATE.R) 
   and run through line 62, where you'll specify the park (PARK), years (YEARS), 
   and species (SPECIES) of interest. 

2. Save this script as: 
   src/multi-season-models/PARK/spOccupancy-PARK-SPECIES-FIRSTYEAR-LASTYEAR.R.
  
3. Continue running this script from line 63 through the section labeled "Run 
   models and compare fit". In this part of the script, you'll specify 
   covariates for a set of candidate models, run the models using the 
   spOccupancy package, and create a table comparing model fit.

4. Run the rest of the script, beginning with the section "Look at results and 
   predictions from best model". Here, you'll select a model for inference,
   assess convergence and goodness-of-fit, produce tables with parameter 
   estimates, generate figures with predicted occurrence probabilities in
   the first and last year, and calculate marginal covariate effects.

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
  + single-season-models: R scripts used to process/format data, run 
  single-season occupancy models, and interpret/visualize results.
    + 2022: scripts used to run occupancy models for select species in various
    parks in 2022. 
    + JAGS: scripts used to run the models in JAGS (no longer part of the 
    workflow).
  + functions.R: an R script that defines custom functions to be used in any 
  script in the repo.    
  + map-detection-data: R scripts used to map detection data.
  + multi-season-models: R scripts used to process/format data, run 
  multi-season occupancy models, and interpret/visualize results.
    + CHIR: scripts used to run multi-season occupancy models for select species
    in Chiricahua NM.
    + JAGS: scripts used to run multi-season models in JAGS (no longer part of
    the workflow.)
    + ORPI: scripts used to run multi-season occupancy models for select species 
    in Organ Pipe NM.
    + SAGW: scripts used to run multi-season occupancy models for select species
    in Saguaro NP.
+ output: Figures and tables summarizing survey effort and species detection
data (_some files may not be under version control, but directory 
structure is_).
   + models-JAGS: R workspaces or objects with output from models run in JAGS
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
   year. (this script calls format-mammal-data.R)

Scripts used to obtain and format covariate data -- found in the 
src/covariate-data folder (should only need to be modified if new covariate
data are available or in the case of climate data, modified each year if running 
multi-season models)

1. [create-multi-layer-rasters.R](src/covariate-data/create-multi-layer-rasters.R):
   Create a multi-layer raster for each park that contains all the spatial
   (time-invariant) covariate data at the same resolution and extent.

2. [get-and-prep-climate-data.R](src/covariate-data/get-climate-data.R): 
   Download gridMET climate data and create rasters.

3. [get-roads-trails-pois-data.R](src/covariate-data/get-roads-trails-pois-data.R): 
   Download feature data from NPS or other federal websites.

4. [prep-spatial-covariate-data.R](src/covariate-data/prep-covaraite-data.R): 
   Process spatial (time-invariant) data to generate rasters for each park.
   
Scripts used to create detection histories and run single-season occupancy models

1. [spOccupancy-TEMPLATE.R](src/single-season-models/spOccupancy-TEMPLATE.R): 
   A template to create a script that will run single-season occupancy models 
   for a selected park, year, and species. (this script calls src/functions.R 
   and several scripts in the src/single-season-models folder [see below])
   
2. [spOccupancy-data-prep.R](src/single-season-models/spOccupancy-data-prep.R):
   Formats detection and covariate data and packages everything into a list 
   needed to run single-season models in the spOccupancy package. (This script
   will typically be called from 
   src/single-season-models/YEAR/spOccupancy-PARK-SPECIES-YEAR.R.)
   
3. [spOccupancy-create-model-formulas.R](src/single-season-models/spOccupancy-create-model-formulas.R):
   Develop formulas for occurrence and detection parts of candidate models. (This 
   script will typically be called from 
   src/single-season-models/YEAR/spOccupancy-PARK-SPECIES-YEAR.R.)
   
4. [spOccupancy-run-candidate-models.R](src/single-season-models/spOccupancy-run-candidate-models.R):
   Run candidate single-season models (with spOccupancy) and gather model 
   diagnostics and statistics for comparisons. (This script will typically be 
   called from src/single-season-models/YEAR/spOccupancy-PARK-SPECIES-YEAR.R.)
   
5. [spOccupancy-predictions.R](src/single-season-models/spOccupancy-predictions.R):
   Calculate predicted occurrence probabilities across a park, and create 
   ggplot objects that depict mean and SDs in each raster cell. (This script
   will typically be called from src/single-season-models/YEAR/spOccupancy-PARK-SPECIES-YEAR.R.)
 
Scripts used to create detection histories and run multi-season occupancy models

1. [spOccupancy-MS-TEMPLATE.R](src/multi-season-models/spOccupancy-MS-TEMPLATE.R): 
   A template to create a script that will run multi-season occupancy models 
   for a selected park, species, and set of years. (this script calls src/functions.R 
   and several scripts in the src/multi-season-models folder [see below])
   
2. [spOccupancy-MS-data-prep.R](src/multi-season-models/spOccupancy-MS-data-prep.R):
   Formats detection and covariate data and packages everything into a list 
   needed to run multi-season models in the spOccupancy package. (This script
   will typically be called from 
   src/multi-season-models/PARK/spOccupancy-PARK-SPECIES-YEARS.R.)
   
3. [spOccupancy-MS-create-model-formulas.R](src/multi-season-models/spOccupancy-MS-create-model-formulas.R):
   Develop formulas for occurrence and detection parts of candidate models. (This 
   script will typically be called from 
   src/multi-season-models/PARK/spOccupancy-PARK-SPECIES-YEARS.R.)
   
4. [spOccupancy-MS-run-candidate-models.R](src/multi-season-models/spOccupancy-MS-run-candidate-models.R):
   Run candidate multi-season models (with spOccupancy) and gather model 
   diagnostics and statistics for comparisons. (This script will typically be 
   called from src/multi-season-models/PARK/spOccupancy-PARK-SPECIES-YEARS.R.)
   
5. [spOccupancy-MS-predictions.R](src/multiseason-models/spOccupancy-MS-predictions.R):
   Calculate predicted occurrence probabilities across a park in the first and
   last year, and create ggplot objects that depict mean and SDs in each raster 
   cell. (This script will typically be called from 
   src/multi-season-models/YPARK/spOccupancy-PARK-SPECIES-YEARS.R.)
   
Scripts used to map detection data

1. [map-detection-data.R](src/map-detection-data/map-detection-data.R): Explore 
   how to create maps visualizing where species were detected in a park and 
   given year.

2. [map-detection-multiseason.R](src/map-detection-data/map-detection-multiseason.R): 
   Explore how to create maps visualizing where species were detected in a park 
   in all years.
