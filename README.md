# SWParks_Mammal_Occupancy
Using camera trap data to assess occupancy of mammals in Southwestern National Parks

**Currently under development**

## Dependencies

+ downloader
+ dplyr
+ flextable (for creating reports)
+ geojsonsf
+ ggplot2
+ gridExtra
+ httr
+ knitr
+ leaflet (mapping spatial data)
+ lubridate
+ ncdf4
+ officedown (for creating reports)
+ officer (for creating reports)
+ raster
+ rmarkdown
+ sf
+ spOccupancy
+ stringr (v 1.5.0 or higher)
+ terra
+ tidyr
+ tidyterra
+ tidyverse
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

2. Save this script as: src/single-season-models/YEAR/spOccupancy-PARK-YEAR-SPECIES.R.

3. Continue running this script from line 68 through the section labeled "Select
   best model of candidate set". In this part of the script, you'll specify 
   covariates for a set of candidate models, run the models using the 
   spOccupancy package, and create a table comparing model fit. Then you'll
   select the best model of the set and remove any extraneous covariates in the
   detection or occurrence part of the model. At the end of this section you'll 
   save the model object to file (in the output/single-season-models folder, the
   contents of which are **not** under version control). 

4. Run the rest of the script, beginning with the section "Evaluate best model
   and look at estimates". Here, you'll assess convergence and goodness-of-fit, 
   produce tables with parameter estimates, generate a figure with predicted 
   occurrence probabilities, and calculate marginal covariate effects.

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
   src/multi-season-models/PARK/spOccupancy-PARK-FIRSTYEAR-LASTYEAR-SPECIES.R.
  
3. Continue running this script from line 63 through the section labeled "Select
   best model". In this part of the script, you'll specify covariates for 
   various sets of candidate models, run the models using the spOccupancy 
   package, and create tables comparing model fit. Ultimately, you'll select the
   "best"" model and remove any extraneous covariates in the detection or 
   occurrence part of the model. At the end of this section you'll save the 
   model object to file (in the output/multi-season-models folder, the
   contents of which are **not** under version control). 

4. Run the rest of the script, beginning with the section "Evaluate best model
   and look at estimates". Here, you'll assess convergence and goodness-of-fit, 
   produce tables with parameter estimates, generate figures with predicted 
   occurrence probabilities in the first and last year, and calculate marginal 
   covariate effects.

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
  + multi-season-models: R scripts used to process/format data, run 
  multi-season occupancy models, and interpret/visualize results.
    + CHIR: scripts used to run multi-season occupancy models for select species
    in Chiricahua NM.
    + ORPI: scripts used to run multi-season occupancy models for select species 
    in Organ Pipe NM.
    + SAGW: scripts used to run multi-season occupancy models for select species
    in Saguaro NP.
    + JAGS: scripts used to run multi-season models in JAGS (no longer part of
    the workflow.)
  + functions.R: an R script that defines custom functions to be used in any 
  script in the repo.    
  + run-single-season-report.R: an R script used to create a word document (via 
  markdown) that summarizes detection data and occurrence models across all 
  species for a given PARK and YEAR. (This script calls various .Rmd files in 
  the src folder)
  + create-single-season-figures.R: an R script used to create and save NPS 
  formatted figures for public-facing versions of single-season reports.
  + map-detection-data: R scripts used to map raw detection data.
+ output: Figures and tables summarizing survey effort and species detection
data (_some files may not be under version control, but directory 
structure is_).
  + single-season-models: folder that contains output from single-season models
  run with the spOccupancy package (eg, output from the "best" model for
  coyotes in SAGW in 2022 would be saved as SAGW-2022-CALA.rds).
  + multi-season-models: folder that contains output from multi-season models
  run with the spOccupancy package (eg, output from the "best" model for
  coyotes in SAGW, 2017-2022 would be saved as SAGW-2017-2022-CALA.rds).
  + single-season-reports: folder that contains word documents that summarize 
  detection data and occurrence models across all species for a given PARK and 
  YEAR (eg, SAGW-2022.docx)
  + NPS-figures
    + single-season: pdf versions of NPS-formatted figures that summarize 
    results from single-season models (occurrence probability maps or covariate 
    marginal effects)
    + multi-season: pdf versions of NPS-formatted figures that summarize 
    results from multi-season analyses (occurrence probability maps, covariate 
    marginal effects, maps with species detections or estimated species 
    richness)
  + models-JAGS: R workspaces or objects with output from models run in JAGS
+ JAGS: JAGS model files (.txt files).

## Scripts

**Scripts used to organize and format photo and location data** -- found in the 
src/photo-data folder (should only need to be modified once per year)
   
1. [format-mammal-data.R](src/photo-data/format-mammal-data.R): Import, 
   organize, and clean up photo observation data, events data, and location data 
   (for all sites, years, species). 
   
2. [delineate-sampling-occasions.R](src/photo-data/delineate-sampling-occasions.R): 
   Create a table summarizing sampling occasions (start and end dates, duration) 
   at each park throughout the duration of the study. (this script calls 
   format-mammal-data.R)

3. [summarize-deployments-photos.R](src/photo-data/summarize-deployments-photos.R): 
   Summarize (and create plots to visualize) when and where cameras were 
   deployed. Create tables with the number of observations of each species, each 
   year. (this script calls format-mammal-data.R)

**Scripts used to obtain and format covariate data** -- found in the 
src/covariate-data folder (should only need to be modified if new covariate
data are available or in the case of climate data, modified each year if running 
multi-season models)

1. [create-multi-layer-rasters.R](src/covariate-data/create-multi-layer-rasters.R):
   Create a multi-layer raster for each park that contains all the spatial
   (time-invariant) covariate data at the same resolution and extent.

2. [get-and-prep-climate-data.R](src/covariate-data/get-and-prep-climate-data.R): 
   Download gridMET climate data and create rasters.

3. [get-roads-trails-pois-data.R](src/covariate-data/get-roads-trails-pois-data.R): 
   Download feature data from NPS or other federal websites.

4. [prep-spatial-covariate-data.R](src/covariate-data/prep-spatial-covariate-data.R): 
   Process spatial (time-invariant) data to generate rasters for each park.
   
**Scripts used to create detection histories and run single-season occupancy models**

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
 
**Scripts used to create detection histories and run multi-season occupancy models**

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
   
5. [spOccupancy-MS-predictions.R](src/multi-season-models/spOccupancy-MS-predictions.R):
   Calculate predicted occurrence probabilities across a park in the first and
   last year, and create ggplot objects that depict mean and SDs in each raster 
   cell. (This script will typically be called from 
   src/multi-season-models/YPARK/spOccupancy-PARK-SPECIES-YEARS.R.)

6. [map-species-detections.R](src/multi-season-models/map-species-detections.R):
   Create maps that depict the number of mammal species detected and the 
   estimated number of species that occur in locations throughout the park.

**Scripts used to create single-season reports and NPS quality figures**

1. [run-single-season-reports.R](src/run-single-season-reports.R): Use markdown
   to create a word document that summarizes detection data and results from 
   occurrence models across all species for a given PARK and YEAR. (This script 
   calls [single-season-report.Rmd](src/single-season-report.Rmd),
   [single-season-report-spp.Rmd](src/single-season-report-spp.Rmd), and
   [single-season-report-marg-figs.Rmd](src/single-season-report-marg-figs.Rmd))

2. [create-single-season-figures.R](src/create-single-season-figures.R): Create 
   and save NPS formatted figures for public-facing versions of single-season
   reports.

**Scripts used to map raw detection data**

1. [map-detection-data.R](src/map-detection-data/map-detection-data.R): Explore 
   how to create maps visualizing where species were detected in a park and 
   given year.

2. [map-detection-multiseason.R](src/map-detection-data/map-detection-multiseason.R): 
   Explore how to create maps visualizing where species were detected in a park 
   in all years.
