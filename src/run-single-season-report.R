################################################################################
# Create a single-season report for a specified park and year

  # For this to work, models for one or more species in the specified park 
  # and year need to be saved here: output/single-single-models/
  # Files should be named: PARK-YEAR-SPECIES.rds
  
  # The resulting word doc will be saved here: output/single-season-reports/
  # File will be named: PARK-YEAR.docx
################################################################################

# Specify the park and year of interest

  park <- "CHIR"
  year <- 2022

# Install any needed packages if they're not installed already
  
  packages <- c("lubridate", "tidyverse", "terra", "spOccupancy", "tidyterra",
                "officer", "officedown", "flextable", "rmarkdown", "knitr")
  
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }

# Create the report

  suppressWarnings(
    rmarkdown::render("src/single-season-report.Rmd",
                      params = list(PARK = park,
                                    YEAR = year,
                                    date = Sys.Date()),
                      output_file = paste0("../output/single-season-reports/", 
                                           park, "-", year))
  )

  # Note that for render(), the default directory for the output file is where 
  # the input file is located (here, src/)
