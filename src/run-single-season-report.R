# Render test

park <- "SAGW"
year <- 2022

rmarkdown::render("src/single-season-report.Rmd",
                  params = list(PARK = park,
                                YEAR = year,
                                date = Sys.Date()),
                  output_dir = paste0(getwd(), "/output/single-season-reports/"), 
                  output_file = paste0(park, "-", year))