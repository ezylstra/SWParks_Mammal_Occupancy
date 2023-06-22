# Create single-season report for a specified park and year

park <- "SAGW"
year <- 2022

suppressWarnings(
  rmarkdown::render("src/single-season-report.Rmd",
                    params = list(PARK = park,
                                  YEAR = year,
                                  date = Sys.Date()),
                    output_file = paste0("../output/single-season-reports/", 
                                         park, "-", year))
)

# Note that the default directory for the output file is where the input 
# file is located (here, src/)
