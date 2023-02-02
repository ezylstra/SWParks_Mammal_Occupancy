################################################################################
# Functions to be used in any script in SWParks_Mammal_Occupancy repo

# ER Zylstra
# Updated 2023-02-02
################################################################################

#------------------------------------------------------------------------------#
# Convert a model formula to a vector with covariate names
#------------------------------------------------------------------------------#

# object: a character string representing a formula for occurrence or 
  # detection probability in an occupancy model (should begin with a ~)
# z_list: a character vector that contains all covariates in a model formula 
  # (or returns a "1" if there are no covariates)

create_cov_list <- function(object) {
  if (!methods::is(object, "character")) {
    stop("Object should be a character.")
  }

  z_list <- object %>%
    str_remove(pattern = "~ ") %>% 
    str_remove_all(pattern = "I[(]") %>%
    str_remove_all(pattern = "[)]") %>%
    str_remove_all(pattern = "\\^2") %>%
    str_split_1(pattern = " [+] ")

    return(z_list)
}
  
  
  