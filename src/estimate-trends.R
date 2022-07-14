################################################################################
# Estimate trends based on a multi-season dynamic model
# SAGW, LECA (black-tailed jackrabbit)

# ER Zylstra
# Updated 2022-07-14
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(jagsUI)
library(ggplot2)

# Load photo, location, events, species data 
source("src/format-mammal-data.R")

# dat = information about each photo (date, time, species, location)
# events = information about each camera deployment (dates, location, duration)
# event_mat = camera location x day matrix with 1/0 indicating whether camera
#             was deployed or not
# locs = information about each camera location (park, lat/long, name)
# species = table with species observed (species code, common name, # of obs)

# Load sampling occasion data (park, year, start/end, duration)
occasions <- read.csv("data/occasions/occasions-all-parks.csv")

park <- "SAGW"
species <- "LECA"

# Will eventually need some character string to identify the model with 
# particular covariates, information contained in the name of the model file.
# For now, will just use "MS-test" for model created in sagw-leca-multiseason.R

# Load JAGS model 
model_file <- paste0("output/models/",
                     tolower(park), "-",
                     tolower(species), "-",
                     "MS-test.rds")
jags_model <- readRDS(file = model_file)

# Extract posterior samples
samples <- jags_model$samples
samples <- do.call(rbind, samples)
pao <- samples[,grep("PAO", colnames(samples))]

# For each MCMC iteration, estimate a linear trend in occupancy
# (should probably do this on the logit scale)
# (note: this is a trend for surveyed locations only)
pao_logit <- log(pao / (1 - pao))
year_trend <- 0:(ncol(pao_logit) - 1)
trends <- data.frame(iter = 1:nrow(pao_logit), int = NA, slope = NA)
for (i in 1:nrow(pao_logit)) {
  m <- lm(pao_logit[i,] ~ year_trend)
  trends[i,2:3] <- coef(m)
}

# Summarize trend estimates
summary(trends$slope)

# Plot trends on the logit scale (each gray line represents one MCMC iteration)
trends <- trends %>%
  mutate(yr2022 = int + 5 * slope)

ggplot() +
  geom_segment(trends,
               mapping = aes(x = 2017, xend = 2022, y = int, yend = yr2022),
               size = 0.3, col = "gray") +
  geom_segment(trends,
               mapping = aes(x = 2017, xend = 2022, y = mean(int), yend = mean(yr2022)),
               size = 0.8, col = "dodgerblue3") +
  labs(x = "Year", y = "logit(Proportion of sites occupied)")


# TO DO:
# Plot (logit-linear) trends on the probability scale.



