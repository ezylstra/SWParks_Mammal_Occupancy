################################################################################
# Estimate trends based on a multi-season dynamic model
# SAGW, LECA (black-tailed jackrabbit)

# ER Zylstra
# Updated 2022-07-21
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(jagsUI)
library(ggplot2)

rm(list = ls())

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

# Extract 1000 (of 3000) posterior samples for PAO estimates
samples <- jags_model$samples
samples <- do.call(rbind, samples)
samples <- samples[seq(1, 3000, by = 3),]
pao <- samples[,grep("PAO", colnames(samples))]

# For each MCMC iteration, estimate a linear trend in logit(occupancy)
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
hist(trends$slope, breaks = 50)

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

# Prep data to plot (logit-linear) trends on the probability scale.
  # Create a sequence of values that spans yr_trend
  yr_predict <- seq(0, max(year_trend), length = 100)
  # Create a matrix, with a column of 1s (for the intercept) and yr_predict
  yr_predict_matrix <- rbind(1, yr_predict)
  # Predict values: logit(occupancy) = beta0 + slope*yr_predict (matrix math)
  preds_logit <- as.matrix(trends[,c("int", "slope")]) %*% yr_predict_matrix
  # Convert predictions to probability scale. 
  preds_prob <- exp(preds_logit)/(1 + exp(preds_logit))
  colnames(preds_prob) <- yr_predict
  # Calculate the median value across iterations for each value of yr_predict
  preds_median <- data.frame(Year = yr_predict + 2017, 
                             Occupancy = apply(preds_prob, 2, median))
  # Add column to identify MCMC iteration
  preds_prob <- cbind(preds_prob, iter = 1:nrow(preds_prob))
  # Convert predictions to long form for ggplot
  preds_long <- preds_prob %>%
    as.data.frame %>%
    pivot_longer(cols = !iter,
                 names_to = "Year",
                 values_to = "Occupancy") %>%
    mutate(Year = as.numeric(Year) + 2017) %>%
    as.data.frame

# Plot trends on the probability scale (each line represents one MCMC iteration)
ggplot() +
  geom_line(preds_long,
            mapping = aes(x = Year, y = Occupancy, group = iter),
            col = "gray") + 
  geom_line(preds_median,
            mapping = aes(x = Year, y = Occupancy),
            size = 0.8, col = "dodgerblue3") +
  labs(y = "Proportion of sites occupied")
  