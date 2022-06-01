################################################################################
# SODN -- Camera trap data, 2016-2022
# Data exploration

# ER Zylstra
# 2022-05-19
################################################################################

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# rm(list = ls())

#-------------------------------------------------------------------------------#
# Run script to import, format data
#-------------------------------------------------------------------------------#

source("format-mammal-data.R")

#-------------------------------------------------------------------------------#
# Visualize when and where cameras were deployed
#-------------------------------------------------------------------------------#

# Create a park-specific index for camera location (1-65)
events <- events %>% 
  group_by(Park) %>% 
  mutate(locnum = as.numeric(as.factor(StdLocName))) %>%
  as.data.frame()

# Plot events at three main parks, excluding 2016 (different locs sampled at CHIR)
ggplot() + 
  geom_segment(filter(events, Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 0.3, color = "dodgerblue3") +
  geom_segment(filter(events, Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = d_date, y = locnum+0.4, yend = locnum-0.4),
               size = 0.3, color = "dodgerblue3") + 
  geom_segment(filter(events, Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = r_date, xend = r_date, y = locnum+0.4, yend = locnum-0.4),
               size = 0.3, color = "dodgerblue3") + 
  labs(x = 'Date', y = "Camera number") + 
  facet_grid(rows = vars(Park))
# Save plot in output/folder
# ggsave("output/SamplingEvents_3Parks.jpg",
#        width = 6.5, height = 6.5,
#        units = "in")

# Plot events at four smaller parks
ggplot() + 
  geom_segment(filter(events, !Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
               size = 0.5, color = "dodgerblue3") +
  geom_segment(filter(events, !Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = d_date, xend = d_date, y = locnum+0.4, yend = locnum-0.4),
               size = 0.3, color = "dodgerblue3") + 
  geom_segment(filter(events, !Park %in% c("CHIR", "ORPI", "SAGW") & d_yr > 2016),
               mapping = aes(x = r_date, xend = r_date, y = locnum+0.4, yend = locnum-0.4),
               size = 0.3, color = "dodgerblue3") + 
  labs(x = 'Date', y = "Camera number") + 
  facet_grid(rows = vars(Park))
# Save plot in output/folder
# ggsave("output/SamplingEvents_OtherParks.jpg",
#        width = 6.5, height = 6.5,
#        units = "in")

#-------------------------------------------------------------------------------#
# Visualize deployments and photos
#-------------------------------------------------------------------------------#

# Functions to create histograms with photo dates and charts with deployments
ggcust_hist <- function(dat, park, year, xlims){
  dd <- filter(dat, Park == park & yr == year)
  ggplot(dd, aes(obsdate)) + 
    geom_histogram(binwidth = 1) +
    labs(x = "", y = "Number of photos") + 
    coord_cartesian(xlim = xlims) + 
    theme(text=element_text(size = 8),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin = margin(0.2, 0.2, 0, 0.2, unit = "cm"))
}

ggcust_segs <- function(dat, park, year, xlims){
  ggplot() +
    geom_segment(filter(dat, Park == park & d_yr == year),
                 mapping = aes(x = d_date, xend = r_date, y = locnum, yend = locnum),
                 size = 0.3, color = "dodgerblue3") + 
    geom_segment(filter(dat, Park == park & d_yr == year),
                 mapping = aes(x = d_date, xend = d_date, y = locnum+0.4, yend = locnum-0.4),
                 size = 0.3, color = "dodgerblue3") + 
    geom_segment(filter(dat, Park == park & d_yr == year),
                 mapping = aes(x = r_date, xend = r_date, y = locnum+0.4, yend = locnum-0.4),
                 size = 0.3, color = "dodgerblue3") + 
    labs(x = "", y = "Camera number") + 
    coord_cartesian(xlim = xlims) + 
    theme(text=element_text(size = 8),
          axis.title.x=element_blank(),
          plot.margin = margin(0.1, 0.2, 0.2, 0.2, unit = "cm"))
}


# SAGW
s17h <- ggcust_hist(dat, "SAGW", 2017, xlims = as.Date(c("2017-01-01", "2017-03-31"))) +
  annotate("text", x = as.Date("2017-03-31"), y = Inf, hjust=1, vjust=1.5, label = "2017", size = 3.5)
s17s <- ggcust_segs(events, "SAGW", 2017, xlims = as.Date(c("2017-01-01", "2017-03-31")))

s18h <- ggcust_hist(dat, "SAGW", 2018, xlims = as.Date(c("2018-01-01", "2018-03-31"))) +
  annotate("text", x = as.Date("2018-03-31"), y = Inf, hjust=1, vjust=1.5, label = "2018", size = 3.5)
s18s <- ggcust_segs(events, "SAGW", 2018, xlims = as.Date(c("2018-01-01", "2018-03-31")))

s20h <- ggcust_hist(dat, "SAGW", 2020, xlims = as.Date(c("2020-01-01", "2020-03-31"))) +
  annotate("text", x = as.Date("2020-03-31"), y = Inf, hjust=1, vjust=1.5, label = "2020", size = 3.5)
s20s <- ggcust_segs(events, "SAGW", 2020, xlims = as.Date(c("2020-01-01", "2020-03-31")))

s21h <- ggcust_hist(dat, "SAGW", 2021, xlims = as.Date(c("2021-01-01", "2021-03-31"))) +
  annotate("text", x = as.Date("2021-03-31"), y = Inf, hjust=1, vjust=1.5, label = "2021", size = 3.5)
s21s <- ggcust_segs(events, "SAGW", 2021, xlims = as.Date(c("2021-01-01", "2021-03-31")))

s22h <- ggcust_hist(dat, "SAGW", 2022, xlims = as.Date(c("2022-01-01", "2022-03-31"))) +
  annotate("text", x = as.Date("2022-03-31"), y = Inf, hjust=1, vjust=1.5, label = "2022", size = 3.5)
s22s <- ggcust_segs(events, "SAGW", 2022, xlims = as.Date(c("2022-01-01", "2022-03-31")))

s17 <- grid.arrange(s17h, s17s, nrow = 2)
s18 <- grid.arrange(s18h, s18s, nrow = 2)
s20 <- grid.arrange(s20h, s20s, nrow = 2)
s21 <- grid.arrange(s21h, s21s, nrow = 2)
s22 <- grid.arrange(s22h, s22s, nrow = 2)

sALL <- grid.arrange(s17, s18, s20, s21, s22, nrow = 2)
# Save plot in output/folder
# ggsave("output/Events&Obs_SAGW.jpg",
#        sALL,
#        width = 6.5, height = 6.5,
#        units = "in")

# ORPI
o17h <- ggcust_hist(dat, "ORPI", 2017, xlims = as.Date(c("2017-03-01", "2017-11-30"))) +
  annotate("text", x = as.Date("2017-11-30"), y = Inf, hjust=1, vjust=1.5, label = "2017", size = 3.5)
o17s <- ggcust_segs(events, "ORPI", 2017, xlims = as.Date(c("2017-03-01", "2017-11-30")))

o18h <- ggcust_hist(dat, "ORPI", 2018, xlims = as.Date(c("2018-03-01", "2018-11-30"))) +
  annotate("text", x = as.Date("2018-11-30"), y = Inf, hjust=1, vjust=1.5, label = "2018", size = 3.5)
o18s <- ggcust_segs(events, "ORPI", 2018, xlims = as.Date(c("2018-03-01", "2018-11-30")))

o20h <- ggcust_hist(dat, "ORPI", 2020, xlims = as.Date(c("2020-03-01", "2020-11-30"))) +
  annotate("text", x = as.Date("2020-11-30"), y = Inf, hjust=1, vjust=1.5, label = "2020", size = 3.5)
o20s <- ggcust_segs(events, "ORPI", 2020, xlims = as.Date(c("2020-03-01", "2020-11-30")))

o21h <- ggcust_hist(dat, "ORPI", 2021, xlims = as.Date(c("2021-03-01", "2021-11-30"))) +
  annotate("text", x = as.Date("2021-11-30"), y = Inf, hjust=1, vjust=1.5, label = "2021", size = 3.5)
o21s <- ggcust_segs(events, "ORPI", 2021, xlims = as.Date(c("2021-03-01", "2021-11-30")))

o17 <- grid.arrange(o17h, o17s, nrow = 2)
o18 <- grid.arrange(o18h, o18s, nrow = 2)
o20 <- grid.arrange(o20h, o20s, nrow = 2)
o21 <- grid.arrange(o21h, o21s, nrow = 2)

oALL <- grid.arrange(o17, o18, o20, o21, nrow = 2)
# Save plot in output/folder
# ggsave("output/Events&Obs_ORPI.jpg",
#        oALL,
#        width = 6.5, height = 6.5,
#        units = "in")

# CHIR
c17h <- ggcust_hist(dat, "CHIR", 2017, xlims = as.Date(c("2017-01-01", "2017-12-31"))) +
  annotate("text", x = as.Date("2017-12-31"), y = Inf, hjust=1, vjust=1.5, label = "2017", size = 3.5)
c17s <- ggcust_segs(events, "CHIR", 2017, xlims = as.Date(c("2017-01-01", "2017-12-31")))

c18h <- ggcust_hist(dat, "CHIR", 2018, xlims = as.Date(c("2018-01-01", "2018-12-31"))) +
  annotate("text", x = as.Date("2018-12-31"), y = Inf, hjust=1, vjust=1.5, label = "2018", size = 3.5)
c18s <- ggcust_segs(events, "CHIR", 2018, xlims = as.Date(c("2018-01-01", "2018-12-31")))

#Need to adjust start dates for 2019, since cameras were up before Jan 1
c19h <- ggcust_hist(dat, "CHIR", 2019, xlims = as.Date(c("2019-01-01", "2019-12-31"))) +
  annotate("text", x = as.Date("2019-12-31"), y = Inf, hjust=1, vjust=1.5, label = "2019", size = 3.5)
c19s <- ggplot() +
  geom_segment(filter(events, Park == "CHIR" & d_yr == 2019),
               mapping = aes(x = as.Date("2019-01-01"), xend = r_date, 
                             y = locnum, yend = locnum),
               size = 0.3, color = "dodgerblue3") +   
  geom_segment(filter(events, Park == "CHIR" & d_yr == 2019),
               mapping = aes(x = d_date, xend = d_date, y = locnum+0.4, yend = locnum-0.4),
               size = 0.3, color = "dodgerblue3") + 
  geom_segment(filter(events, Park == "CHIR" & d_yr == 2019),
               mapping = aes(x = r_date, xend = r_date, y = locnum+0.4, yend = locnum-0.4),
               size = 0.3, color = "dodgerblue3") +   
  labs(x = "", y = "Camera number") + 
  coord_cartesian(xlim = as.Date(c("2019-01-01", "2019-12-31"))) + 
  theme(text=element_text(size = 8),
        axis.title.x=element_blank(),
        plot.margin = margin(0.1, 0.2, 0.2, 0.2, unit = "cm"))

c21h <- ggcust_hist(dat, "CHIR", 2021, xlims = as.Date(c("2021-01-01", "2021-12-31"))) +
  annotate("text", x = as.Date("2021-12-31"), y = Inf, hjust=1, vjust=1.5, label = "2021", size = 3.5)
c21s <- ggcust_segs(events, "CHIR", 2021, xlims = as.Date(c("2021-01-01", "2021-12-31")))

c17 <- grid.arrange(c17h, c17s, nrow = 2)
c18 <- grid.arrange(c18h, c18s, nrow = 2)
c19 <- grid.arrange(c19h, c19s, nrow = 2)
c21 <- grid.arrange(c21h, c21s, nrow = 2)

cALL <- grid.arrange(c17, c18, c19, c21, nrow = 2)
# Save plot in output/folder
# ggsave("output/Events&obs_CHIR.jpg",
#        cALL,
#        width = 6.5, height = 6.5,
#        units = "in")

#-------------------------------------------------------------------------------#
# Detections by species
#-------------------------------------------------------------------------------#

(species <- arrange(species, desc(n)))

park_yr <- dat %>%
  group_by(Park, yr) %>%
  summarize(all_spp = length(Species_code)) %>%
  as.data.frame()

spp_park_yr <- dat %>%
  group_by(Park, yr, Species_code) %>%
  summarize(n_photos = length(Species_code)) %>%
  as.data.frame()

spp_park_yr <- left_join(spp_park_yr, park_yr)
spp_park_yr <- mutate(spp_park_yr, prop_photos = n_photos / all_spp)
spp_park_yr <- arrange(spp_park_yr, Park, yr, desc(prop_photos))

# Tables for SAGW, where cameras deployed during same period each year
# Number of photos
n_sagw_yr <- spp_park_yr %>% 
  filter(Park == "SAGW") %>%
  select(yr, Species_code, n_photos) %>%
  spread(key = yr, value = n_photos, fill = 0)
oldnames <- names(n_sagw_yr)[2:6]
n_sagw_yr <-n_sagw_yr %>% 
  rename_with(~ paste0("n_",oldnames), all_of(oldnames)) %>%
  arrange(desc(n_2017))
# Proportion of photos
p_sagw_yr <- spp_park_yr %>% 
  filter(Park == "SAGW") %>%
  select(yr, Species_code, prop_photos) %>%
  spread(key = yr, value = prop_photos, fill = 0) %>%
  mutate(across(2:6, round, 3))
oldnames <- names(p_sagw_yr)[2:6]
p_sagw_yr <- p_sagw_yr %>% 
  rename_with(~ paste0("p_",oldnames), all_of(oldnames)) %>%
  arrange(desc(p_2017))

spp_sagw_yr <- left_join(n_sagw_yr, p_sagw_yr)
spp_sagw_yr
rm(n_sagw_yr, p_sagw_yr)

# Table for ORPI, mostly deployments in April-ish (but stayed out longer in 2020)
# Number of photos
n_orpi_yr <- spp_park_yr %>% 
  filter(Park == "ORPI") %>%
  select(yr, Species_code, n_photos) %>%
  spread(key = yr, value = n_photos, fill = 0)
oldnames <- names(n_orpi_yr)[2:6]
n_orpi_yr <-n_orpi_yr %>% 
  rename_with(~ paste0("n_",oldnames), all_of(oldnames)) %>%
  arrange(desc(n_2017))
# Number of photos
p_orpi_yr <- spp_park_yr %>% 
  filter(Park == "ORPI") %>%
  select(yr, Species_code, prop_photos) %>%
  spread(key = yr, value = prop_photos, fill = 0) %>%
  mutate(across(2:6, round, 3))
oldnames <- names(p_orpi_yr)[2:6]
p_orpi_yr <- p_orpi_yr %>% 
  rename_with(~ paste0("p_",oldnames), all_of(oldnames)) %>%
  arrange(desc(p_2017))

spp_orpi_yr <- left_join(n_orpi_yr, p_orpi_yr)
spp_orpi_yr
rm(n_orpi_yr, p_orpi_yr)

# Table for CHIR (but note that deployments are at diff times of year)
# Number of photos
n_chir_yr <- spp_park_yr %>% 
  filter(Park == "CHIR" & yr != 2020) %>%
  select(yr, Species_code, n_photos) %>%
  spread(key = yr, value = n_photos, fill = 0)
oldnames <- names(n_chir_yr)[2:5]
n_chir_yr <-n_chir_yr %>% 
  rename_with(~ paste0("n_",oldnames), all_of(oldnames)) %>%
  arrange(desc(n_2017))
# Number of photos
p_chir_yr <- spp_park_yr %>% 
  filter(Park == "CHIR" & yr != 2020) %>%
  select(yr, Species_code, prop_photos) %>%
  spread(key = yr, value = prop_photos, fill = 0) %>%
  mutate(across(2:5, round, 3))
oldnames <- names(p_chir_yr)[2:5]
p_chir_yr <- p_chir_yr %>% 
  rename_with(~ paste0("p_",oldnames), all_of(oldnames)) %>%
  arrange(desc(p_2017))

spp_chir_yr <- left_join(n_chir_yr, p_chir_yr)
spp_chir_yr
rm(n_chir_yr, p_chir_yr)

#-------------------------------------------------------------------------------#
# Calculate time between photos of the same species at same location
#-------------------------------------------------------------------------------#
 
dat <- arrange(dat, Park, loc_short, Species_code, yr, obsdate, obstime)
dat$time_btw <- NA
for (i in 2:nrow(dat)) {
  if(dat$Park[i] == dat$Park[i-1] &
     dat$loc_short[i] == dat$loc_short[i-1] & 
     dat$Species_code[i] == dat$Species_code[i-1] &
     dat$yr[i] == dat$yr[i-1]) {
    dat$time_btw[i] <- as.double(difftime(dat$datetime[i], 
                                          dat$datetime[i-1], 
                                          units = "mins"))
  } else {
    dat$time_btw[i] <- NA
  }
}
# checks:
head(dat[,c(4:5,8,14,20)], 20)
tail(dat[,c(4:5,8,14,20)], 20)

# species average
tb_spp <- dat %>%
  group_by(Species_code) %>%
  summarize(n_diffs = sum(!is.na(time_btw)),
            sec_btw_md = round(median(time_btw, na.rm = TRUE),2),
            secs_btw_mn = round(mean(time_btw, na.rm = TRUE),2)) %>%
  arrange(desc(n_diffs)) %>%
  as.data.frame()
tb_spp

