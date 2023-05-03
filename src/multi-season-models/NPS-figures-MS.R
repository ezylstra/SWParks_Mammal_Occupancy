################################################################################
# Prepare plots for publication in NPS format from from single-season occupancy models
#    marginal effect of covariates on occupancy probability
#    marginal effect of covariates on detection probability    
#    maps of occupancy probability across a park
#
# CL McIntyre
# Updated 2023-04-21
################################################################################

#library(grDevices) #grDevices should automatically load

# this script will save figures for all occupancy and detection co-variates and
# maps of occupancy probability for the first and last year 
# if you don't want to save all of them, you can tweak the code

#------------------------------------------------------------------------------#
# Create custom theme for plots used in NPS publications (including web)
#------------------------------------------------------------------------------#

# use Frutiger font - NPS standard - whenever possible
# NPS employees can download Fruiter to their NPS computer from https://www.nps.gov/subjects/hfc/nps-typefaces.htm
# must be on the vpn to do the download

# add the font to R
windowsFonts("Frutiger LT Std 55 Roman" = windowsFont("Frutiger LT Std 55 Roman"))

# make custom theme using Frutiger and defining font sizes
theme_NPS <- ggplot2::theme_classic() + 
  theme(legend.title = element_text(size=10, color="black")) +
  theme(legend.text = element_text(size=9,color="black")) +
  theme(axis.title = element_text(size=10,color="black")) + 
  theme(axis.text = element_text(size=9,color="black")) +
  theme(plot.title = element_text(size=12, color="black", hjust = 0.5)) +
  theme(axis.ticks = element_line(color='black')) + 
  theme(axis.line = element_line(color='black')) +
  theme(plot.subtitle = element_text(size=10, color="black", hjust = 0.5)) +
  theme(text = element_text(family = "Frutiger LT Std 55 Roman", face="plain"))

# get longer park name for use in plots
park_name <- ifelse(PARK=="CHIR","Chiricahua NM",ifelse(PARK=="SAGW","Saguaro NP","Organ Pipe Cactus NM"))


#------------------------------------------------------------------------------#
# Apply custom theme to trend plots (other co-variates at their means) 
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

if ("years" %in% psi_covs) {
  trend_NPS <- trend
  title <- paste("Change over time in proportion of area used by ",species$Common_name[species$Species_code==SPECIES], sep="")
  subtitle <- paste(park_name, ", ", min(YEARS),"-", max(YEARS), sep="")
  trend_NPS <- trend_NPS + theme_NPS + ggtitle(str_wrap(title, 60), subtitle)
  fig_file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"-", max(YEARS),"_trend_4NPS.pdf",sep="")
  ggsave(trend_NPS, file = fig_file, device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  print(trend_NPS)
}


#------------------------------------------------------------------------------#
# Apply custom theme to marginal effect of covariates on occupancy probability 
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

for (fig in str_subset(ls(), "marginal_psi_")) {
  fig_NPS <- get(fig)
  title <- paste("Occurrence probability of ",species$Common_name[species$Species_code==SPECIES], " vs. ", tolower(fig_NPS$labels[[1]]), sep="")
  subtitle <- paste(park_name, ", ", min(YEARS),"-", max(YEARS), sep="")
  fig_NPS <- fig_NPS + theme_NPS + ggtitle(str_wrap(title, 60), subtitle)
  fig_name <- fig
  fig_file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"-", max(YEARS),"_",fig_name,"_4NPS.pdf",sep="")
  ggsave(fig_NPS, file = fig_file, device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  print(fig_NPS)
}

#------------------------------------------------------------------------------#
# Apply custom theme to marginal effect of covariates on detection probability 
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

for (fig in str_subset(ls(), "marginal_p_")) {
  fig_NPS <- get(fig)
  title <- paste("Detection probability of ",species$Common_name[species$Species_code==SPECIES], " vs. ", tolower(fig_NPS$labels[[1]]), sep="")
  subtitle <- paste(park_name, ", ", min(YEARS),"-", max(YEARS), sep="")
  fig_NPS <- fig_NPS + theme_NPS + ggtitle(str_wrap(title, 60), subtitle)
  fig_name <- fig
  fig_file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"-", max(YEARS),"_",fig_name,"_4NPS.pdf",sep="")
  ggsave(fig_NPS, file = fig_file, device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  print(fig_NPS)
}

#------------------------------------------------------------------------------#
# Apply custom theme to spatial visualizations of mean occurrence probability
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

if (length(psi_spatcovs) > 0) { 
  mn_title <-  paste("Mean occurrence probability of ",species$Common_name[species$Species_code==SPECIES], sep="")
  subtitle_lastyr <- paste(park_name, ", ", max(YEARS), sep="")
  subtitle_firstyr <- paste(park_name, ", ", min(YEARS), sep="")
  
  # reassign min and max values to be the same for first year and last year
  # so that scales are the same for both figures
  preds_mn_min <- min(min(preds_mn_firstyr[], na.rm=TRUE), min(preds_mn_lastyr[], na.rm=TRUE))
  preds_mn_max <- max(max(preds_mn_firstyr[], na.rm=TRUE), max(preds_mn_lastyr[], na.rm=TRUE))
  
  # view and save spatial prediction for first  year as individual file 
   plot_preds_firstyr_NPS <- ggplot() + 
    geom_spatraster(data = preds_mn_firstyr, mapping = aes(fill = mean_firstyr)) + 
    scale_fill_viridis_c(na.value = 'transparent', limits=c(preds_mn_min,preds_mn_max)) + 
    labs(fill = '', title = mn_title, subtitle = subtitle_firstyr) +
    theme_NPS +
    theme(axis.line = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank())
  print(plot_preds_firstyr_NPS)
  ggsave(plot_preds_firstyr_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"_","MeanOccurrenceMap_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  
  # view and save spatial prediction for first  year as individual file 
  plot_preds_lastyr_NPS <- ggplot() + 
    geom_spatraster(data = preds_mn_lastyr, mapping = aes(fill = mean_lastyr)) + 
    scale_fill_viridis_c(na.value = 'transparent', limits=c(preds_mn_min,preds_mn_max)) +
    labs(fill = '', title = mn_title, subtitle = subtitle_lastyr) +
    theme_NPS +
    theme(axis.line = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank())
  print(plot_preds_lastyr_NPS)
  ggsave(plot_preds_lastyr_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",max(YEARS),"_","MeanOccurrenceMap_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  
  # with both years stacked in a single figure
  plot_preds_firstlast_NPS <- grid.arrange(plot_preds_firstyr_NPS, plot_preds_lastyr_NPS+theme(plot.title=element_blank()), nrow = 2)
  print(plot_preds_firstlast_NPS)
  ggsave(plot_preds_firstlast_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"-",max(YEARS),"_","MeanOccurrenceMaps_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 6, units="in")
}

