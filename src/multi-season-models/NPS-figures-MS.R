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
park <- ifelse(PARK=="CHIR","Chiricahua NM",ifelse(PARK=="SAGW","Saguaro NP","Organ Pipe Cactus NM"))
SpeciesName <- ifelse()

#------------------------------------------------------------------------------#
# Apply custom theme to trend plots (other co-variates at their means) 
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

if ("years" %in% psi_covs) {
  trend_NPS <- trend
  title <- paste("Change over time in proportion of area used by ",species$Common_name[species$Species_code==SPECIES], sep="")
  subtitle <- paste(park, ", ", min(YEARS),"-", max(YEARS), sep="")
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
  subtitle <- paste(park, ", ", min(YEARS),"-", max(YEARS), sep="")
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
  subtitle <- paste(park, ", ", min(YEARS),"-", max(YEARS), sep="")
  fig_NPS <- fig_NPS + theme_NPS + ggtitle(str_wrap(title, 60), subtitle)
  fig_name <- fig
  fig_file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"-", max(YEARS),"_",fig_name,"_4NPS.pdf",sep="")
  ggsave(fig_NPS, file = fig_file, device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  print(fig_NPS)
}

#------------------------------------------------------------------------------#
# Apply custom theme to visualizations of occurrence probability (mean + SD)
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

mn_title <-  paste("Mean occurrence probability of ",species$Common_name[species$Species_code==SPECIES], sep="")
subtitle_lastyr <- paste(park, ", ", max(YEARS), sep="")
subtitle_firstyr <- paste(park, ", ", min(YEARS), sep="")
  
# view and save spatial prediction for first year
if (length(psi_spatcovs) > 0) { 
  plot_preds_firstyr_NPS <- plot_preds_mn_firstyr + 
    #geom_spatraster(data = preds_mn, mapping = aes(fill = mean)) + 
    #scale_fill_viridis_c(na.value = 'transparent') +
    #geom_spatvector(data=park_boundary, fill=NA, color="black", size=0.5) + 
    labs(fill = '', title = mn_title, subtitle = subtitle_firstyr) +
    theme_NPS +
    theme(axis.line = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank())
  plot_preds_firstyr_NPS
  ggsave(plot_preds_firstyr_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",min(YEARS),"_","MeanOccurrenceMap_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
}

# view and save spatial prediction for last year
if (length(psi_spatcovs) > 0) { 
  plot_preds_lastyr_NPS <- plot_preds_mn_lastyr + 
    #geom_spatraster(data = preds_mn, mapping = aes(fill = mean)) + 
    #scale_fill_viridis_c(na.value = 'transparent') +
    #geom_spatvector(data=park_boundary, fill=NA, color="black", size=0.5) + 
    labs(fill = '', title = mn_title, subtitle = subtitle_lastyr) +
    theme_NPS +
    theme(axis.line = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank())
  plot_preds_lastyr_NPS
  ggsave(plot_preds_lastyr_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",max(YEARS),"_","MeanOccurrenceMap_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
}
