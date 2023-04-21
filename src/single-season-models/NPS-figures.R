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


#------------------------------------------------------------------------------#
# Apply custom theme to marginal effect of covariates on occupancy probability 
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

for (fig in str_subset(ls(), "marginal_psi_")) {
  fig_NPS <- get(fig)
  title <- paste("Occupancy probability of ",species$Common_name[species$Species_code==SPECIES], " vs. ", fig_NPS$labels[[1]], sep="")
  subtitle <- paste(park, ", ", YEAR, sep="")
  fig_NPS <- fig_NPS + theme_NPS + ggtitle(str_wrap(title, 60), subtitle)
  fig_name <- fig
  fig_file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",YEAR,"_",fig_name,"_4NPS.pdf",sep="")
  ggsave(fig_NPS, file = fig_file, device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  print(fig_NPS)
}

#------------------------------------------------------------------------------#
# Apply custom theme to marginal effect of covariates on detection probability 
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

for (fig in str_subset(ls(), "marginal_p_")) {
  fig_NPS <- get(fig)
  title <- paste("Detection probability of ",species$Common_name[species$Species_code==SPECIES], " vs. ", fig_NPS$labels[[1]], sep="")
  subtitle <- paste(park, ", ", YEAR, sep="")
  fig_NPS <- fig_NPS + theme_NPS + ggtitle(str_wrap(title, 60), subtitle)
  fig_name <- fig
  fig_file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",YEAR,"_",fig_name,"_4NPS.pdf",sep="")
  ggsave(fig_NPS, file = fig_file, device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
  print(fig_NPS)
}

#------------------------------------------------------------------------------#
# Apply custom theme to visualizations of occurrence probability (mean + SD)
# and save as 4x6" PDFs
#------------------------------------------------------------------------------#

# Use tidyterra to create plots using ggplot syntax
mn_title <-  paste("Mean occurrence probability of ",species$Common_name[species$Species_code==SPECIES], sep="")
sd_title <-  paste("Standard deviation of occurrence probability of ",species$Common_name[species$Species_code==SPECIES], sep="")
subtitle <- paste(park, ", ", YEAR, sep="")
  
plot_preds_mn_NPS <- ggplot() + 
  geom_spatraster(data = preds_mn, mapping = aes(fill = mean)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  #geom_spatvector(data=park_boundary, fill=NA, color="black", size=0.5) + 
  labs(fill = '', title = mn_title, subtitle = subtitle) +
  theme_NPS +
  theme(axis.line = element_blank()) + 
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank())
plot_preds_mn_NPS
ggsave(plot_preds_mn_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",YEAR,"_","MeanOccurrenceMap_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")

plot_preds_sd_NPS <- ggplot() + 
  geom_spatraster(data = preds_sd, mapping = aes(fill = sd)) + 
  scale_fill_viridis_c(na.value = 'transparent') +
  #geom_spatvector(data=park_boundary, fill=NA, color="black", size=0.5) + 
  labs(fill = '', title = sd_title, subtitle = subtitle) +
  theme_NPS +
  theme(axis.line = element_blank()) + 
  theme(axis.text = element_blank()) + 
  theme(axis.ticks = element_blank())
plot_preds_sd_NPS
ggsave(plot_preds_sd_NPS, file = paste(getwd(),"/output/NPS-figures/",PARK,"-",SPECIES,"-",YEAR,"_","sdOccurrenceMap_4NPS.pdf",sep=""), device = cairo_pdf, dpi=300, width = 6, height = 4, units="in")
