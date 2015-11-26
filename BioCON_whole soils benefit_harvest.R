library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\2015_NSF_LaPierre\\data\\La Pierre and Simms data\\2015 whole soil benefit')

###import data and planned replicates
harvest <- read.csv('La Pierre_2015_CDR biocon_whole soil_harvest_corrected.csv')

plan <- read.csv('La Pierre_2015_CDR biocon_whole soil_plots.csv')%>%
  gather(key=replicate, value=value, rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8)%>%
  select(-replicate)
names(plan)[names(plan)=="value"] <- "replicate"

#merge plan and harvest data to see what reps are missing (i.e., didn't grow in the field); also used for data check
harvestAll <- merge(plan, harvest, by=c('plot', 'species', 'replicate'), all=T)%>%
  group_by(plot, species, replicate)%>%
  summarise(length=length(shoot_mass))
