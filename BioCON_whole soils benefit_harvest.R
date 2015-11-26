library(tidyr)
library(dplyr)

setwd('C:\\Users\\Kim\\Dropbox\\2015_NSF_LaPierre\\data\\La Pierre and Simms data\\2015 whole soil benefit')

###import data and planned replicates
harvest <- read.csv('La Pierre_2015_CDR biocon_whole soil_harvest_corrected.csv')

plan <- read.csv('La Pierre_2015_CDR biocon_whole soil_plots.csv')%>%
  gather(key=replicate, value=value, rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8)%>%
  select(-replicate)
names(plan)[names(plan)=="value"] <- "replicate"

trt <- read.csv('e141_treatments.csv')

#merge plan and harvest data to see what reps are missing (i.e., didn't grow in the field); also used for data check
harvestAll <- merge(plan, harvest, by=c('plot', 'species', 'replicate'), all=T)

#merge harvest data with treatment information
harvestTrt <- merge(harvestAll, trt, by=c('ring', 'plot'), all=T)%>%
  filter(year!='NA')

###calculate the relative biomass compared to uninoculated controls (percent difference)

#subset out the soil inoculation controls and inoculated plants
harvestCtl <- harvestTrt%>%filter(is.na(experiment))
names(harvestCtl)[names(harvestCtl)=='shoot_mass'] <- 'shoot_mass_ctl'
names(harvestCtl)[names(harvestCtl)=='root_mass'] <- 'root_mass_ctl'
harvestCtlShoot <- harvestCtl%>%select(-root_mass_ctl)%>%
  group_by(plot, species)%>%
  summarise(shoot_mass_ctl=mean(shoot_mass_ctl))
  harvestCtlShoot$CO2_trt <- with(harvestCtlShoot, ifelse(plot=='aCO2aN', 'Camb', ifelse(plot=='aCO2eN', 'Camb', 'Cenrich')))
  harvestCtlShoot$N_trt <- with(harvestCtlShoot, ifelse(plot=='aCO2aN', 'Namb', ifelse(plot=='eCO2aN', 'Namb', 'Nenrich')))
harvestCtlRoot <- harvestCtl%>%select(-shoot_mass_ctl)%>%
  group_by(plot, species)%>%
  summarise(root_mass_ctl=mean(root_mass_ctl))
  harvestCtlRoot$CO2_trt <- with(harvestCtlRoot, ifelse(plot=='aCO2aN', 'Camb', ifelse(plot=='aCO2eN', 'Camb', 'Cenrich')))
  harvestCtlRoot$N_trt <- with(harvestCtlRoot, ifelse(plot=='aCO2aN', 'Namb', ifelse(plot=='eCO2aN', 'Namb', 'Nenrich')))
harvestCtlAll <- merge(harvestCtlShoot, harvestCtlRoot, by=c('plot', 'species', 'CO2_trt', 'N_trt'))%>%
  select(-plot)

#subset out inoculated plants
harvestIno <- harvestTrt%>%filter(!is.na(experiment))

#merge uninoculated controls back with inoculated plants, calculate relative cover, drop unneccesary columns
harvestRel <- merge(harvestIno, harvestCtlAll, by=c('CO2_trt', 'N_trt', 'species'), all=T)%>%
  mutate(shoot_mass_rel=((shoot_mass-shoot_mass_ctl)/shoot_mass_ctl), root_mass_rel=((root_mass-root_mass_ctl)/root_mass_ctl))
  
  
  
  
  
  
  
  
  


