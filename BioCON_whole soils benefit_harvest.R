library(lme4)
library(ggplot2)
library(tidyr)
library(dplyr)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

setwd('C:\\Users\\Kim\\Dropbox\\2015_NSF_LaPierre\\data\\La Pierre and Simms data\\2015 whole soil benefit')

###import data and planned replicates
harvest <- read.csv('La Pierre_2015_CDR biocon_whole soil_harvest_corrected.csv')%>%
  #calculate total nodule number
  mutate(total_nod = pink_nod + white_nod)%>%
  select(-notes)

plan <- read.csv('La Pierre_2015_CDR biocon_whole soil_plots.csv')%>%
  gather(key=replicate, value=value, rep1, rep2, rep3, rep4, rep5, rep6, rep7, rep8)%>%
  select(-replicate)
names(plan)[names(plan)=="value"] <- "replicate"

trt <- read.csv('e141_treatments.csv')%>%
  #get total number of legumes in each plot
  mutate(tot_legume=Amorpha.canescens+Lespedeza.capitata+Lupinus.perennis+Petalostemum.villosum)

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
  mutate(shoot_mass_rel=((shoot_mass-shoot_mass_ctl)/shoot_mass_ctl), root_mass_rel=((root_mass-root_mass_ctl)/root_mass_ctl))%>%
  select(-shoot_mass_ctl, -root_mass_ctl)%>%
  #average across replicates
  group_by(CO2_trt, N_trt, species, ring, plot, year, spp_count, group_count, experiment, monospecies, monogroup, Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, C.3, C.4, Forb, Legume, tot_legume)%>%
  summarise(shoot_mass=mean(shoot_mass), root_mass=mean(root_mass), pink_nod=mean(pink_nod), white_nod=mean(white_nod), total_nod=mean(total_nod), shoot_mass_rel=mean(shoot_mass_rel), root_mass_rel=mean(root_mass_rel))


###soil feedbacks under control conditions
#subset out CO2 and N treatments
#separate models for each species

#shoot mass
summary(lmer(shoot_mass_rel ~ Lespedeza.capitata + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(shoot_mass_rel ~ Lupinus.perennis + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#root mass
summary(lmer(root_mass_rel ~ Lespedeza.capitata + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(root_mass_rel ~ Lupinus.perennis + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#pink nodules
summary(lmer(pink_nod ~ Lespedeza.capitata + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(pink_nod ~ Lupinus.perennis + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#white nodules
summary(lmer(white_nod ~ Lespedeza.capitata + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(white_nod ~ Lupinus.perennis + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#total nodules
summary(lmer(total_nod ~ Lespedeza.capitata + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(total_nod ~ Lupinus.perennis + tot_legume + (1|spp_count), data=subset(harvestRel, species=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#plot results
harvestRel$orderLECA <- factor(harvestRel$Lespedeza.capitata, levels=c('0', '1'))
harvestRel$orderLUPE <- factor(harvestRel$Lupinus.perennis, levels=c('0', '1'))
harvestRel$orderTotLeg <- factor(harvestRel$tot_legume, levels=c('1', '4'))

feedbackLECAfig <- ggplot(barGraphStats(data=subset(harvestRel, species=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich'), variable='shoot_mass_rel', byFactorNames=c('orderTotLeg', 'orderLECA')), aes(x=orderTotLeg, y=mean, fill=orderLECA)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Number of legumes in soil') + ylab('Relative shoot mass (g)')
  
feedbackLUPEfig <- ggplot(barGraphStats(data=subset(harvestRel, species=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich'), variable='shoot_mass_rel', byFactorNames=c('orderTotLeg', 'orderLUPE')), aes(x=orderTotLeg, y=mean, fill=orderLUPE)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Number of legumes in soil') + ylab('Relative shoot mass (g)')

###quick ANOVA for CO2 and N trts
summary(response <- lmer(shoot_mass ~ CO2_trt + N_trt + CO2_trt*N_trt + species + species*CO2_trt + species*N_trt + species*CO2_trt*N_trt + (1|ring), data=harvestRel))
















