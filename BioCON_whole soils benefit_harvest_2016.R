library(lme4)
library(boot)
library(ggplot2)
library(grid)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

#cbind fill function
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

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

setwd('C:\\Users\\Kim\\Dropbox\\NSF BioCON rhizobia\\data\\La Pierre and Simms data\\2016 whole soil benefit')

###import data and planned replicates
#envelope data
envelope_plant <- read.csv('La Pierre_2016_whole soil_harvest.csv')

#split out data on number of plants
plants <- envelope_plant%>%
  select(pot, num_plants, notes)

#split out envelopes
envelope <- envelope_plant%>%
  select(pot, shoot_env, root_env)%>%
  filter(!is.na(root_env))%>%
  gather(key=root.shoot_harvest, value=envelope, shoot_env, root_env, na.rm=F)%>%
  mutate(root.shoot_harvest=ifelse(root.shoot_harvest=='shoot_env', 's', 'r'))

#pot data
plan <- read.csv('La Pierre_2016_whole soil_pot setup.csv')%>%
  select(-growth_06132016)

#biomass data
harvest <- read.csv('La Pierre_2016_BioCON_biomass.csv')

#split out nodule number
nodules <- harvest%>%
  select(envelope, f_nodules, n_nodules)%>%
  #calculate total nodule number and total biomass
  mutate(total_nod=f_nodules+n_nodules)%>%
  filter(!is.na(total_nod))%>%
  merge(envelope, by=c('envelope'))%>%
  select(f_nodules, n_nodules, total_nod, pot)

#split out biomass
biomass <- harvest%>%
  select(envelope, mass_g, root.shoot, notes)

#trt data
trt <- read.csv('e141_treatments.csv')%>%
  #get total number of legumes in each plot
  mutate(tot_legume=Amorpha.canescens+Lespedeza.capitata+Lupinus.perennis+Petalostemum.villosum)


#merge pot, biomass, nodule, and number of plants data
harvestAll <- merge(envelope, biomass, by=c('envelope'), all=T)%>%
  # mutate(check=ifelse(root.shoot==root.shoot_harvest, 1, 0))%>%
  filter(pot!='is.na')%>%
  select(pot, mass_g, root.shoot_harvest)%>%
  spread(key=root.shoot_harvest, value=mass_g)%>%
  #calculate total mass and rename root and shoot mass
  mutate(mass_total=r+s, mass_root=r, mass_shoot=s)%>%
  select(-r, -s)%>%
  merge(nodules, by=c('pot'), all=T)%>%
  merge(plants, by=c('pot'), all=T)


#merge harvest data with treatment information
harvestTrt <- merge(harvestAll, plan, by=c('pot'), all=T)%>%
  merge(trt, by=c('ring', 'plot'), all=T)%>%
  mutate(LECA_trt=as.factor(paste(Lespedeza.capitata, tot_legume, sep='::')), LUPE_trt=as.factor(paste(Lupinus.perennis, tot_legume, sep='::')), AMCA_trt=as.factor(paste(Amorpha.canescens, tot_legume, sep='::')))%>%
  #remove plants that never grew
  filter(num_plants>0)%>%
  #calcualte mass and nodules per plant (many pots had more than one plant in them)
  mutate(mass_total_ind=mass_total/num_plants, total_nod_ind=total_nod/num_plants, func_nod_ind=f_nodules/num_plants)%>%
  #remove outliers until they can be checked
  filter(pot!=200004060 & pot!=200003217 & pot!=200001213 & pot!=200003064 & pot!=200003072 & pot!=200003018)

harvestTrt$LECA_trt2 <- as.factor(with(harvestTrt, ifelse(LECA_trt=='0::1' & spp_count==1, 'Other Legume in Monoculture', ifelse(LECA_trt=='0::1' & spp_count==4, 'Other Legume in Polyculture', ifelse(LECA_trt=='1::1' & spp_count==1, 'Self in Monoculture', ifelse(LECA_trt=='1::1' & spp_count==4, 'Self in Polyculture', ifelse(LECA_trt=='1::4' & spp_count==4, 'All Legume Polyculture', ifelse(spp_count==16, '16 Species Polyculture', 'NA'))))))))

harvestTrt$LUPE_trt2 <- as.factor(with(harvestTrt, ifelse(LUPE_trt=='0::1' & spp_count==1, 'Other Legume in Monoculture', ifelse(LUPE_trt=='0::1' & spp_count==4, 'Other Legume in Polyculture', ifelse(LUPE_trt=='1::1' & spp_count==1, 'Self in Monoculture', ifelse(LUPE_trt=='1::1' & spp_count==4, 'Self in Polyculture', ifelse(LUPE_trt=='1::4' & spp_count==4, 'All Legume Polyculture', ifelse(spp_count==16, '16 Species Polyculture', 'NA'))))))))

harvestTrt$AMCA_trt2 <- as.factor(with(harvestTrt, ifelse(AMCA_trt=='0::1' & spp_count==1, 'Other Legume in Monoculture', ifelse(AMCA_trt=='0::1' & spp_count==4, 'Other Legume in Polyculture', ifelse(AMCA_trt=='1::1' & spp_count==1, 'Self in Monoculture', ifelse(AMCA_trt=='1::1' & spp_count==4, 'Self in Polyculture', ifelse(AMCA_trt=='1::4' & spp_count==4, 'All Legume Polyculture', ifelse(spp_count==16, '16 Species Polyculture', 'NA'))))))))

harvestTrt <- harvestTrt%>% select(-LUPE_trt, -LECA_trt, -AMCA_trt)



  



###calculate the relative biomass compared to uninoculated controls (percent difference)
#subset out the soil inoculation controls and inoculated plants
harvestCtl <- harvestTrt%>%
  filter(plot=='control')%>%
  #filter the few plants that had nodules (ring 4 AMCA=6; ring 1 AMCA=3; ring 5 LUPE=1)
  filter(total_nod==0)
names(harvestCtl)[names(harvestCtl)=='mass_shoot'] <- 'mass_shoot_ctl'
names(harvestCtl)[names(harvestCtl)=='mass_root'] <- 'mass_root_ctl'
names(harvestCtl)[names(harvestCtl)=='mass_total'] <- 'mass_total_ctl'
harvestCtlShoot <- harvestCtl%>% select(-mass_root_ctl, -mass_total_ctl)%>%
  group_by(ctl_N_trt, ctl_CO2_trt, spp)%>%
  summarise(mass_shoot_ctl=mean(mass_shoot_ctl))%>%
  ungroup()%>%
  mutate(CO2_trt=ctl_CO2_trt, N_trt=ctl_N_trt)%>%
  select(-ctl_CO2_trt, -ctl_N_trt)
harvestCtlRoot <- harvestCtl%>% select(-mass_shoot_ctl, -mass_total_ctl)%>%
  group_by(ctl_N_trt, ctl_CO2_trt, spp)%>%
  summarise(mass_root_ctl=mean(mass_root_ctl))%>%
  ungroup()%>%
  mutate(CO2_trt=ctl_CO2_trt, N_trt=ctl_N_trt)%>%
  select(-ctl_CO2_trt, -ctl_N_trt)  
harvestCtlTotBio <- harvestCtl%>% select(-mass_shoot_ctl, -mass_root_ctl)%>%
  group_by(ctl_N_trt, ctl_CO2_trt, spp)%>%
  summarise(mass_total_ctl=mean(mass_total_ctl))%>%
  ungroup()%>%
  mutate(CO2_trt=ctl_CO2_trt, N_trt=ctl_N_trt)%>%
  select(-ctl_CO2_trt, -ctl_N_trt)
harvestCtlAll <- merge(harvestCtlShoot, harvestCtlRoot, by=c('spp', 'CO2_trt', 'N_trt'))%>%
  merge(harvestCtlTotBio, by=c('spp', 'CO2_trt', 'N_trt'))

#subset out inoculated plants
harvestIno <- harvestTrt%>% filter(plot!='control')

#merge uninoculated controls back with inoculated plants, calculate relative cover, drop unneccesary columns
harvestRel <- merge(harvestIno, harvestCtlAll, by=c('CO2_trt', 'N_trt', 'spp'), all=T)%>%
  mutate(mass_shoot_rel=((mass_shoot-mass_shoot_ctl)/mass_shoot_ctl), mass_root_rel=((mass_root-mass_root_ctl)/mass_root_ctl), mass_total_rel=((mass_total-mass_total_ctl)/mass_total_ctl))%>%
  select(-mass_shoot_ctl, -mass_root_ctl, -mass_total_ctl)%>%
  #average across replicates
  group_by(CO2_trt, N_trt, spp, ring, plot, spp_count, group_count, experiment, monospecies, monogroup, Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, C.3, C.4, Forb, Legume, tot_legume, LECA_trt2, LUPE_trt2, AMCA_trt2)%>%
  summarise(mass_shoot=mean(mass_shoot), mass_root=mean(mass_root), mass_total=mean(mass_total), f_nodules=mean(f_nodules), n_nodules=mean(n_nodules), total_nod=mean(total_nod), mass_shoot_rel=mean(mass_shoot_rel), mass_root_rel=mean(mass_root_rel), mass_total_rel=mean(mass_total_rel))%>%
  mutate(trt=paste(CO2_trt, N_trt, sep='_'))


###soil feedbacks under control conditions
#subset out CO2 and N treatments
#separate models for each species

#shoot mass
summary(lmer(mass_shoot_rel ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_shoot_rel ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_shoot_rel ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#root mass
summary(lmer(mass_root_rel ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_root_rel ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_root_rel ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#total biomass
summary(lmer(mass_total_rel ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_total_rel ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_total_rel ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#functional nodules
summary(lmer(f_nodules ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(f_nodules ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(f_nodules ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#non-functional nodules
summary(lmer(n_nodules ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(n_nodules ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(n_nodules ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#total nodules
summary(lmer(total_nod ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(total_nod ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(total_nod ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#relativized biomass response (to uninoculated plants) and averaged across all pot reps
feedbackLECAbio <- ggplot(barGraphStats(data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich'), variable='mass_total_rel', byFactorNames=c('LECA_trt2')), aes(x=LECA_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('Relative biomass (g)') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,2)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=2, label='(a) LECA Biomass', size=9, hjust=0) #+
#   annotate('text', x=1, y=0.44, label='b', size=7) + 
#   annotate('text', x=2, y=0.75, label='b', size=7) +
#   annotate('text', x=3, y=0.07, label='a', size=7)

feedbackLUPEbio <- ggplot(barGraphStats(data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich'), variable='mass_total_rel', byFactorNames=c('LUPE_trt2')), aes(x=LUPE_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0, 1.5)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=1.5, label='(b) LUPE Biomass', size=9, hjust=0)

feedbackAMCAbio <- ggplot(barGraphStats(data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich'), variable='mass_total_rel', byFactorNames=c('AMCA_trt2')), aes(x=AMCA_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0, 1.5)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=1.5, label='(c) AMCA Biomass', size=9, hjust=0)

pushViewport(viewport(layout=grid.layout(1,3)))
print(feedbackLECAbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(feedbackLUPEbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(feedbackAMCAbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))



#nodule response (per individual plant in a pot), not averaged across replicate pots
feedbackLECAnod <- ggplot(barGraphStats(data=subset(harvestTrt, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich' & total_nod!='is.na' & ring!='is.na'), variable='func_nod_ind', byFactorNames=c('LECA_trt2')), aes(x=LECA_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Nodule Number') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,12)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=12, label='(a) LECA Nodules', size=9, hjust=0)

feedbackLUPEnod <- ggplot(barGraphStats(data=subset(harvestTrt, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich' & total_nod!='is.na' & ring!='is.na'), variable='func_nod_ind', byFactorNames=c('LUPE_trt2')), aes(x=LUPE_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,12)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=12, label='(b) LUPE Nodules', size=9, hjust=0)

feedbackAMCAnod <- ggplot(barGraphStats(data=subset(harvestTrt, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich' & ring!='is.na'), variable='func_nod_ind', byFactorNames=c('AMCA_trt2')), aes(x=AMCA_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,12)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=12, label='(c) AMCA Nodules', size=9, hjust=0)

pushViewport(viewport(layout=grid.layout(1,3)))
print(feedbackLECAnod, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(feedbackLUPEnod, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(feedbackAMCAnod, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))



###not relativized by control biomass and not summarized across reps
# #shoot mass
# summary(lmer(mass_shoot ~ LECA_trt2 + (1|num_plants), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))
# 
# summary(lmer(mass_shoot ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))
# 
# summary(lmer(mass_shoot ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))
# 
# #root mass
# summary(lmer(mass_root ~ LECA_trt2 + (1|ring), data=subset(harvestRel, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))
# 
# summary(lmer(mass_root ~ LUPE_trt2 + (1|ring), data=subset(harvestRel, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))
# 
# summary(lmer(mass_root ~ AMCA_trt2 + (1|ring), data=subset(harvestRel, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

#total biomass
summary(lmer(mass_total_ind ~ LECA_trt2 + (1|ring), data=subset(harvestTrt, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_total_ind ~ LUPE_trt2 + (1|ring), data=subset(harvestTrt, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

summary(lmer(mass_total_ind ~ AMCA_trt2 + (1|ring), data=subset(harvestTrt, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich')))

feedbackLECAbio <- ggplot(barGraphStats(data=subset(harvestTrt, spp=='LECA' & CO2_trt!='Cenrich' & N_trt!='Nenrich' & ring!='is.na'), variable='mass_total_ind', byFactorNames=c('LECA_trt2')), aes(x=LECA_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('Total biomass (g)') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,0.12)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=0.12, label='(a) LECA Biomass', size=9, hjust=0) #+
#   annotate('text', x=1, y=0.44, label='b', size=7) + 
#   annotate('text', x=2, y=0.75, label='b', size=7) +
#   annotate('text', x=3, y=0.07, label='a', size=7)

feedbackLUPEbio <- ggplot(barGraphStats(data=subset(harvestTrt, spp=='LUPE' & CO2_trt!='Cenrich' & N_trt!='Nenrich' & mass_total_ind!='is.na' & ring!='is.na'), variable='mass_total_ind', byFactorNames=c('LUPE_trt2')), aes(x=LUPE_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,0.3)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=0.3, label='(b) LUPE Biomass', size=9, hjust=0)

feedbackAMCAbio <- ggplot(barGraphStats(data=subset(harvestTrt, spp=='AMCA' & CO2_trt!='Cenrich' & N_trt!='Nenrich' & ring!='is.na'), variable='mass_total_ind', byFactorNames=c('AMCA_trt2')), aes(x=AMCA_trt2, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0, 0.07)) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  annotate('text', x=0.5, y=0.07, label='(c) AMCA Biomass', size=9, hjust=0)

pushViewport(viewport(layout=grid.layout(1,3)))
print(feedbackLECAbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(feedbackLUPEbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(feedbackAMCAbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))




###rhizobial specialization (ratio of each sp's nodulation in its own monoculture soil vs other spp's monoculture soils)
#needs bootstrapping for CI
rhizSpecializationOtherFunction <- function(d, i){
  d2 <- d[i,]
  return(d2$other_nods)
}

rhizSpecializationSelfFunction <- function(d, i){
  d2 <- d[i,]
  return(d2$self_nods)
}

#LECA
lecaNods <- subset(harvestTrt, subset=(spp=='LECA'))%>%
  filter(LECA_trt2=='Self in Polyculture'|LECA_trt2=='Other Legume in Polyculture')%>%
  mutate(LECA_trt3=ifelse(LECA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Namb'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, LECA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LECA_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

lecaSelfNods <- subset(lecaNods, subset=(LECA_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
lecaOtherNods <- subset(lecaNods, subset=(LECA_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

lecaNodsSelfBootModel <- boot(lecaSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(lecaNodsSelfBootModel)
lecaNodsSelfBoot <- as.data.frame(lecaNodsSelfBootModel$t)
lecaNodsSelfBootMean <- as.data.frame(rowMeans(lecaNodsSelfBoot[1:1000,]))
names(lecaNodsSelfBootMean)[names(lecaNodsSelfBootMean) == 'rowMeans(lecaNodsSelfBoot[1:1000, ])'] <- 'self_nods'

lecaNodsOtherBootModel <- boot(lecaOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(lecaNodsOtherBootModel)
lecaNodsOtherBoot <- as.data.frame(lecaNodsOtherBootModel$t)
lecaNodsOtherBootMean <- as.data.frame(rowMeans(lecaNodsOtherBoot[1:1000,]))
names(lecaNodsOtherBootMean)[names(lecaNodsOtherBootMean) == 'rowMeans(lecaNodsOtherBoot[1:1000, ])'] <- 'other_nods'

lecaSpecializationBoot <- cbind(lecaNodsSelfBootMean, lecaNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='LECA')

#LUPE
lupeNods <- subset(harvestTrt, subset=(spp=='LUPE'))%>%
  filter(LUPE_trt2=='Self in Polyculture'|LUPE_trt2=='Other Legume in Polyculture')%>%
  mutate(LUPE_trt3=ifelse(LUPE_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Namb'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, LUPE_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LUPE_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

lupeSelfNods <- subset(lupeNods, subset=(LUPE_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
lupeOtherNods <- subset(lupeNods, subset=(LUPE_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

lupeNodsSelfBootModel <- boot(lupeSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(lupeNodsSelfBootModel)
lupeNodsSelfBoot <- as.data.frame(lupeNodsSelfBootModel$t)
lupeNodsSelfBootMean <- as.data.frame(rowMeans(lupeNodsSelfBoot[1:1000,]))
names(lupeNodsSelfBootMean)[names(lupeNodsSelfBootMean) == 'rowMeans(lupeNodsSelfBoot[1:1000, ])'] <- 'self_nods'

lupeNodsOtherBootModel <- boot(lupeOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(lupeNodsOtherBootModel)
lupeNodsOtherBoot <- as.data.frame(lupeNodsOtherBootModel$t)
lupeNodsOtherBootMean <- as.data.frame(rowMeans(lupeNodsOtherBoot[1:1000,]))
names(lupeNodsOtherBootMean)[names(lupeNodsOtherBootMean) == 'rowMeans(lupeNodsOtherBoot[1:1000, ])'] <- 'other_nods'

lupeSpecializationBoot <- cbind(lupeNodsSelfBootMean, lupeNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='LUPE')

#AMCA
amcaNods <- subset(harvestTrt, subset=(spp=='AMCA'))%>%
  filter(AMCA_trt2=='Self in Polyculture'|AMCA_trt2=='Other Legume in Polyculture')%>%
  mutate(AMCA_trt3=ifelse(AMCA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Namb'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, AMCA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, AMCA_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

amcaSelfNods <- subset(amcaNods, subset=(AMCA_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
amcaOtherNods <- subset(amcaNods, subset=(AMCA_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

amcaNodsSelfBootModel <- boot(amcaSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(amcaNodsSelfBootModel)
amcaNodsSelfBoot <- as.data.frame(amcaNodsSelfBootModel$t)
amcaNodsSelfBootMean <- as.data.frame(rowMeans(amcaNodsSelfBoot[1:1000,]))
names(amcaNodsSelfBootMean)[names(amcaNodsSelfBootMean) == 'rowMeans(amcaNodsSelfBoot[1:1000, ])'] <- 'self_nods'

amcaNodsOtherBootModel <- boot(amcaOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(amcaNodsOtherBootModel)
amcaNodsOtherBoot <- as.data.frame(amcaNodsOtherBootModel$t)
amcaNodsOtherBootMean <- as.data.frame(rowMeans(amcaNodsOtherBoot[1:1000,]))
names(amcaNodsOtherBootMean)[names(amcaNodsOtherBootMean) == 'rowMeans(amcaNodsOtherBoot[1:1000, ])'] <- 'other_nods'

amcaSpecializationBoot <- cbind(amcaNodsSelfBootMean, amcaNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='AMCA')%>%
  #a few go to infinity because so few self pots
  filter(specialization!='Inf')

#combine spp specializations
allSpecializationBoot <- rbind(lecaSpecializationBoot, lupeSpecializationBoot, amcaSpecializationBoot)%>%
  group_by(spp)%>%
  summarize(specialization_mean=mean(specialization), specialization_sd=sd(specialization))%>%
  ungroup()%>%
  mutate(specialization_CI=1.96*specialization_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=specialization_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=specialization_mean-specialization_CI, ymax=specialization_mean+specialization_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Rhizobial Specialization') +
  annotate('text', x=1, y=0.48, label='a', size=10) +
  annotate('text', x=2, y=0.35, label='a', size=10) +
  annotate('text', x=3, y=1.35, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500




###biomass Relative Yield (ratio of each sp's biomass in its own polyculture soil vs other spp's polyculture soils)
#needs bootstrapping for CI
biomassOtherFunction <- function(d, i){
  d2 <- d[i,]
  return(d2$other_mass)
}

biomassSelfFunction <- function(d, i){
  d2 <- d[i,]
  return(d2$self_mass)
}

#LECA
lecaBio <- subset(harvestTrt, subset=(spp=='LECA'))%>%
  filter(LECA_trt2=='Self in Polyculture'|LECA_trt2=='All Legume Polyculture')%>%
  mutate(LECA_trt3=ifelse(LECA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Namb'&mass_total_ind!='is.na')%>%
  select(ring, plot, pot, mass_total_ind, LECA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LECA_trt3)%>%
  summarise(mass_total_ind2=mean(mass_total_ind))%>%
  ungroup()

lecaSelfMass <- subset(lecaBio, subset=(LECA_trt3=='self'))%>%
  mutate(self_mass=mass_total_ind2)%>%
  select(self_mass)
lecaOtherMass <- subset(lecaBio, subset=(LECA_trt3=='other'))%>%
  mutate(other_mass=mass_total_ind2)%>%
  select(other_mass)

lecaMassSelfBootModel <- boot(lecaSelfMass, biomassSelfFunction, R=1000)
plot(lecaMassSelfBootModel)
lecaMassSelfBoot <- as.data.frame(lecaMassSelfBootModel$t)
lecaMassSelfBootMean <- as.data.frame(rowMeans(lecaMassSelfBoot[1:1000,]))
names(lecaMassSelfBootMean)[names(lecaMassSelfBootMean) == 'rowMeans(lecaMassSelfBoot[1:1000, ])'] <- 'self_mass'

lecaMassOtherBootModel <- boot(lecaOtherMass, biomassOtherFunction, R=1000)
plot(lecaMassOtherBootModel)
lecaMassOtherBoot <- as.data.frame(lecaMassOtherBootModel$t)
lecaMassOtherBootMean <- as.data.frame(rowMeans(lecaMassOtherBoot[1:1000,]))
names(lecaMassOtherBootMean)[names(lecaMassOtherBootMean) == 'rowMeans(lecaMassOtherBoot[1:1000, ])'] <- 'other_mass'

lecaSpecializationBoot <- cbind(lecaMassSelfBootMean, lecaMassOtherBootMean)%>%
  mutate(specialization=other_mass/self_mass, spp='LECA')

#LUPE
lupeBio <- subset(harvestTrt, subset=(spp=='LUPE'))%>%
  filter(LUPE_trt2=='Self in Polyculture'|LUPE_trt2=='All Legume Polyculture')%>%
  mutate(LUPE_trt3=ifelse(LUPE_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Namb'&mass_total_ind!='is.na')%>%
  select(ring, plot, pot, mass_total_ind, LUPE_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LUPE_trt3)%>%
  summarise(mass_total_ind2=mean(mass_total_ind))%>%
  ungroup()

lupeSelfMass <- subset(lupeBio, subset=(LUPE_trt3=='self'))%>%
  mutate(self_mass=mass_total_ind2)%>%
  select(self_mass)
lupeOtherMass <- subset(lupeBio, subset=(LUPE_trt3=='other'))%>%
  mutate(other_mass=mass_total_ind2)%>%
  select(other_mass)

lupeMassSelfBootModel <- boot(lupeSelfMass, biomassSelfFunction, R=1000)
plot(lupeMassSelfBootModel)
lupeMassSelfBoot <- as.data.frame(lupeMassSelfBootModel$t)
lupeMassSelfBootMean <- as.data.frame(rowMeans(lupeMassSelfBoot[1:1000,]))
names(lupeMassSelfBootMean)[names(lupeMassSelfBootMean) == 'rowMeans(lupeMassSelfBoot[1:1000, ])'] <- 'self_mass'

lupeMassOtherBootModel <- boot(lupeOtherMass, biomassOtherFunction, R=1000)
plot(lupeMassOtherBootModel)
lupeMassOtherBoot <- as.data.frame(lupeMassOtherBootModel$t)
lupeMassOtherBootMean <- as.data.frame(rowMeans(lupeMassOtherBoot[1:1000,]))
names(lupeMassOtherBootMean)[names(lupeMassOtherBootMean) == 'rowMeans(lupeMassOtherBoot[1:1000, ])'] <- 'other_mass'

lupeSpecializationBoot <- cbind(lupeMassSelfBootMean, lupeMassOtherBootMean)%>%
  mutate(specialization=other_mass/self_mass, spp='LUPE')

#AMCA
amcaBio <- subset(harvestTrt, subset=(spp=='AMCA'))%>%
  filter(AMCA_trt2=='Self in Polyculture'|AMCA_trt2=='All Legume Polyculture')%>%
  mutate(AMCA_trt3=ifelse(AMCA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Namb'&mass_total_ind!='is.na')%>%
  select(ring, plot, pot, mass_total_ind, AMCA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, AMCA_trt3)%>%
  summarise(mass_total_ind2=mean(mass_total_ind))%>%
  ungroup()

amcaSelfMass <- subset(amcaBio, subset=(AMCA_trt3=='self'))%>%
  mutate(self_mass=mass_total_ind2)%>%
  select(self_mass)
amcaOtherMass <- subset(amcaBio, subset=(AMCA_trt3=='other'))%>%
  mutate(other_mass=mass_total_ind2)%>%
  select(other_mass)

amcaMassSelfBootModel <- boot(amcaSelfMass, biomassSelfFunction, R=1000)
plot(amcaMassSelfBootModel)
amcaMassSelfBoot <- as.data.frame(amcaMassSelfBootModel$t)
amcaMassSelfBootMean <- as.data.frame(rowMeans(amcaMassSelfBoot[1:1000,]))
names(amcaMassSelfBootMean)[names(amcaMassSelfBootMean) == 'rowMeans(amcaMassSelfBoot[1:1000, ])'] <- 'self_mass'

amcaMassOtherBootModel <- boot(amcaOtherMass, biomassOtherFunction, R=1000)
plot(amcaMassOtherBootModel)
amcaMassOtherBoot <- as.data.frame(amcaMassOtherBootModel$t)
amcaMassOtherBootMean <- as.data.frame(rowMeans(amcaMassOtherBoot[1:1000,]))
names(amcaMassOtherBootMean)[names(amcaMassOtherBootMean) == 'rowMeans(amcaMassOtherBoot[1:1000, ])'] <- 'other_mass'

amcaSpecializationBoot <- cbind(amcaMassSelfBootMean, amcaMassOtherBootMean)%>%
  mutate(specialization=other_mass/self_mass, spp='AMCA')

#combine spp specializations
allSpecializationBoot <- rbind(lecaSpecializationBoot, lupeSpecializationBoot, amcaSpecializationBoot)%>%
  group_by(spp)%>%
  summarize(specialization_mean=mean(specialization), specialization_sd=sd(specialization))%>%
  ungroup()%>%
  mutate(specialization_CI=1.96*specialization_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=specialization_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=specialization_mean-specialization_CI, ymax=specialization_mean+specialization_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Relative Yield in Pots') +
  # annotate('text', x=1, y=0.5, label='a', size=10) +
  # annotate('text', x=2, y=0.4, label='a', size=10) +
  # annotate('text', x=3, y=1.3, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500








###quick ANOVA for CO2 and N trts
summary(response <- lmer(mass_total_rel ~ CO2_trt + N_trt + CO2_trt*N_trt + spp + spp*CO2_trt + spp*N_trt + spp*CO2_trt*N_trt + (1|ring), data=harvestRel))

feedbackLECAbio <- ggplot(barGraphStats(data=subset(harvestRel, spp=='LECA'), variable='mass_total_rel', byFactorNames=c('LECA_trt2', 'CO2_trt', 'N_trt')), aes(x=LECA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('Relative biomass (g)') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(-0.25,2.0)) +
  annotate('text', x=0.5, y=2.0, label='(a) LECA Biomass', size=9, hjust=0)  +
  theme(axis.text.x=element_text(angle=45, hjust=1)) #+
#   annotate('text', x=1, y=0.44, label='b', size=7) + 
#   annotate('text', x=2, y=0.75, label='b', size=7) +
#   annotate('text', x=3, y=0.07, label='a', size=7)

feedbackLUPEbio <- ggplot(barGraphStats(data=subset(harvestRel, spp=='LUPE'), variable='mass_total_rel', byFactorNames=c('LUPE_trt2', 'CO2_trt', 'N_trt')), aes(x=LUPE_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(-0.25, 1.5)) +
  annotate('text', x=0.5, y=1.5, label='(b) LUPE Biomass', size=9, hjust=0) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

feedbackAMCAbio <- ggplot(barGraphStats(data=subset(harvestRel, spp=='AMCA'), variable='mass_total_rel', byFactorNames=c('AMCA_trt2', 'CO2_trt', 'N_trt')), aes(x=AMCA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(-0.25, 3.2)) +
  annotate('text', x=0.5, y=3.2, label='(c) AMCA Biomass', size=9, hjust=0) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


feedbackLECAnod <- ggplot(barGraphStats(data=subset(harvestRel, spp=='LECA'), variable='total_nod', byFactorNames=c('LECA_trt2', 'CO2_trt', 'N_trt')), aes(x=LECA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Nodule Number') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,15)) +
  annotate('text', x=0.5, y=15, label='(c) LECA Nodules', size=9, hjust=0) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

feedbackLUPEnod <- ggplot(barGraphStats(data=subset(harvestRel, spp=='LUPE'), variable='total_nod', byFactorNames=c('LUPE_trt2', 'CO2_trt', 'N_trt')), aes(x=LUPE_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,12)) +
  annotate('text', x=0.5, y=12, label='(d) LUPE Nodules', size=9, hjust=0) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

feedbackAMCAnod <- ggplot(barGraphStats(data=subset(harvestRel, spp=='AMCA'), variable='total_nod', byFactorNames=c('AMCA_trt2', 'CO2_trt', 'N_trt')), aes(x=AMCA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('') +
  scale_x_discrete(limits=c('Other Legume in Monoculture', 'Other Legume in Polyculture', 'Self in Monoculture', 'Self in Polyculture', 'All Legume Polyculture', '16 Species Polyculture')) +
  coord_cartesian(ylim=c(0,22)) +
  annotate('text', x=0.5, y=22, label='(e) AMCA Nodules', size=9, hjust=0) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

pushViewport(viewport(layout=grid.layout(2,3)))
print(feedbackLECAbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(feedbackLUPEbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(feedbackAMCAbio, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(feedbackLECAnod, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(feedbackLUPEnod, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(feedbackAMCAnod, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))


#nodulation with broad bins of self vs all
ggplot(barGraphStats(data=subset(harvestTrt, spp=='LECA' & LECA_trt2!='Other Legume in Monoculture' & LECA_trt2!='Self in Monoculture' & LECA_trt2!='16 Species Polyculture' & total_nod_ind!='NA'), variable='total_nod_ind', byFactorNames=c('LECA_trt2', 'CO2_trt', 'N_trt')), aes(x=LECA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Total Nodules') +
  scale_x_discrete(labels=c('Other Legume in Polyculture'='Other Legume', 'Self in Polyculture'='Grown Alone', 'All Legume Polyculture'='All Legumes'), limits=c('Other Legume in Polyculture', 'Self in Polyculture', 'All Legume Polyculture')) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  coord_cartesian(ylim=c(0,12)) +
  annotate('text', x=0.5, y=12, label='(a) LECA Nodules', size=9, hjust=0)

ggplot(barGraphStats(data=subset(harvestTrt, spp=='LUPE' & LUPE_trt2!='Other Legume in Monoculture' & LUPE_trt2!='Self in Monoculture' & LUPE_trt2!='16 Species Polyculture'), variable='total_nod_ind', byFactorNames=c('LUPE_trt2', 'CO2_trt', 'N_trt')), aes(x=LUPE_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Total Nodules') +
  scale_x_discrete(labels=c('Other Legume in Polyculture'='Other Legume', 'Self in Polyculture'='Grown Alone', 'All Legume Polyculture'='All Legumes'), limits=c('Other Legume in Polyculture', 'Self in Polyculture', 'All Legume Polyculture')) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  coord_cartesian(ylim=c(0,12)) +
  annotate('text', x=0.5, y=12, label='(b) LUPE Nodules', size=9, hjust=0)

ggplot(barGraphStats(data=subset(harvestTrt, spp=='AMCA' & AMCA_trt2!='Other Legume in Monoculture' & AMCA_trt2!='Self in Monoculture' & AMCA_trt2!='16 Species Polyculture' & total_nod_ind!='NA'), variable='total_nod_ind', byFactorNames=c('AMCA_trt2', 'CO2_trt', 'N_trt')), aes(x=AMCA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Total Nodules') +
  scale_x_discrete(labels=c('Other Legume in Polyculture'='Other Legume', 'Self in Polyculture'='Grown Alone', 'All Legume Polyculture'='All Legumes'), limits=c('Other Legume in Polyculture', 'Self in Polyculture', 'All Legume Polyculture')) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  coord_cartesian(ylim=c(0,12)) +
  annotate('text', x=0.5, y=12, label='(c) AMCA Nodules', size=9, hjust=0)


#biomass with broad bins of self vs all
ggplot(barGraphStats(data=subset(harvestTrt, spp=='LECA' & LECA_trt2!='Other Legume in Monoculture' & LECA_trt2!='Self in Monoculture' & LECA_trt2!='16 Species Polyculture' & mass_total_ind!='NA'), variable='mass_total_ind', byFactorNames=c('LECA_trt2', 'CO2_trt', 'N_trt')), aes(x=LECA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Total Biomass') +
  scale_x_discrete(labels=c('Other Legume in Polyculture'='Other Legume', 'Self in Polyculture'='Grown Alone', 'All Legume Polyculture'='All Legumes'), limits=c('Other Legume in Polyculture', 'Self in Polyculture', 'All Legume Polyculture')) +
  scale_y_continuous(breaks=seq(0,0.2,0.05)) +
  coord_cartesian(ylim=c(0,0.2)) +
  annotate('text', x=0.5, y=0.2, label='(a) LECA Biomass', size=9, hjust=0)

ggplot(barGraphStats(data=subset(harvestTrt, spp=='LUPE' & LUPE_trt2!='Other Legume in Monoculture' & LUPE_trt2!='Self in Monoculture' & LUPE_trt2!='16 Species Polyculture' & mass_total_ind!='NA'), variable='mass_total_ind', byFactorNames=c('LUPE_trt2', 'CO2_trt', 'N_trt')), aes(x=LUPE_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Total Biomass') +
  scale_x_discrete(labels=c('Other Legume in Polyculture'='Other Legume', 'Self in Polyculture'='Grown Alone', 'All Legume Polyculture'='All Legumes'), limits=c('Other Legume in Polyculture', 'Self in Polyculture', 'All Legume Polyculture')) +
  scale_y_continuous(breaks=seq(0,0.4,0.1)) +
  coord_cartesian(ylim=c(0,0.4)) +
  annotate('text', x=0.5, y=0.4, label='(b) LUPE Biomass', size=9, hjust=0)

ggplot(barGraphStats(data=subset(harvestTrt, spp=='AMCA' & AMCA_trt2!='Other Legume in Monoculture' & AMCA_trt2!='Self in Monoculture' & AMCA_trt2!='16 Species Polyculture'), variable='mass_total', byFactorNames=c('AMCA_trt2', 'CO2_trt', 'N_trt')), aes(x=AMCA_trt2, y=mean, fill=interaction(CO2_trt, N_trt))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(.9), width=0.2) +
  xlab('Soil Origin') + ylab('Total Biomass') +
  scale_x_discrete(labels=c('Other Legume in Polyculture'='Other Legume', 'Self in Polyculture'='Grown Alone', 'All Legume Polyculture'='All Legumes'), limits=c('Other Legume in Polyculture', 'Self in Polyculture', 'All Legume Polyculture')) +
  scale_y_continuous(breaks=seq(0,0.15,0.05)) +
  coord_cartesian(ylim=c(0,0.15)) +
  annotate('text', x=0.5, y=0.15, label='(c) AMCA Biomass', size=9, hjust=0)




##eCO2 specialization
#LECA
lecaNods <- subset(harvestTrt, subset=(spp=='LECA'))%>%
  filter(LECA_trt2=='Self in Polyculture'|LECA_trt2=='Other Legume in Polyculture')%>%
  mutate(LECA_trt3=ifelse(LECA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Cenrich'&N_trt=='Namb'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, LECA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LECA_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

lecaSelfNods <- subset(lecaNods, subset=(LECA_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
lecaOtherNods <- subset(lecaNods, subset=(LECA_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

lecaNodsSelfBootModel <- boot(lecaSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(lecaNodsSelfBootModel)
lecaNodsSelfBoot <- as.data.frame(lecaNodsSelfBootModel$t)
lecaNodsSelfBootMean <- as.data.frame(rowMeans(lecaNodsSelfBoot[1:1000,]))
names(lecaNodsSelfBootMean)[names(lecaNodsSelfBootMean) == 'rowMeans(lecaNodsSelfBoot[1:1000, ])'] <- 'self_nods'

lecaNodsOtherBootModel <- boot(lecaOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(lecaNodsOtherBootModel)
lecaNodsOtherBoot <- as.data.frame(lecaNodsOtherBootModel$t)
lecaNodsOtherBootMean <- as.data.frame(rowMeans(lecaNodsOtherBoot[1:1000,]))
names(lecaNodsOtherBootMean)[names(lecaNodsOtherBootMean) == 'rowMeans(lecaNodsOtherBoot[1:1000, ])'] <- 'other_nods'

lecaSpecializationBoot <- cbind(lecaNodsSelfBootMean, lecaNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='LECA')

#LUPE
lupeNods <- subset(harvestTrt, subset=(spp=='LUPE'))%>%
  filter(LUPE_trt2=='Self in Polyculture'|LUPE_trt2=='Other Legume in Polyculture')%>%
  mutate(LUPE_trt3=ifelse(LUPE_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Cenrich'&N_trt=='Namb'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, LUPE_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LUPE_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

lupeSelfNods <- subset(lupeNods, subset=(LUPE_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
lupeOtherNods <- subset(lupeNods, subset=(LUPE_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

lupeNodsSelfBootModel <- boot(lupeSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(lupeNodsSelfBootModel)
lupeNodsSelfBoot <- as.data.frame(lupeNodsSelfBootModel$t)
lupeNodsSelfBootMean <- as.data.frame(rowMeans(lupeNodsSelfBoot[1:1000,]))
names(lupeNodsSelfBootMean)[names(lupeNodsSelfBootMean) == 'rowMeans(lupeNodsSelfBoot[1:1000, ])'] <- 'self_nods'

lupeNodsOtherBootModel <- boot(lupeOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(lupeNodsOtherBootModel)
lupeNodsOtherBoot <- as.data.frame(lupeNodsOtherBootModel$t)
lupeNodsOtherBootMean <- as.data.frame(rowMeans(lupeNodsOtherBoot[1:1000,]))
names(lupeNodsOtherBootMean)[names(lupeNodsOtherBootMean) == 'rowMeans(lupeNodsOtherBoot[1:1000, ])'] <- 'other_nods'

lupeSpecializationBoot <- cbind(lupeNodsSelfBootMean, lupeNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='LUPE')

#AMCA
amcaNods <- subset(harvestTrt, subset=(spp=='AMCA'))%>%
  filter(AMCA_trt2=='Self in Polyculture'|AMCA_trt2=='Other Legume in Polyculture')%>%
  mutate(AMCA_trt3=ifelse(AMCA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Cenrich'&N_trt=='Namb'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, AMCA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, AMCA_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

amcaSelfNods <- subset(amcaNods, subset=(AMCA_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
amcaOtherNods <- subset(amcaNods, subset=(AMCA_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

amcaNodsSelfBootModel <- boot(amcaSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(amcaNodsSelfBootModel)
amcaNodsSelfBoot <- as.data.frame(amcaNodsSelfBootModel$t)
amcaNodsSelfBootMean <- as.data.frame(rowMeans(amcaNodsSelfBoot[1:1000,]))
names(amcaNodsSelfBootMean)[names(amcaNodsSelfBootMean) == 'rowMeans(amcaNodsSelfBoot[1:1000, ])'] <- 'self_nods'

amcaNodsOtherBootModel <- boot(amcaOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(amcaNodsOtherBootModel)
amcaNodsOtherBoot <- as.data.frame(amcaNodsOtherBootModel$t)
amcaNodsOtherBootMean <- as.data.frame(rowMeans(amcaNodsOtherBoot[1:1000,]))
names(amcaNodsOtherBootMean)[names(amcaNodsOtherBootMean) == 'rowMeans(amcaNodsOtherBoot[1:1000, ])'] <- 'other_nods'

amcaSpecializationBoot <- cbind(amcaNodsSelfBootMean, amcaNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='AMCA')%>%
  #a few go to infinity because so few self pots
  filter(specialization!='Inf')

#combine spp specializations
allSpecializationBoot <- rbind(lecaSpecializationBoot, lupeSpecializationBoot, amcaSpecializationBoot)%>%
  group_by(spp)%>%
  summarize(specialization_mean=mean(specialization), specialization_sd=sd(specialization))%>%
  ungroup()%>%
  mutate(specialization_CI=1.96*specialization_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=specialization_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=specialization_mean-specialization_CI, ymax=specialization_mean+specialization_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Rhizobial Specialization') +
  annotate('text', x=1, y=0.48, label='a', size=10) +
  annotate('text', x=2, y=0.35, label='a', size=10) +
  annotate('text', x=3, y=1.35, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500


##eN specialization
#LECA
lecaNods <- subset(harvestTrt, subset=(spp=='LECA'))%>%
  filter(LECA_trt2=='Self in Polyculture'|LECA_trt2=='Other Legume in Polyculture')%>%
  mutate(LECA_trt3=ifelse(LECA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Nenrich'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, LECA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LECA_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

lecaSelfNods <- subset(lecaNods, subset=(LECA_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
lecaOtherNods <- subset(lecaNods, subset=(LECA_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

lecaNodsSelfBootModel <- boot(lecaSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(lecaNodsSelfBootModel)
lecaNodsSelfBoot <- as.data.frame(lecaNodsSelfBootModel$t)
lecaNodsSelfBootMean <- as.data.frame(rowMeans(lecaNodsSelfBoot[1:1000,]))
names(lecaNodsSelfBootMean)[names(lecaNodsSelfBootMean) == 'rowMeans(lecaNodsSelfBoot[1:1000, ])'] <- 'self_nods'

lecaNodsOtherBootModel <- boot(lecaOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(lecaNodsOtherBootModel)
lecaNodsOtherBoot <- as.data.frame(lecaNodsOtherBootModel$t)
lecaNodsOtherBootMean <- as.data.frame(rowMeans(lecaNodsOtherBoot[1:1000,]))
names(lecaNodsOtherBootMean)[names(lecaNodsOtherBootMean) == 'rowMeans(lecaNodsOtherBoot[1:1000, ])'] <- 'other_nods'

lecaSpecializationBoot <- cbind(lecaNodsSelfBootMean, lecaNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='LECA')

#LUPE
lupeNods <- subset(harvestTrt, subset=(spp=='LUPE'))%>%
  filter(LUPE_trt2=='Self in Polyculture'|LUPE_trt2=='Other Legume in Polyculture')%>%
  mutate(LUPE_trt3=ifelse(LUPE_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Nenrich'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, LUPE_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, LUPE_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

lupeSelfNods <- subset(lupeNods, subset=(LUPE_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
lupeOtherNods <- subset(lupeNods, subset=(LUPE_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

lupeNodsSelfBootModel <- boot(lupeSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(lupeNodsSelfBootModel)
lupeNodsSelfBoot <- as.data.frame(lupeNodsSelfBootModel$t)
lupeNodsSelfBootMean <- as.data.frame(rowMeans(lupeNodsSelfBoot[1:1000,]))
names(lupeNodsSelfBootMean)[names(lupeNodsSelfBootMean) == 'rowMeans(lupeNodsSelfBoot[1:1000, ])'] <- 'self_nods'

lupeNodsOtherBootModel <- boot(lupeOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(lupeNodsOtherBootModel)
lupeNodsOtherBoot <- as.data.frame(lupeNodsOtherBootModel$t)
lupeNodsOtherBootMean <- as.data.frame(rowMeans(lupeNodsOtherBoot[1:1000,]))
names(lupeNodsOtherBootMean)[names(lupeNodsOtherBootMean) == 'rowMeans(lupeNodsOtherBoot[1:1000, ])'] <- 'other_nods'

lupeSpecializationBoot <- cbind(lupeNodsSelfBootMean, lupeNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='LUPE')

#AMCA
amcaNods <- subset(harvestTrt, subset=(spp=='AMCA'))%>%
  filter(AMCA_trt2=='Self in Polyculture'|AMCA_trt2=='Other Legume in Polyculture')%>%
  mutate(AMCA_trt3=ifelse(AMCA_trt2=='Self in Polyculture', 'self', 'other'))%>%
  filter(CO2_trt=='Camb'&N_trt=='Nenrich'&func_nod_ind!='is.na')%>%
  select(ring, plot, pot, func_nod_ind, AMCA_trt3)%>%
  #remove duplicate row - FIX THIS
  group_by(ring, plot, pot, AMCA_trt3)%>%
  summarise(total_nod_ind2=mean(func_nod_ind))%>%
  ungroup()

amcaSelfNods <- subset(amcaNods, subset=(AMCA_trt3=='self'))%>%
  mutate(self_nods=total_nod_ind2)%>%
  select(self_nods)
amcaOtherNods <- subset(amcaNods, subset=(AMCA_trt3=='other'))%>%
  mutate(other_nods=total_nod_ind2)%>%
  select(other_nods)

amcaNodsSelfBootModel <- boot(amcaSelfNods, rhizSpecializationSelfFunction, R=1000)
plot(amcaNodsSelfBootModel)
amcaNodsSelfBoot <- as.data.frame(amcaNodsSelfBootModel$t)
amcaNodsSelfBootMean <- as.data.frame(rowMeans(amcaNodsSelfBoot[1:1000,]))
names(amcaNodsSelfBootMean)[names(amcaNodsSelfBootMean) == 'rowMeans(amcaNodsSelfBoot[1:1000, ])'] <- 'self_nods'

amcaNodsOtherBootModel <- boot(amcaOtherNods, rhizSpecializationOtherFunction, R=1000)
plot(amcaNodsOtherBootModel)
amcaNodsOtherBoot <- as.data.frame(amcaNodsOtherBootModel$t)
amcaNodsOtherBootMean <- as.data.frame(rowMeans(amcaNodsOtherBoot[1:1000,]))
names(amcaNodsOtherBootMean)[names(amcaNodsOtherBootMean) == 'rowMeans(amcaNodsOtherBoot[1:1000, ])'] <- 'other_nods'

amcaSpecializationBoot <- cbind(amcaNodsSelfBootMean, amcaNodsOtherBootMean)%>%
  mutate(specialization=other_nods/self_nods, spp='AMCA')%>%
  #a few go to infinity because so few self pots
  filter(specialization!='Inf')

#combine spp specializations
allSpecializationBoot <- rbind(lecaSpecializationBoot, lupeSpecializationBoot, amcaSpecializationBoot)%>%
  group_by(spp)%>%
  summarize(specialization_mean=mean(specialization), specialization_sd=sd(specialization))%>%
  ungroup()%>%
  mutate(specialization_CI=1.96*specialization_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=specialization_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=specialization_mean-specialization_CI, ymax=specialization_mean+specialization_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Rhizobial Specialization') +
  annotate('text', x=1, y=0.48, label='a', size=10) +
  annotate('text', x=2, y=0.35, label='a', size=10) +
  annotate('text', x=3, y=1.35, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500