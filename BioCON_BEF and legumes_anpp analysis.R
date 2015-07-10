library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(nlme) 
library(lsmeans)
library(lavaan)

source('C:\\Users\\Kim\\Desktop\\general R code\\general-R-code\\bar graph summary stats.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\general R code\\general-R-code\\ggplot_theme set.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_treatment data.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_climate data.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_anpp data.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_cover data.r', chdir=T)

setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

#merge ANPP and climate data
anppClimate <- merge(anppRY, climate)
polyClimate <- merge(anppPolyRY, climate)
legume4to1Climate <- merge(legume4to1complete, climate)



##############################################
##############################################

#Diversity analyses only (i.e., control plots, no N or CO2 additions)

##############################################
##############################################

#mixed model
#IMPORTANT: 4 spp polycultures
#IMPORTANT: using growing year climate variables
RYTmixedCategorical <- lme(logRYT ~ gy_precip_cm*leg_num_spp, random=~1|plot, data=subset(anppClimate, spp_count==4 & trt=='Camb_Namb'))
summary(RYTmixedCategorical)
anova(RYTmixedCategorical)
lsmeans(RYTmixedCategorical, cld~leg_num_spp)

anppClimate$order <- factor(anppClimate$leg_num_spp, levels=c('0 legumes', 'AMCA', 'LECA', 'LUPE', 'PEVI', '2 legumes', '3 legumes', '4 legumes'))

ggplot(subset(anppClimate, spp_count==4 & trt=='Camb_Namb'), aes(x=gy_precip_cm, y=logRYT, colour=order)) +
  geom_point() +
  geom_smooth(method=lm) +
  ylab('log Relative Yield Total') +
  xlab('Total Annual Precipitation (cm)')


#ANOVA model
#IMPORTANT: 4 spp polycultures
legsppnumANOVA <- aov(logRYT ~ leg_num_spp, data=subset(anppClimate, spp_count==4 & trt=='Camb_Namb'))
summary(legsppnumANOVA)
lsmeans(legsppnumANOVA, cld~leg_num_spp)

ggplot(barGraphStats(data=subset(anppClimate, spp_count==4 & trt=='Camb_Namb'), variable='logRYT', byFactorNames=c('order')), aes(x=order, y=mean)) +
  geom_bar(stat='identity', fill='white', colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  ylab('log Relative Yield Total') +
  theme(axis.title.x=element_blank()) +
  annotate('text', x=c(1,2,5,8), y=c(0.14,0.14,0.08,0.17), label='a', size=8) +
  annotate('text', x=c(3,4,6), y=c(0.33,0.36,0.39), label='b', size=8) +
  annotate('text', x=7, y=0.3, label='ab', size=8) +
  coord_cartesian(ylim=c(0,0.45))


#mixed model
#IMPORTANT: 4 spp polycultures, only included plots with one of the four legumes (intra-functional group differences)
RYTmixedCategoricalLegs <- lme(logRYT ~ gy_precip_cm*leg_num_spp, random=~1|plot, data=subset(anppClimate, spp_count==4 & trt=='Camb_Namb' & leg_num_spp!='mix' & leg_num_spp!='none' & leg_num_spp!='all 4'))
summary(RYTmixedCategoricalLegs)
anova(RYTmixedCategoricalLegs)
lsmeans(RYTmixedCategoricalLegs, cld~leg_num_spp)

color <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")

ggplot(subset(anppClimate, spp_count==4 & trt=='Camb_Namb' & legume_num==1), aes(x=gy_precip_cm, y=logRYT, colour=order)) +
  geom_point(size=2) +
  geom_smooth(method=lm, size=1) +
  ylab('log Relative Yield Total') +
  xlab('Total Annual Precipitation (cm)') +
  scale_colour_manual(values=color)


#mixed model
#effect of other legumes on RY of each legume species
RYmodel <- lme(logRY ~ gy_precip_cm*other_legs, random=~1|plot, data=subset(polyClimate, spp_count==4 & trt=='Camb_Namb' & legume_num!=0 & legume_num!=2 & legume_num!=3))
summary(RYmodel)
anova(RYmodel)

polyClimate$order <- factor(polyClimate$other_legs, levels=c('0 legumes', 'single legume', '2 legumes', '3 legumes', '4 legumes'))
ggplot(barGraphStats(data=subset(polyClimate, spp_count==4 & trt=='Camb_Namb' & leg_num_spp!='2 legumes' & leg_num_spp!='3 legumes' & Legume==1 & spp_type=='legume'), variable='logRY', byFactorNames=c('species', 'order', 'gy_precip_cm')), aes(x=species, y=mean, fill=order)) +
  geom_bar(stat='identity', position=position_dodge(), colour='black') +
  scale_fill_manual(values=c("white", "grey")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('log Relative Yield') +
  facet_wrap(~gy_precip_cm, ncol=4) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())


#difference between 1 and 4 legume plots
RYmodelPrecip <- lm(diff ~ gy_precip_cm*species, data=subset(legume4to1Climate, trt=='Camb_Namb'))
summary(RYmodelPrecip)
anova(RYmodelPrecip)

color <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")

ggplot(subset(legume4to1Climate, trt=='Camb_Namb'), aes(x=gy_precip_cm, y=diff, colour=species)) +
  geom_point(size=2) +
  geom_smooth(method=lm, size=1, se=F) +
  xlab('Total Annual Precipitation (cm)') +
  ylab('Difference in Relative Yield') +
  scale_colour_manual(values=color,
                      labels=c('Amorpha canescens', 'Lespedeza capitata', 'Lupinus perennis', 'Petalostemum villosum'))

