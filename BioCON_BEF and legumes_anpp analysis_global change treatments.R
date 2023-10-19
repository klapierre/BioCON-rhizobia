library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(nlme) 
library(lsmeans)
library(lavaan)
library(tidyverse)

source('C:\\Users\\kjkomatsu\\Desktop\\general R code\\general-R-code\\bar graph summary stats.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\general R code\\general-R-code\\ggplot_theme set.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\BioCON_treatment data.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\BioCON_climate data.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\BioCON_anpp data.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\BioCON_cover data.r', chdir=T)

setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

#merge ANPP and climate data
anppClimate <- merge(anppRY, climate)
polyClimate <- merge(anppPolyRY, climate)
RRclimate <- merge(anppRR, climate)

#mixed model
#do legumes or non-legumes respond more to CO2 and N?
#this is just species RY in 4 spp polyculture, not considering what type of polyculture they were in (BIC and AIC were lower when year and/or gy_precip were left out of model)
trtType4sppMixed <- lme(RYscaled_log ~ CO2_trt*N_trt*spp_type, random=~1|plot, data=subset(polyClimate, spp_count==4))
summary(trtType4sppMixed)
anova(trtType4sppMixed)

color <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")
ggplot(barGraphStats(data=subset(polyClimate, spp_count==4), variable='RYscaled_log', byFactorNames=c('spp_type', 'trt')), aes(x=spp_type, y=mean, fill=trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=color) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('log Relative Yield') +
  theme(axis.title.x=element_blank())

#does hedgesd vary with gy precip?
summary(lm(hedgesd_CO2 ~ gy_precip_cm*spp_type, data=RRclimate)) #yes (marginal)
summary(lm(hedgesd_N ~ gy_precip_cm*spp_type, data=RRclimate)) #yes
summary(lm(hedgesd_CO2_N ~ gy_precip_cm*spp_type, data=RRclimate)) #yes

#mixed model with hedges d
anppHedgesDa <- RRclimate[,-c(5:16,20:25)]
anppHedgesD <- melt(anppHedgesDa, id.vars=c('year', 'gy_precip_cm', 'gy_temp', 'spp_type'), variable.name='trt', value.name='hedges_d')
hedgesdModel <- lm(hedges_d ~ year*trt*spp_type , data=anppHedgesD)
summary(hedgesdModel)
anova(hedgesdModel)

CO2plot <- ggplot(anppRR, aes(x=year, y=hedgesd_CO2, shape=spp_type)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, aes(linetype=spp_type), colour='black') +
  geom_errorbar(aes(ymin=hedgesd_CO2-hedgesd_ci_CO2, ymax=hedgesd_CO2+hedgesd_ci_CO2, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position=c(0.8,0.95), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-2.5,12), xlim=c(1997.5,2005.5)) +
  scale_y_continuous(breaks=seq(-2,12,2)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  annotate('text', x=1997.75, y=11.5, label='(a) CO2', size=7, hjust=0)
Nplot <- ggplot(anppRR, aes(x=year, y=hedgesd_N, shape=spp_type)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, aes(linetype=spp_type), colour='black') +
  geom_errorbar(aes(ymin=hedgesd_N-hedgesd_ci_N, ymax=hedgesd_N+hedgesd_ci_N, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-2.5,12), xlim=c(1997.5,2005.5)) +
  scale_y_continuous(breaks=seq(-2,12,2)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  annotate('text', x=1997.75, y=11.5, label='(b) N', size=7, hjust=0)
CO2Nplot <- ggplot(anppRR, aes(x=year, y=hedgesd_CO2_N, shape=spp_type)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, aes(linetype=spp_type), colour='black') +
  geom_errorbar(aes(ymin=hedgesd_CO2_N-hedgesd_ci_CO2_N, ymax=hedgesd_CO2_N+hedgesd_ci_CO2_N, width=0.2)) +
  ylab('Hedges d') +
  coord_cartesian(ylim=c(-2.5,12), xlim=c(1997.5,2005.5)) +
  scale_y_continuous(breaks=seq(-2,12,2)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45))+
  annotate('text', x=1997.75, y=11.5, label='(c) CO2 + N', size=7, hjust=0)

#put hedges d panels together
pushViewport(viewport(layout=grid.layout(1,3))) 
print(CO2plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(Nplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(CO2Nplot, vp=viewport(layout.pos.row=1, layout.pos.col=3))



#mixed model with hedges d
hedgesdModelSpp <- lm(hedges_d ~ year*trt*species , data=anppHedgesDSpp)
summary(hedgesdModelSpp)
anova(hedgesdModelSpp)

amcaplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Amorpha.canescens'), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position=c(0.8,0.85), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  scale_colour_discrete(labels=c('eCO2', 'eN', 'eCO2 + eN')) +
  scale_shape_discrete(labels=c('eCO2', 'eN', 'eCO2 + eN')) +
  annotate('text', x=1997.75, y=18, label='(a) A. canescens', size=7, hjust=0)
lecaplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Lespedeza.capitata'), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  annotate('text', x=1997.75, y=18, label='(b) Le. capitata', size=7, hjust=0)
lupeplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Lupinus.perennis'), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  annotate('text', x=1997.75, y=18, label='(c) Lu. perennis', size=7, hjust=0)
peviplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Petalostemum.villosum' & year<2005), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  annotate('text', x=1997.75, y=18, label='(d) P. villosum', size=7, hjust=0)

#put hedges d panels together
pushViewport(viewport(layout=grid.layout(2,2))) 
print(amcaplotSpp, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lecaplotSpp, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(lupeplotSpp, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(peviplotSpp, vp=viewport(layout.pos.row=2, layout.pos.col=2))



#only look at 2000-2005 data for ESA talk, averaged
anppHedgesDSpp2 <- anppHedgesDSpp%>%
  mutate(spp2=ifelse(species=='Amorpha.canescens', 'AMCA', ifelse(species=='Lespedeza.capitata', 'LECA', ifelse(species=='Lupinus.perennis','LUPE', 'PEVI'))))

ggplot(data=barGraphStats(data=subset(anppHedgesDSpp2, year>1999 & trt!='hedgesd_CO2_N' & species!='Petalostemum.villosum'), variable="hedges_d", byFactorNames=c("spp2", "trt")), aes(x=spp2, y=mean, fill=interaction(trt,spp2))) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#00760a', '#00760a', '#e46c0a', '#e46c0a')) +
  ylab('Hedges D Aboveground Biomass') +
  xlab('') +
  theme(legend.position='none') +
  geom_hline(yintercept=0)







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
  filter(CO2_trt=='Cenrich'&N_trt=='Namb'&mass_total_ind!='is.na')%>%
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
  filter(CO2_trt=='Cenrich'&N_trt=='Namb'&mass_total_ind!='is.na')%>%
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
  filter(CO2_trt=='Cenrich'&N_trt=='Namb'&mass_total_ind!='is.na')%>%
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












