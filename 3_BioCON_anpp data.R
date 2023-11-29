################################################################################
##  3_BioCON_anpp data.R: Structuring ANPP data from BioCON experimental plots.
##
##  Author: Kimberly Komatsu
################################################################################

# library(boot)
library(grid)
library(nlme) 
library(lsmeans)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

#### Source home-built functions ####
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\general-R-code\\bar graph summary stats.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\general-R-code\\ggplot_theme set.r', chdir=T)

#### Source treatment data ####
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\1_BioCON_treatment data.r', chdir=T)

#### Read in ANPP data from 1998-2020
anpp <- read.table('e141_Plant aboveground biomass data_1998-2020.txt', sep="\t", header=TRUE) %>%
  rename(plot=Plot, ring=Ring) %>% 
  select(-Sampling.., -CO2.Treatment, -Nitrogen.Treatment, -CountOfSpecies, -CountOfGroup, -Experiment, -monospecies, -Monogroup) %>% 
  #convert date column to something useful
  left_join(trt) %>% 
  mutate(date=as.Date(as.character(Date), format='%m/%d/%Y'),
         year=as.numeric(format(date, '%Y')),
         month=as.factor(format(date, '%B'))) %>% 
  #remove trailing spaces, fixing capitalization, replacing . with space
  mutate(Species=str_squish(Species),
         Species=str_to_sentence(Species),
         monospecies=str_squish(monospecies),
         monospecies=str_to_sentence(monospecies),
         monospecies=chartr(".", " ", monospecies)) %>% 
  #remove non-vascular and dead biomass
  filter(!(Species %in% c('Moss', 'Oak leaves', 'Bare ground', 'Mosses & lichens', 'Miscellaneous litter'))) %>%   
  #rename species to monospecies if labeled as unsorted/green biomass in a monoculture
  mutate(Species=ifelse(spp_count==1 & Species %in% c('Unsorted biomass', 'Green biomass'),
                        monospecies, Species)) %>% 
  #rename species to other if not one of the focal species for the given plot
  rowwise() %>% 
  mutate(species2=ifelse(grepl(Species, spp_trt), Species, "other")) %>% 
  #drop 0 spp plots and terraCON plots
  filter(spp_count>0,
         !(Water.Treatment. %in% c('H2Oamb','H2Oneg')),
         !(Temp.Treatment. %in% c('HTamb','HTelev'))) %>% 
  #get max value across spring and fall clipping
  group_by(plot, ring, CO2_trt, N_trt, spp_count, group_count, experiment, monospecies, monogroup, year, species2) %>% 
  summarise(anpp=max(Aboveground.Biomass..g.m.2.)) %>% 
  ungroup()


#### Calculate Relative Yield and Relative Yield Total for each plot ####

#subset monocultures and drop non-focal species biomass
anppMono <- anpp %>% 
  filter(spp_count==1, species2!='other') %>% 
  #calculate avearge species monoculture biomass for each year and CO2/N trt
  group_by(year, CO2_trt, N_trt, species2) %>% 
  summarise(mono_anpp_mean=mean(anpp)) %>% 
  ungroup() %>% 
  #remove species with 0 biomass in monoculture
  filter(mono_anpp_mean>0)

#subset polycultures and calculate Relative Yield
anppPoly <- anpp %>% 
  filter(spp_count>1, species2!='other') %>% 
  #remove species without biomass
  filter(anpp>0) %>% 
  #join monoculture data
  left_join(anppMono) %>% 
  #calculate relative yield for each species
  mutate(RY=ifelse(spp_count==4, anpp/(mono_anpp_mean/4),
            ifelse(spp_count==9, anpp/(mono_anpp_mean/9),
                   anpp/(mono_anpp_mean/16)))) %>% 
  left_join(trt)

#calculate Relative Yield Total and total anpp per plot
RYT <- anppPoly %>% 
  group_by(plot, ring, CO2_trt, N_trt, spp_count, group_count, experiment, monospecies, monogroup, year) %>% 
  summarise(RYT=sum(RY), #calculate RYT
            total_anpp=sum(anpp)) %>% #calculate total anpp per plot 
  ungroup() %>% 
  left_join(trt)

#check RYT and total anpp for normality within plots
shapiro.test(anppPoly$RY)
qqnorm(anppPoly$RY)

shapiro.test(RYT$RYT)
qqnorm(RYT$RYT)

shapiro.test(RYT$total_anpp)
qqnorm(RYT$total_anpp)

#log10 transform for normality
anppPoly$logRY <- log10(anppPoly$RY)
shapiro.test(anppPoly$logRY)
qqnorm(anppPoly$logRY)

RYT$logRYT <- log10(RYT$RYT)
shapiro.test(RYT$logRYT)
qqnorm(RYT$logRYT)

RYT$log_total_anpp <- log10(RYT$total_anpp)
shapiro.test(RYT$log_total_anpp)
qqnorm(RYT$log_total_anpp)


#### Choose data subset ####
anppPolySubset <- subset(anppPoly, trt=='Camb_Namb' & year==2005)
RYTSubset <- subset(RYT, trt=='Camb_Namb' & year==2005)


#### Statistical Models - Relative Yield of legumes only ####
#note: subset only 2005 because that is the last year the plots were sorted to species in the dataset; drop Petalostemum because we don't have pot experiment data for that species (and it grows poorly in the plots as well)
summary(RY4spp <- lme(logRY ~ species2,
                  random=~1|plot,
                  data=anppPolySubset))
anova(RY4spp)
lsmeans(RY4spp, ~species2)

#plot rhizobial specialization
ggplot(data=barGraphStats(data=subset(anppPoly, monogroup=='Legume' & trt=='Camb_Namb' & year==2005 & species2!='Petalostemum villosum'), 
                          variable="RY", byFactorNames=c("species2")),
       aes(x=species2, y=mean, fill=species2)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Relative Yield in Field') +
  annotate('text', x=1, y=0.7, label='a', size=10) +
  annotate('text', x=2, y=0.6, label='a', size=10) +
  annotate('text', x=3, y=3, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('Amorpha canescens','Lespedeza capitata', 'Lupinus perennis'), labels=c('specialist', 'specialist', 'generalist'))
#export at 700x500



#with eCO2
#AMCA
amcaBio <- subset(anppRY05, subset=(spp=='AMCA_RY' & trt=='Cenrich_Namb'))
amcaRYBootModel <- boot(amcaBio, RYfunction, R=1000)
plot(amcaRYBootModel)
amcaRYBoot <- as.data.frame(amcaRYBootModel$t)
amcaRYBootMean <- as.data.frame(rowMeans(amcaRYBoot[1:1000,]))%>%
  mutate(spp='AMCA')
names(amcaRYBootMean)[names(amcaRYBootMean) == 'rowMeans(amcaRYBoot[1:1000, ])'] <- 'RY'

#LECA
lecaBio <- subset(anppRY05, subset=(spp=='LECA_RY' & trt=='Cenrich_Namb'))
lecaRYBootModel <- boot(lecaBio, RYfunction, R=1000)
plot(lecaRYBootModel)
lecaRYBoot <- as.data.frame(lecaRYBootModel$t)
lecaRYBootMean <- as.data.frame(rowMeans(lecaRYBoot[1:1000,]))%>%
  mutate(spp='LECA')
names(lecaRYBootMean)[names(lecaRYBootMean) == 'rowMeans(lecaRYBoot[1:1000, ])'] <- 'RY'

#LUPE
lupeBio <- subset(anppRY05, subset=(spp=='LUPE_RY' & trt=='Cenrich_Namb'))
lupeRYBootModel <- boot(lupeBio, RYfunction, R=1000)
plot(lupeRYBootModel)
lupeRYBoot <- as.data.frame(lupeRYBootModel$t)
lupeRYBootMean <- as.data.frame(rowMeans(lupeRYBoot[1:1000,]))%>%
  mutate(spp='LUPE')
names(lupeRYBootMean)[names(lupeRYBootMean) == 'rowMeans(lupeRYBoot[1:1000, ])'] <- 'RY'

#combine spp specializations
allSpecializationBoot <- rbind(amcaRYBootMean, lecaRYBootMean, lupeRYBootMean)%>%
  group_by(spp)%>%
  summarize(RY_mean=mean(RY), RY_sd=sd(RY))%>%
  ungroup()%>%
  mutate(RY_CI=1.96*RY_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=RY_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RY_mean-RY_CI, ymax=RY_mean+RY_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Relative Yield in Field') +
  annotate('text', x=1, y=0.7, label='a', size=10) +
  annotate('text', x=2, y=3.5, label='ab', size=10) +
  annotate('text', x=3, y=3.6, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500


##with eN
#AMCA
amcaBio <- subset(anppRY05, subset=(spp=='AMCA_RY' & trt=='Camb_Nenrich'))
amcaRYBootModel <- boot(amcaBio, RYfunction, R=1000)
plot(amcaRYBootModel)
amcaRYBoot <- as.data.frame(amcaRYBootModel$t)
amcaRYBootMean <- as.data.frame(rowMeans(amcaRYBoot[1:1000,]))%>%
  mutate(spp='AMCA')
names(amcaRYBootMean)[names(amcaRYBootMean) == 'rowMeans(amcaRYBoot[1:1000, ])'] <- 'RY'

#LECA
lecaBio <- subset(anppRY05, subset=(spp=='LECA_RY' & trt=='Camb_Nenrich'))
lecaRYBootModel <- boot(lecaBio, RYfunction, R=1000)
plot(lecaRYBootModel)
lecaRYBoot <- as.data.frame(lecaRYBootModel$t)
lecaRYBootMean <- as.data.frame(rowMeans(lecaRYBoot[1:1000,]))%>%
  mutate(spp='LECA')
names(lecaRYBootMean)[names(lecaRYBootMean) == 'rowMeans(lecaRYBoot[1:1000, ])'] <- 'RY'

#LUPE
lupeBio <- subset(anppRY05, subset=(spp=='LUPE_RY' & trt=='Camb_Nenrich'))
lupeRYBootModel <- boot(lupeBio, RYfunction, R=1000)
plot(lupeRYBootModel)
lupeRYBoot <- as.data.frame(lupeRYBootModel$t)
lupeRYBootMean <- as.data.frame(rowMeans(lupeRYBoot[1:1000,]))%>%
  mutate(spp='LUPE')
names(lupeRYBootMean)[names(lupeRYBootMean) == 'rowMeans(lupeRYBoot[1:1000, ])'] <- 'RY'

#combine spp specializations
allSpecializationBoot <- rbind(amcaRYBootMean, lecaRYBootMean, lupeRYBootMean)%>%
  group_by(spp)%>%
  summarize(RY_mean=mean(RY), RY_sd=sd(RY))%>%
  ungroup()%>%
  mutate(RY_CI=1.96*RY_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=RY_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RY_mean-RY_CI, ymax=RY_mean+RY_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Relative Yield in Field') +
  annotate('text', x=1, y=0.8, label='a', size=10) +
  annotate('text', x=2, y=1.3, label='a', size=10) +
  annotate('text', x=3, y=8.4, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500





##############################################
##############################################

#RY of legumes in monoculture vs polyculture, with and without other legumes
anppPolyRY <- merge(anppPoly, anppMonoAvg, all=T)
anppPolyRY$RY <- with(anppPolyRY, ifelse(spp_count==4, (anpp/(mono_anpp/4)),
                                         ifelse(spp_count==9, (anpp/(mono_anpp/9)), (anpp/(mono_anpp/16)))))

anppPolyRY <- anppPolyRY[complete.cases(anppPolyRY),]

#make column identifying if a species is a legume or non-legume
anppPolyRY$spp_type <- with(anppPolyRY, ifelse(species=='Amorpha.canescens', 'legume',
                                               ifelse(species=='Lespedeza.capitata', 'legume',
                                                      ifelse(species=='Lupinus.perennis', 'legume',
                                                             ifelse(species=='Petalostemum.villosum', 'legume', 'non-legume')))))

anppPolyRY$other_legs <- with(anppPolyRY, as.factor(ifelse(leg_num_spp=='0 legumes', '0 legumes',
                                                           ifelse(leg_num_spp=='2 legumes', '2 legumes',
                                                                  ifelse(leg_num_spp=='4 legumes', '4 legumes', 'single legume')))))

#normalize data

qqnorm(anppPolyRY$RY)
anppPolyRY$logRY <- log10(anppPolyRY$RY)
qqnorm(anppPolyRY$logRY)

##############################################
##############################################

#scale up RY to spp_num*RY
anppPolyRY$RYscaled <- with(anppPolyRY, spp_count*RY)
anppPolyRY$RYscaled_log <- log10(anppPolyRY$RYscaled)

#normality for anpp
qqnorm(anppPolyRY$anpp)
anppPolyRY$log_anpp <- log10(anppPolyRY$anpp)
qqnorm(anppPolyRY$log_anpp)

##############################################
##############################################

#difference between 1 and 4 legume plots
legume4 <- subset(anppPolyRY, spp_count==4 & other_legs=='4 legumes')
legume4mean <- with(legume4, aggregate(RY, list(species=species, year=year, CO2_trt=CO2_trt, N_trt=N_trt, trt=trt), mean))
names(legume4mean)[names(legume4mean)=='x'] <- 'RY4'
legume4mean$logRY4 <- log10(legume4mean$RY4)

legume1 <- subset(anppPolyRY, spp_count==4 & other_legs=='single legume')
legume1mean <- with(legume1, aggregate(RY, list(species=species, year=year, CO2_trt=CO2_trt, N_trt=N_trt, trt=trt), mean))
names(legume1mean)[names(legume1mean)=='x'] <- 'RY1'
legume1mean$logRY1 <- log10(legume1mean$RY1)

legume4to1 <- merge(legume4mean, legume1mean, all=T)
legume4to1complete <- legume4to1[complete.cases(legume4to1[,c(7,9)]),]
legume4to1complete$diff <- legume4to1complete$RY4-legume4to1complete$RY1
legume4to1complete$diff_log <- legume4to1complete$logRY4-legume4to1complete$logRY1

shapiro.test(legume4to1complete$diff)
qqnorm(legume4to1complete$diff)

shapiro.test(legume4to1complete$diff_log)
qqnorm(legume4to1complete$diff_log)

##############################################
##############################################

#get anpp hedges d
anppTrtMeans <- ddply(subset(anppPolyRY, spp_count==4), c('year', 'spp_type', 'trt'), summarise,
                      anpp_mean=mean(anpp),
                      anpp_sd=sd(anpp),
                      anpp_n=length(anpp))

anppRRmean <- dcast(anppTrtMeans, year + spp_type ~ trt, value.var='anpp_mean')
anppRRsd <- dcast(anppTrtMeans, year + spp_type ~ trt, value.var='anpp_sd')
names(anppRRsd)[names(anppRRsd)=='Camb_Namb'] <- 'ctl_sd'
names(anppRRsd)[names(anppRRsd)=='Camb_Nenrich'] <- 'N_sd'
names(anppRRsd)[names(anppRRsd)=='Cenrich_Namb'] <- 'CO2_sd'
names(anppRRsd)[names(anppRRsd)=='Cenrich_Nenrich'] <- 'CO2_N_sd'
anppRRn <- dcast(anppTrtMeans, year + spp_type ~ trt, value.var='anpp_n')
names(anppRRn)[names(anppRRn)=='Camb_Namb'] <- 'ctl_n'
names(anppRRn)[names(anppRRn)=='Camb_Nenrich'] <- 'N_n'
names(anppRRn)[names(anppRRn)=='Cenrich_Namb'] <- 'CO2_n'
names(anppRRn)[names(anppRRn)=='Cenrich_Nenrich'] <- 'CO2_N_n'
anppRRmeansd <- merge(anppRRmean, anppRRsd)
anppRR <- merge(anppRRmeansd, anppRRn)
anppRR$hedgesd_CO2 <- with(anppRR, ((Cenrich_Namb-Camb_Namb)/sqrt((CO2_sd*(CO2_n-1)+ctl_sd*(ctl_n-1))/(CO2_n+ctl_n-2)))*(1-(3/(4*(CO2_n+ctl_n-2)-1))))
anppRR$hedgesd_N <- with(anppRR, ((Camb_Nenrich-Camb_Namb)/sqrt((N_sd*(N_n-1)+ctl_sd*(ctl_n-1))/(N_n+ctl_n-2)))*(1-(3/(4*(N_n+ctl_n-2)-1))))
anppRR$hedgesd_CO2_N <- with(anppRR, ((Cenrich_Nenrich-Camb_Namb)/sqrt((CO2_N_sd*(CO2_N_n-1)+ctl_sd*(ctl_n-1))/(CO2_N_n+ctl_n-2)))*(1-(3/(4*(CO2_N_n+ctl_n-2)-1))))
anppRR$hedgesd_var_CO2 <- with(anppRR, ((CO2_n+ctl_n)/(CO2_n*ctl_n))+((hedgesd_CO2^2)/(2*(CO2_n+ctl_n))))
anppRR$hedgesd_var_N <- with(anppRR, ((N_n+ctl_n)/(N_n*ctl_n))+((hedgesd_N^2)/(2*(N_n+ctl_n))))
anppRR$hedgesd_var_CO2_N <- with(anppRR, ((CO2_N_n+ctl_n)/(CO2_N_n*ctl_n))+((hedgesd_CO2_N^2)/(2*(CO2_N_n+ctl_n))))
anppRR$hedgesd_ci_CO2 <- 1.96*anppRR$hedgesd_var_CO2
anppRR$hedgesd_ci_N <- 1.96*anppRR$hedgesd_var_N
anppRR$hedgesd_ci_CO2_N <- 1.96*anppRR$hedgesd_var_CO2_N

##############################################
##############################################

#get anpp hedges d by legume spp
anppTrtMeansSpp <- ddply(subset(anppPolyRY, spp_count==4 & spp_type=='legume'), c('year', 'species', 'trt'), summarise,
                         anpp_mean=mean(anpp),
                         anpp_sd=sd(anpp),
                         anpp_n=length(anpp))

anppRRmeanSpp <- dcast(anppTrtMeansSpp, year + species ~ trt, value.var='anpp_mean')
anppRRsdSpp <- dcast(anppTrtMeansSpp, year + species ~ trt, value.var='anpp_sd')
names(anppRRsdSpp)[names(anppRRsdSpp)=='Camb_Namb'] <- 'ctl_sd'
names(anppRRsdSpp)[names(anppRRsdSpp)=='Camb_Nenrich'] <- 'N_sd'
names(anppRRsdSpp)[names(anppRRsdSpp)=='Cenrich_Namb'] <- 'CO2_sd'
names(anppRRsdSpp)[names(anppRRsdSpp)=='Cenrich_Nenrich'] <- 'CO2_N_sd'
anppRRnSpp <- dcast(anppTrtMeansSpp, year + species ~ trt, value.var='anpp_n')
names(anppRRnSpp)[names(anppRRnSpp)=='Camb_Namb'] <- 'ctl_n'
names(anppRRnSpp)[names(anppRRnSpp)=='Camb_Nenrich'] <- 'N_n'
names(anppRRnSpp)[names(anppRRnSpp)=='Cenrich_Namb'] <- 'CO2_n'
names(anppRRnSpp)[names(anppRRnSpp)=='Cenrich_Nenrich'] <- 'CO2_N_n'
anppRRmeansdSpp <- merge(anppRRmeanSpp, anppRRsdSpp)
anppRRSpp <- merge(anppRRmeansdSpp, anppRRnSpp)
anppRRSpp$hedgesd_CO2 <- with(anppRRSpp, ((Cenrich_Namb-Camb_Namb)/sqrt((CO2_sd*(CO2_n-1)+ctl_sd*(ctl_n-1))/(CO2_n+ctl_n-2)))*(1-(3/(4*(CO2_n+ctl_n-2)-1))))
anppRRSpp$hedgesd_N <- with(anppRRSpp, ((Camb_Nenrich-Camb_Namb)/sqrt((N_sd*(N_n-1)+ctl_sd*(ctl_n-1))/(N_n+ctl_n-2)))*(1-(3/(4*(N_n+ctl_n-2)-1))))
anppRRSpp$hedgesd_CO2_N <- with(anppRRSpp, ((Cenrich_Nenrich-Camb_Namb)/sqrt((CO2_N_sd*(CO2_N_n-1)+ctl_sd*(ctl_n-1))/(CO2_N_n+ctl_n-2)))*(1-(3/(4*(CO2_N_n+ctl_n-2)-1))))
anppRRSpp$hedgesd_var_CO2 <- with(anppRRSpp, ((CO2_n+ctl_n)/(CO2_n*ctl_n))+((hedgesd_CO2^2)/(2*(CO2_n+ctl_n))))
anppRRSpp$hedgesd_var_N <- with(anppRRSpp, ((N_n+ctl_n)/(N_n*ctl_n))+((hedgesd_N^2)/(2*(N_n+ctl_n))))
anppRRSpp$hedgesd_var_CO2_N <- with(anppRRSpp, ((CO2_N_n+ctl_n)/(CO2_N_n*ctl_n))+((hedgesd_CO2_N^2)/(2*(CO2_N_n+ctl_n))))
anppRRSpp$hedgesd_ci_CO2 <- 1.96*anppRRSpp$hedgesd_var_CO2
anppRRSpp$hedgesd_ci_N <- 1.96*anppRRSpp$hedgesd_var_N
anppRRSpp$hedgesd_ci_CO2_N <- 1.96*anppRRSpp$hedgesd_var_CO2_N

anppHedgesDaSpp <- anppRRSpp[,-c(3:14,18:23)]
anppHedgesDSppd <- melt(anppHedgesDaSpp, id.vars=c('year', 'species'), variable.name='trt', value.name='hedges_d')
anppHedgesDbSpp <- anppRRSpp[,-c(3:20)]
anppHedgesDSppci <- melt(anppHedgesDbSpp, id.vars=c('year', 'species'), variable.name='trt', value.name='hedges_d_ci')
anppHedgesDSppci$trt2 <- with(anppHedgesDSppci, ifelse(trt=='hedgesd_ci_CO2', 'hedgesd_CO2',
                                                       ifelse(trt=='hedgesd_ci_N', 'hedgesd_N', 'hedgesd_CO2_N')))
anppHedgesDSppci$trt <- NULL
names(anppHedgesDSppci)[names(anppHedgesDSppci)=='trt2'] <- 'trt'
anppHedgesDSpp <- merge(anppHedgesDSppd, anppHedgesDSppci)







#boostrap RY values
RYfunction <- function(d, i){
  d2 <- d[i,]
  return(d2$RY)
}

#AMCA
amcaBio <- subset(anppRY05, subset=(spp=='AMCA_RY'))
amcaRYBootModel <- boot(amcaBio, RYfunction, R=1000)
plot(amcaRYBootModel)
amcaRYBoot <- as.data.frame(amcaRYBootModel$t)
amcaRYBootMean <- as.data.frame(rowMeans(amcaRYBoot[1:1000,]))%>%
  mutate(spp='AMCA')
names(amcaRYBootMean)[names(amcaRYBootMean) == 'rowMeans(amcaRYBoot[1:1000, ])'] <- 'RY'

#LECA
lecaBio <- subset(anppRY05, subset=(spp=='LECA_RY'))
lecaRYBootModel <- boot(lecaBio, RYfunction, R=1000)
plot(lecaRYBootModel)
lecaRYBoot <- as.data.frame(lecaRYBootModel$t)
lecaRYBootMean <- as.data.frame(rowMeans(lecaRYBoot[1:1000,]))%>%
  mutate(spp='LECA')
names(lecaRYBootMean)[names(lecaRYBootMean) == 'rowMeans(lecaRYBoot[1:1000, ])'] <- 'RY'

#LUPE
lupeBio <- subset(anppRY05, subset=(spp=='LUPE_RY'))
lupeRYBootModel <- boot(lupeBio, RYfunction, R=1000)
plot(lupeRYBootModel)
lupeRYBoot <- as.data.frame(lupeRYBootModel$t)
lupeRYBootMean <- as.data.frame(rowMeans(lupeRYBoot[1:1000,]))%>%
  mutate(spp='LUPE')
names(lupeRYBootMean)[names(lupeRYBootMean) == 'rowMeans(lupeRYBoot[1:1000, ])'] <- 'RY'

#combine spp specializations
allSpecializationBoot <- rbind(amcaRYBootMean, lecaRYBootMean, lupeRYBootMean)%>%
  group_by(spp)%>%
  summarize(RY_mean=mean(RY), RY_sd=sd(RY))%>%
  ungroup()%>%
  mutate(RY_CI=1.96*RY_sd)

#plot rhizobial specialization
ggplot(data=allSpecializationBoot, aes(x=spp, y=RY_mean, fill=spp)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RY_mean-RY_CI, ymax=RY_mean+RY_CI), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Relative Yield in Field') +
  annotate('text', x=1, y=0.6, label='a', size=10) +
  annotate('text', x=2, y=1.2, label='a', size=10) +
  annotate('text', x=3, y=4.7, label='b', size=10) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('AMCA', 'LUPE'), labels=c('specialist', 'generalist'))
#export at 700x500







##############################################
##############################################

#clean up workspace
rm(list=c('anpp', 'anpp16', 'anppAll', 'anppInitial', 'anppInitialYear', 'anppLast', 'anppLastLong', 'anppLastLongComplete', 'anppMono', 'anppMonoAvgWide', 'anppMonoAvg', 'anppPoly',  'anppNoTrt', 'anppSum', 'anppTrt', 'anppTrue16', 'anppTrue16Trt', 'anppTrue16TrtNonzero', 'trt', 'anppSub', 'legume1', 'legume1mean', 'legume4', 'legume4mean', 'legume4to1', 'anppTrtMeans', 'anppRRmean', 'anppRRmeansd', 'anppRRn', 'anppRRsd', 'anppHedgesDaSpp', 'anppHedgesDbSpp', 'anppHedgesDSppci', 'anppHedgesDSppd', 'anppRRmeansdSpp', 'anppRRmeanSpp', 'anppRRnSpp', 'anppRRSpp', 'anppTrtMeansSpp'))
