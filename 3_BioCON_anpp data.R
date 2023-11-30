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

#get list of species that should be present in each plot
sppList <- trt %>% 
  pivot_longer(Achillea.millefolium:Sorghastrum.nutans, names_to='species', values_to='species2') %>% 
  filter(species2!=0) %>% 
  select(plot, ring, species2, spp_count) %>% 
  unique()

#subset polycultures and calculate Relative Yield
anppExpected <- sppList %>% 
  left_join(trt) %>% 
  select(plot, ring, CO2_trt, N_trt, species2, spp_count) %>% 
  filter(spp_count>1) %>% 
  right_join(anppMono, relationship='many-to-many') %>% 
  #calculate relative yield for each species
  mutate(expected_spp_biomass=ifelse(spp_count==4, mono_anpp_mean/4,
                                     ifelse(spp_count==9, mono_anpp_mean/9,
                                            mono_anpp_mean/16))) %>% 
  group_by(plot, ring, year) %>% 
  summarize(expected_biomass=sum(expected_spp_biomass)) %>% 
  ungroup()

#calculate Relative Yield Total and total anpp per plot
RYT <- anpp %>% 
  filter(spp_count>1) %>% 
  group_by(plot, ring, year) %>% 
  summarize(total_biomass=sum(anpp)) %>% 
  ungroup() %>% 
  left_join(anppExpected) %>% 
  mutate(RYT=total_biomass/expected_biomass) %>% 
  left_join(trt)

#check RYT and total anpp for normality within plots
# shapiro.test(anppPoly$RY)
qqnorm(anppPoly$RY)

shapiro.test(RYT$RYT)
qqnorm(RYT$RYT)

shapiro.test(RYT$total_biomass)
qqnorm(RYT$total_biomass)

#log10 transform for normality
anppPoly$logRY <- log10(anppPoly$RY)
# shapiro.test(anppPoly$logRY)
# qqnorm(anppPoly$logRY)

RYT$logRYT <- log10(RYT$RYT)
shapiro.test(RYT$logRYT)
qqnorm(RYT$logRYT)

RYT$total_biomass <- log10(RYT$total_biomass)
shapiro.test(RYT$total_biomass)
qqnorm(RYT$total_biomass)


#### Choose data subset ####
anppPolySubset <- subset(anppPoly, trt=='Camb_Namb')
RYTSubset <- subset(RYT, trt=='Camb_Namb')


#### Statistical Models - Relative Yield of legumes only ####
#note: subset only 2005 because that is the last year the plots were sorted to species in the dataset; drop Petalostemum because we don't have pot experiment data for that species (and it grows poorly in the plots as well)
summary(RYmodel <- lme(log(RY) ~ species2,
                       random=~1|plot,
                       data=subset(anppPolySubset, year==2005 &
                                     species2 %in% c('Amorpha canescens', 
                                                     'Lespedeza capitata', 
                                                     'Lupinus perennis'))))
anova(RYmodel)
lsmeans(RYmodel, ~species2)

# #plot rhizobial specialization by year
# ggplot(data=barGraphStats(data=subset(anppPolySubset, monogroup=='Legume' &
#                                         species2 %in% c('Amorpha canescens','Lespedeza capitata','Lupinus perennis')), 
#                           variable="RY", byFactorNames=c("year", "species2")),
#        aes(x=year, y=mean, color=species2)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#   geom_hline(aes(yintercept=1), linetype="dashed") +
#   xlab('Year') +
#   ylab('Relative Yield in Field')
# 
ggplot(data=barGraphStats(data=subset(anppPolySubset, monogroup=='Legume' & year==2005 &
                                      species2 %in% c('Amorpha canescens','Lespedeza capitata','Lupinus perennis')),
                          variable="RY", byFactorNames=c("species2")),
       aes(x=species2, y=mean, fill=species2)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Species') +
  ylab('Relative Yield in Field') +
  annotate('text', x=1, y=0.8, label='a', size=10) +
  annotate('text', x=2, y=0.78, label='a', size=10) +
  annotate('text', x=3, y=3.6, label='b', size=10) +
  scale_x_discrete(breaks=c('Amorpha canescens','Lespedeza capitata', 'Lupinus perennis'), labels=c('AMCA', 'LECA', 'LUPE')) +
  scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('Amorpha canescens','Lespedeza capitata', 'Lupinus perennis'), labels=c('specialist', 'specialist', 'generalist'))
#export at 700x500


#### Statistical Models - Relative Yield Total ####
#note: subset only 2005 because that is the last year the plots were sorted to species in the dataset

#4 spp polycultures only
summary(RYTmodel <- lme(log10(RYT) ~ as.factor(legume_num),
                    random=~1|year/plot,
                    data=subset(RYTSubset, legume_num %in% c(0,1,4) & spp_count==4)))
anova(RYTmodel)
lsmeans(RYTmodel, ~legume_num)

# #plot rhizobial specialization by year
# ggplot(data=barGraphStats(data=subset(RYTSubset, legume_num %in% c(0,1,4) & spp_count==4),
#                           variable="RYT", byFactorNames=c("year", "legume_num")),
#        aes(x=year, y=mean, color=as.factor(legume_num))) +
#   geom_point(position=position_dodge(0.9)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
#   geom_hline(aes(yintercept=1), linetype="dashed") +
#   xlab('Year') +
#   ylab('Relative Yield in Field')

ggplot(data=barGraphStats(data=subset(RYTSubset, legume_num %in% c(0,1,4) & spp_count==4), variable="RYT", byFactorNames=c("legume_num")),
       aes(x=as.factor(legume_num), y=mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2) +
  geom_hline(aes(yintercept=1), linetype="dashed") +
  xlab('Legume Number') +
  ylab('Relative Yield Total') +
  annotate('text', x=1, y=2.0, label='a', size=10) +
  annotate('text', x=2, y=2.0, label='a', size=10) +
  annotate('text', x=3, y=1.75, label='b', size=10)
#export at 700x500