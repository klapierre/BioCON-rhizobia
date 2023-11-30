################################################################################
##  4_BioCON_cover data.R: Structuring cover data from BioCON experimental plots.
##
##  Author: Kimberly Komatsu
################################################################################

library(car)
library(grid)
library(nlme) 
library(lsmeans)
library(tidyverse)


setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

#### Source home-built functions ####
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\general-R-code\\bar graph summary stats.r', chdir=T)
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\general-R-code\\ggplot_theme set.r', chdir=T)

#### Source anpp data ####
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\3_BioCON_anpp data.r', chdir=T)

#cleanup
rm(list=setdiff(ls(), "anpp"))

#### Source treatment data ####
source('C:\\Users\\kjkomatsu\\Desktop\\R files\\BioCON-rhizobia\\1_BioCON_treatment data.r', chdir=T)

#cover data
cover <- read.table('e141_Plant species percent cover data_1999-2019.txt', sep="\t", header=TRUE) %>%
  rename(plot=Plot, ring=Ring, year=Year) %>% 
  select(-Sampling.number, -CO2.Treatment, -Nitrogen.Treatment, -CountOfSpecies, -CountOfGroup, -Experiment, -Monospecies, -Monogroup) %>% 
  #convert date column to something useful
  left_join(trt) %>% 
  #remove trailing spaces, fixing capitalization, replacing . with space
  mutate(Species=str_squish(Species),
         Species=str_to_sentence(Species)) %>% 
  #remove non-vascular and dead biomass
  filter(!(Species %in% c('Moss', 'Oak leaves', 'Mushrooms', 'Mushroom', 'Bareground', 'Bare ground', 'Mosses & lichens', 'Miscellaneous litter', 'Micellaneous litter', 'Miscelaneous liter', 'Miscellaneous llitter', 'Fungi', ''))) %>% 
  filter(Percent.cover>0) %>% 
  #rename species to other if not one of the focal species for the given plot
  rowwise() %>% 
  mutate(species2=ifelse(grepl(Species, spp_trt), Species, "other")) %>% 
  #drop 0 spp plots and terraCON plots
  filter(spp_count>0,
         !(Water.Treatment %in% c('H2Oamb','H2Oneg')),
         !(Temp.Treatment %in% c('HTamb','HTelev'))) %>% 
  #get max value across spring and fall clipping
  group_by(plot, ring, CO2_trt, N_trt, spp_count, group_count, experiment, monospecies, monogroup, year, species2) %>% 
  summarise(cover=max(Percent.cover)) %>% 
  ungroup()


#### Merge cover and anpp to test correlations ####
coverANPP <- cover %>% 
  left_join(anpp)

#linear model comparing cover and ANPP across species
summary(coverANPPmodel <- lm(anpp~cover, data=coverANPP))
outlierTest(coverANPPmodel)
plot(coverANPP$anpp~coverANPP$cover)

#subset each species
coverANPPacmi <- subset(coverANPP, species2=='Achillea millefolium')
coverANPPacmiModel <- lm(anpp~cover, data=coverANPPacmi)
outlierTest(coverANPPacmiModel)
summary(coverANPPacmiModel)
plot(coverANPPacmi$anpp~coverANPPacmi$cover)

coverANPPastu <- subset(coverANPP, species2=='Asclepias tuberosa')
coverANPPastuModel <- lm(anpp~cover, data=coverANPPastu)
outlierTest(coverANPPastuModel)
summary(coverANPPastuModel)
plot(coverANPPastu$anpp~coverANPPastu$cover)

coverANPPkocr <- subset(coverANPP, species2=='Koeleria cristata')
coverANPPkocrModel <- lm(anpp~cover, data=coverANPPkocr)
outlierTest(coverANPPkocrModel)
summary(coverANPPkocrModel)
plot(coverANPPkocr$anpp~coverANPPkocr$cover)

coverANPPlupe <- subset(coverANPP, species2=='Lupinus perennis')
coverANPPlupeModel <- lm(anpp~cover, data=coverANPPlupe)
outlierTest(coverANPPlupeModel)
summary(coverANPPlupeModel)
plot(coverANPPlupe$anpp~coverANPPlupe$cover)

coverANPPpopr <- subset(coverANPP, species2=='Poa pratensis')
coverANPPpoprModel <- lm(anpp~cover, data=coverANPPpopr)
outlierTest(coverANPPpoprModel)
summary(coverANPPpoprModel)
plot(coverANPPpopr$anpp~coverANPPpopr$cover)

coverANPPsonu <- subset(coverANPP, species2=='Sorghastrum nutans') #NOTE: very poor fit
coverANPPsonuModel <- lm(anpp~cover, data=coverANPPsonu)
outlierTest(coverANPPsonuModel)
summary(coverANPPsonuModel)
plot(coverANPPsonu$anpp~coverANPPsonu$cover)

coverANPPagre <- subset(coverANPP, species2=='Agropyron repens')
coverANPPagreModel <- lm(anpp~cover, data=coverANPPagre)
outlierTest(coverANPPagreModel)
summary(coverANPPagreModel)
plot(coverANPPagre$anpp~coverANPPagre$cover)

coverANPPange <- subset(coverANPP, species2=='Andropogon gerardi')
coverANPPangeModel <- lm(anpp~cover, data=coverANPPange)
outlierTest(coverANPPangeModel)
summary(coverANPPangeModel)
plot(coverANPPange$anpp~coverANPPange$cover)

coverANPPbrin <- subset(coverANPP, species2=='Bromus inermis')
coverANPPbrinModel <- lm(anpp~cover, data=coverANPPbrin)
outlierTest(coverANPPbrinModel)
summary(coverANPPbrinModel)
plot(coverANPPbrin$anpp~coverANPPbrin$cover)

coverANPPpevi <- subset(coverANPP, species2=='Petalostemum villosum')
coverANPPpeviModel <- lm(anpp~cover, data=coverANPPpevi)
outlierTest(coverANPPpeviModel)
summary(coverANPPpeviModel)
plot(coverANPPpevi$anpp~coverANPPpevi$cover)

coverANPPamca <- subset(coverANPP, species2=='Amorpha canescens')
coverANPPamcaModel <- lm(anpp~cover, data=coverANPPamca)
outlierTest(coverANPPamcaModel)
summary(coverANPPamcaModel)
plot(coverANPPamca$anpp~coverANPPamca$cover)

coverANPPbogr <- subset(coverANPP, species2=='Bouteloua gracilis')
coverANPPbogrModel <- lm(anpp~cover, data=coverANPPbogr)
outlierTest(coverANPPbogrModel)
summary(coverANPPbogrModel)
plot(coverANPPbogr$anpp~coverANPPbogr$cover)

coverANPPscsc <- subset(coverANPP, species2=='Schizachyrium scoparium')
coverANPPscscModel <- lm(anpp~cover, data=coverANPPscsc)
outlierTest(coverANPPscscModel)
summary(coverANPPscscModel)
plot(coverANPPscsc$anpp~coverANPPscsc$cover)

coverANPPsori <- subset(coverANPP, species2=='Solidago rigida')
coverANPPsoriModel <- lm(anpp~cover, data=coverANPPsori)
outlierTest(coverANPPsoriModel)
summary(coverANPPsoriModel)
plot(coverANPPsori$anpp~coverANPPsori$cover)

coverANPPancy <- subset(coverANPP, species2=='Anemone cylindrica')
coverANPPancyModel <- lm(anpp~cover, data=coverANPPancy)
outlierTest(coverANPPancyModel)
summary(coverANPPancyModel)
plot(coverANPPancy$anpp~coverANPPancy$cover)

coverANPPleca <- subset(coverANPP, species2=='Lespedeza capitata')
coverANPPlecaModel <- lm(anpp~cover, data=coverANPPleca)
outlierTest(coverANPPlecaModel)
summary(coverANPPlecaModel)
plot(coverANPPleca$anpp~coverANPPleca$cover)

#total cover and anpp
coverANPPtotal <- coverANPP %>% 
  group_by(year, plot, ring, CO2_trt, N_trt, spp_count, group_count, monogroup, monospecies) %>% 
  summarise(anpp_total=sum(anpp),
            cover_total=sum(cover)) %>% 
  ungroup()

coverANPPtotalModel <- lm(anpp_total~cover_total, data=coverANPPtotal)
summary(coverANPPtotalModel)
plot(coverANPPtotal$anpp_total~coverANPPtotal$cover_total)

#####################################################################
#### Overall, cover is a poor proxy for biomass for most species ####
#####################################################################


# #### Calculate Relative Yield and Relative Yield Total for each plot ####
# 
# #subset monocultures and drop non-focal species biomass
# coverMono <- cover %>% 
#   filter(spp_count==1, species2!='other') %>% 
#   #calculate avearge species monoculture biomass for each year and CO2/N trt
#   group_by(year, CO2_trt, N_trt, species2) %>% 
#   summarise(mono_cover_mean=mean(cover)) %>% 
#   ungroup() %>% 
#   #remove species with 0 biomass in monoculture
#   filter(mono_cover_mean>0)
# 
# #get list of species that should be present in each plot
# sppList <- trt %>% 
#   pivot_longer(Achillea.millefolium:Sorghastrum.nutans, names_to='species', values_to='species2') %>% 
#   filter(species2!=0) %>% 
#   select(plot, ring, species2, spp_count) %>% 
#   unique()
# 
# #subset polycultures and calculate Relative Yield
# coverExpected <- sppList %>% 
#   left_join(trt) %>% 
#   select(plot, ring, CO2_trt, N_trt, species2, spp_count) %>% 
#   filter(spp_count>1) %>% 
#   right_join(coverMono, relationship='many-to-many') %>% 
#   #calculate relative yield for each species
#   mutate(expected_spp_cover=ifelse(spp_count==4, mono_cover_mean/4,
#                             ifelse(spp_count==9, mono_cover_mean/9,
#                                    mono_cover_mean/16))) %>% 
#   group_by(plot, ring, year) %>% 
#   summarize(expected_cover=sum(expected_spp_cover)) %>% 
#   ungroup()
# 
# coverPoly <- cover %>% 
#   select(plot, ring, year, CO2_trt, N_trt, species2, spp_count, cover) %>% 
#   filter(spp_count>1) %>% 
#   right_join(coverMono) %>% 
#   #calculate relative yield for each species
#   mutate(expected_spp_cover=ifelse(spp_count==4, mono_cover_mean/4,
#                             ifelse(spp_count==9, mono_cover_mean/9,
#                                    mono_cover_mean/16))) %>% 
#   mutate(RY=(cover/expected_spp_cover)) %>% 
#   left_join(trt)
# 
# #calculate Relative Yield Total and total anpp per plot
# RYT <- cover %>% 
#   filter(spp_count>1) %>% 
#   group_by(plot, ring, year) %>% 
#   summarize(total_cover=sum(cover)) %>% 
#   ungroup() %>% 
#   left_join(coverExpected) %>% 
#   mutate(RYT=total_cover/expected_cover) %>% 
#   left_join(trt)
# 
# #check RYT and total anpp for normality within plots
# # shapiro.test(anppPoly$RY)
# qqnorm(coverPoly$RY)
# 
# shapiro.test(RYT$RYT)
# qqnorm(RYT$RYT)
# 
# shapiro.test(RYT$total_cover)
# qqnorm(RYT$total_cover)
# 
# #log10 transform for normality
# coverPoly$logRY <- log10(coverPoly$RY)
# # shapiro.test(anppPoly$logRY)
# # qqnorm(anppPoly$logRY)
# 
# RYT$logRYT <- log10(RYT$RYT)
# shapiro.test(RYT$logRYT)
# qqnorm(RYT$logRYT)
# 
# RYT$log_total_cover <- log10(RYT$total_cover)
# shapiro.test(RYT$log_total_cover)
# qqnorm(RYT$log_total_cover)
# 
# 
# #### Choose data subset ####
# coverPolySubset <- subset(coverPoly, trt=='Camb_Namb')
# RYTSubset <- subset(RYT, trt=='Camb_Namb')
# 
# 
# #### Statistical Models - Relative Yield of legumes only ####
# #note: subset only 2005 because that is the last year the plots were sorted to species in the dataset; drop Petalostemum because we don't have pot experiment data for that species (and it grows poorly in the plots as well)
# summary(RYmodel <- lme(log(RY) ~ species2,
#                        random=~1|plot,
#                        data=subset(coverPolySubset, year==2005 &
#                                    species2 %in% c('Amorpha canescens', 
#                                                    'Lespedeza capitata', 
#                                                    'Lupinus perennis'))))
# anova(RYmodel)
# lsmeans(RYmodel, ~species2)
# 
# #through time
# summary(RYmodel <- lme(log(RY) ~ species2,
#                        random=~1|year/plot,
#                        data=subset(coverPolySubset, species2 %in% c('Amorpha canescens',
#                                                                     'Lespedeza capitata',
#                                                                     'Lupinus perennis'))))
# anova(RYmodel)
# lsmeans(RYmodel, ~species2)
# 
# #plot rhizobial specialization by year
# ggplot(data=barGraphStats(data=subset(coverPolySubset, monogroup=='Legume' &
#                                         species2 %in% c('Amorpha canescens','Lespedeza capitata','Lupinus perennis')),
#                           variable="RY", byFactorNames=c("year", "species2")),
#        aes(x=year, y=mean, color=species2)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#   geom_hline(aes(yintercept=1), linetype="dashed") +
#   xlab('Year') +
#   ylab('Relative Yield in Field')
#  
# ggplot(data=barGraphStats(data=subset(coverPolySubset, monogroup=='Legume' & year==2005 &
#                                         species2 %in% c('Amorpha canescens','Lespedeza capitata','Lupinus perennis')),
#                           variable="RY", byFactorNames=c("species2")),
#        aes(x=species2, y=mean, fill=species2)) +
#   geom_bar(stat='identity') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2) +
#   geom_hline(aes(yintercept=1), linetype="dashed") +
#   xlab('Legume Species') +
#   ylab('Relative Yield in Field') +
#   annotate('text', x=1, y=0.8, label='a', size=10) +
#   annotate('text', x=2, y=0.78, label='a', size=10) +
#   annotate('text', x=3, y=3.6, label='b', size=10) +
#   scale_x_discrete(breaks=c('Amorpha canescens','Lespedeza capitata', 'Lupinus perennis'), labels=c('AMCA', 'LECA', 'LUPE')) +
#   scale_fill_manual(values=c('#00760a', '#00760a', '#e46c0a'), breaks=c('Amorpha canescens','Lespedeza capitata', 'Lupinus perennis'), labels=c('specialist', 'specialist', 'generalist'))
# #export at 700x500
# 
# 
# #### Statistical Models - Relative Yield Total ####
# #note: subset only 2005 because that is the last year the plots were sorted to species in the dataset
# 
# #4 spp polycultures only
# summary(RYTmodel <- lme(log10(RYT) ~ as.factor(legume_num),
#                         random=~1|year/plot,
#                         data=subset(RYTSubset, legume_num %in% c(0,1,4) & spp_count==4)))
# anova(RYTmodel)
# lsmeans(RYTmodel, ~legume_num)
# 
# #plot rhizobial specialization by year
# ggplot(data=barGraphStats(data=subset(RYTSubset, legume_num %in% c(0,1,4) & spp_count==4),
#                           variable="RYT", byFactorNames=c("year", "legume_num")),
#        aes(x=year, y=mean, color=as.factor(legume_num))) +
#   geom_point(position=position_dodge(0.9)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
#   geom_hline(aes(yintercept=1), linetype="dashed") +
#   xlab('Year') +
#   ylab('Relative Yield in Field')
# 
# ggplot(data=barGraphStats(data=subset(RYTSubset, legume_num %in% c(0,1,4) & spp_count==4), variable="RYT", byFactorNames=c("legume_num")),
#        aes(x=as.factor(legume_num), y=mean)) +
#   geom_bar(stat='identity') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=0.2) +
#   geom_hline(aes(yintercept=1), linetype="dashed") +
#   xlab('Legume Number') +
#   ylab('Relative Yield Total') +
#   annotate('text', x=1, y=2.0, label='a', size=10) +
#   annotate('text', x=2, y=2.0, label='a', size=10) +
#   annotate('text', x=3, y=1.75, label='b', size=10)
# #export at 700x500