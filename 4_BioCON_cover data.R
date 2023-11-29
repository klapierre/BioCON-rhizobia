################################################################################
##  4_BioCON_cover data.R: Structuring cover data from BioCON experimental plots.
##
##  Author: Kimberly Komatsu
################################################################################

library(car)
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