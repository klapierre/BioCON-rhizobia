library(plyr)
library(reshape2)
library(car)

setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_treatment data.R', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_anpp data.R', chdir=T)

#cover data
coverAll <- read.csv('e141_plant species percent cover.csv')

#remove the water and temperature treatment plots
coverWater <- subset(coverAll, subset=(water_trt!='H2Oamb' & water_trt!='H2Oneg'))

#remove all but the 16 species cover
cover16 <- coverWater[(coverWater$species %in% c("Achillea millefolium","Asclepias tuberosa",
                                             "Koeleria cristata","Lupinus perennis","Poa pratensis",
                                             "Sorghastrum nutans","Agropyron repens","Andropogon gerardi",
                                             "Bromus inermis","Petalostemum villosum","Amorpha canescens",
                                             "Bouteloua gracilis","Schizachyrium scoparium", "Solidago rigida",
                                             "Anemone cylindrica", "Lespedeza capitata")),]

#get max per year (across spring and fall)
coverMax <- ddply(cover16, c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup', 'species'), summarise, cover=max(cover))

coverMax$species <- as.factor(coverMax$species)

#make anpp species names match cover species names
anppMax$species <- as.factor(with(anppMax, ifelse(species=='Achillea.millefolium', 'Achillea millefolium',
                                 ifelse(species=='Agropyron.repens', 'Agropyron repens',
                                 ifelse(species=='Amorpha.canescens', 'Amorpha canescens',
                                 ifelse(species=='Andropogon.gerardi', 'Andropogon gerardi',
                                 ifelse(species=='Anemone.cylindrica', 'Anemone cylindrica',
                                 ifelse(species=='Asclepias.tuberosa', 'Asclepias tuberosa',
                                 ifelse(species=='Bouteloua.gracilis', 'Bouteloua gracilis',
                                 ifelse(species=='Bromus.inermis', 'Bromus inermis',
                                 ifelse(species=='Koeleria.cristata', 'Koeleria cristata',
                                 ifelse(species=='Lespedeza.capitata', 'Lespedeza capitata',
                                 ifelse(species=='Lupinus.perennis', 'Lupinus perennis',
                                 ifelse(species=='Petalostemum.villosum', 'Petalostemum villosum',
                                 ifelse(species=='Poa.pratensis', 'Poa pratensis',
                                 ifelse(species=='Schizachyrium.scoparium', 'Schizachyrium scoparium',
                                 ifelse(species=='Solidago.rigida', 'Solidago rigida', 'Sorghastrum nutans')))))))))))))))))

#merge cover and anpp
coverANPP <- merge(coverMax, anppMax, by=c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup', 'species'))

#linear model comparing cover and ANPP across species
coverANPPmodel <- lm(anpp~cover, data=coverANPP)
outlierTest(coverANPPmodel)
summary(coverANPPmodel)
plot(coverANPP$anpp~coverANPP$cover)

#subset each species
coverANPPacmi <- subset(coverANPP, species=='Achillea millefolium')
coverANPPacmiModel <- lm(anpp~cover, data=coverANPPacmi)
outlierTest(coverANPPacmiModel)
summary(coverANPPacmiModel)
plot(coverANPPacmi$anpp~coverANPPacmi$cover)

coverANPPastu <- subset(coverANPP, species=='Asclepias tuberosa')
coverANPPastuModel <- lm(anpp~cover, data=coverANPPastu)
outlierTest(coverANPPastuModel)
summary(coverANPPastuModel)
plot(coverANPPastu$anpp~coverANPPastu$cover)

coverANPPkocr <- subset(coverANPP, species=='Koeleria cristata')
coverANPPkocrModel <- lm(anpp~cover, data=coverANPPkocr)
outlierTest(coverANPPkocrModel)
summary(coverANPPkocrModel)
plot(coverANPPkocr$anpp~coverANPPkocr$cover)

coverANPPlupe <- subset(coverANPP, species=='Lupinus perennis')
coverANPPlupeModel <- lm(anpp~cover, data=coverANPPlupe)
outlierTest(coverANPPlupeModel)
summary(coverANPPlupeModel)
plot(coverANPPlupe$anpp~coverANPPlupe$cover)

coverANPPpopr <- subset(coverANPP, species=='Poa pratensis')
coverANPPpoprModel <- lm(anpp~cover, data=coverANPPpopr)
outlierTest(coverANPPpoprModel)
summary(coverANPPpoprModel)
plot(coverANPPpopr$anpp~coverANPPpopr$cover)

coverANPPsonu <- subset(coverANPP, species=='Sorghastrum nutans') #NOTE: very poor fit
coverANPPsonuModel <- lm(anpp~cover, data=coverANPPsonu)
outlierTest(coverANPPsonuModel)
summary(coverANPPsonuModel)
plot(coverANPPsonu$anpp~coverANPPsonu$cover)

coverANPPagre <- subset(coverANPP, species=='Agropyron repens')
coverANPPagreModel <- lm(anpp~cover, data=coverANPPagre)
outlierTest(coverANPPagreModel)
summary(coverANPPagreModel)
plot(coverANPPagre$anpp~coverANPPagre$cover)

coverANPPange <- subset(coverANPP, species=='Andropogon gerardi')
coverANPPangeModel <- lm(anpp~cover, data=coverANPPange)
outlierTest(coverANPPangeModel)
summary(coverANPPangeModel)
plot(coverANPPange$anpp~coverANPPange$cover)

coverANPPbrin <- subset(coverANPP, species=='Bromus inermis')
coverANPPbrinModel <- lm(anpp~cover, data=coverANPPbrin)
outlierTest(coverANPPbrinModel)
summary(coverANPPbrinModel)
plot(coverANPPbrin$anpp~coverANPPbrin$cover)

coverANPPpevi <- subset(coverANPP, species=='Petalostemum villosum')
coverANPPpeviModel <- lm(anpp~cover, data=coverANPPpevi)
outlierTest(coverANPPpeviModel)
summary(coverANPPpeviModel)
plot(coverANPPpevi$anpp~coverANPPpevi$cover)

coverANPPamca <- subset(coverANPP, species=='Amorpha canescens')
coverANPPamcaModel <- lm(anpp~cover, data=coverANPPamca)
outlierTest(coverANPPamcaModel)
summary(coverANPPamcaModel)
plot(coverANPPamca$anpp~coverANPPamca$cover)

coverANPPbogr <- subset(coverANPP, species=='Bouteloua gracilis')
coverANPPbogrModel <- lm(anpp~cover, data=coverANPPbogr)
outlierTest(coverANPPbogrModel)
summary(coverANPPbogrModel)
plot(coverANPPbogr$anpp~coverANPPbogr$cover)

coverANPPscsc <- subset(coverANPP, species=='Schizachyrium scoparium')
coverANPPscscModel <- lm(anpp~cover, data=coverANPPscsc)
outlierTest(coverANPPscscModel)
summary(coverANPPscscModel)
plot(coverANPPscsc$anpp~coverANPPscsc$cover)

coverANPPsori <- subset(coverANPP, species=='Solidago rigida')
coverANPPsoriModel <- lm(anpp~cover, data=coverANPPsori)
outlierTest(coverANPPsoriModel)
summary(coverANPPsoriModel)
plot(coverANPPsori$anpp~coverANPPsori$cover)

coverANPPancy <- subset(coverANPP, species=='Anemone cylindrica')
coverANPPancyModel <- lm(anpp~cover, data=coverANPPancy)
outlierTest(coverANPPancyModel)
summary(coverANPPancyModel)
plot(coverANPPancy$anpp~coverANPPancy$cover)

coverANPPleca <- subset(coverANPP, species=='Lespedeza capitata')
coverANPPlecaModel <- lm(anpp~cover, data=coverANPPleca)
outlierTest(coverANPPlecaModel)
summary(coverANPPlecaModel)
plot(coverANPPleca$anpp~coverANPPleca$cover)

#total cover and anpp
coverANPPtotal <- ddply(coverANPP, c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'monospecies', 'monospecies'), summarise,
                        anpp_total=mean(anpp),
                        cover_total=mean(cover))

coverANPPtotalModel <- lm(anpp_total~cover_total, data=coverANPPtotal)
summary(coverANPPtotalModel)
plot(coverANPPtotal$anpp_total~coverANPPtotal$cover_total)



##############################################
##############################################

#clean up workspace
rm(list=c('cover16','coverAll','coverANPPacmi','coverANPPacmiModel','coverANPPagre','coverANPPagreModel','coverANPPamca','coverANPPamcaModel','coverANPPancy','coverANPPancyModel','coverANPPange','coverANPPangeModel','coverANPPastu','coverANPPastuModel','coverANPPbogr','coverANPPbogrModel','coverANPPbrin','coverANPPbrinModel','coverANPPkocr','coverANPPkocrModel','coverANPPleca','coverANPPlecaModel','coverANPPlupe','coverANPPlupeModel','coverANPPmodel','coverANPPpevi','coverANPPpeviModel','coverANPPpopr','coverANPPpoprModel','coverANPPscsc','coverANPPscscModel','coverANPPsonu','coverANPPsonuModel','coverANPPsori','coverANPPsoriModel','coverANPPtotal','coverANPPtotalModel','coverMax','coverWater'))





