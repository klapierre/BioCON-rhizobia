library(plyr)
library(reshape2)

setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_treatment data.R', chdir=T)

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

