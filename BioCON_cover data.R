library(plyr)
library(reshape2)

setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_treatment data.R', chdir=T)

#cover data
coverAll <- read.csv('e141_plant species percent cover.csv')

#get max per year (across spring and fall)
coverMax <- ddply(coverAll, c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup', 'species'), summarise, cover=max(cover))

