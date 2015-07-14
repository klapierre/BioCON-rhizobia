setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

#cover data
coverAll <- read.csv('e141_plant species percent cover.csv')

#get max per year (across spring and fall)
coverMax <- ddply(coverAll, c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup', 'species'), summarise, cover=max(cover))

