setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

#################################################
#################################################

#import data 

#ANPP data (IMPORTANT: be careful that extra spaces are removed from after species names if data is newly downloaded from web)
anppInitial <- read.csv('e141_plant aboveground biomass_1998-2012.csv')

#ANPP data from 2013-2014
anppLast <- read.csv('e141_plant aboveground biomass_2013-2014.csv')

#cover data
coverAll <- read.csv('e141_plant species percent cover.csv')

#treatment data
trt <- read.csv('e141_treatments.csv')

#climate data
climate <- read.csv('e080_Daily climate summary.csv')

