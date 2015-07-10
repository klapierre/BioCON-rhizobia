setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

#ANPP data from 1998-2012 (from web; IMPORTANT: be careful that extra spaces are removed from after species names if data is newly downloaded from web)
anppInitial <- read.csv('e141_plant aboveground biomass_1998-2012.csv')

#ANPP data from 2013-2014 (from Peter Reich)
anppLast <- read.csv('e141_plant aboveground biomass_2013-2014.csv')

#get year and month for 1998-2012 data
anppInitial$date <- as.Date(as.character(anppInitial$date), format='%m/%d/%Y') #make date column into a date
anppInitial$year <- as.numeric(format(anppInitial$date, '%Y')) #year as numeric
anppInitial$month <- as.factor(format(anppInitial$date, '%B')) #month in text format
anppInitialYear <- anppInitial[,-2]

#transpose 2013-2014 data to match 1998-2012 data structure
anppLastLong <- melt(anppLast, id.vars=c('sampling_num', 'year', 'month', 'ring', 'plot', 'CO2_trt', 'N_trt', 'water_trt', 'temp_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup'), variable.name='species', value.name='anpp')
anppLastLongComplete <- anppLastLong[complete.cases(anppLastLong[,16]),]

#append 2013-2014 data to the initial dataset
anppAll <- rbind(anppInitialYear, anppLastLongComplete)

