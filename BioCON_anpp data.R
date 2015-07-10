setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_treatment data.r', chdir=T)

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

#rename Green.Biomass as Unsorted.Biomass (inconsistantly named across years and seasons, but never both used at once)
anppAll$species <- as.character(anppAll$species)
anppAll$species2 <- as.character(ifelse(anppAll$species=='Green.Biomass', 'Unsorted.Biomass', anppAll$species))

#anpp for the 16 BioCON species only (i.e., removing weeds, litter, etc.)
anpp16 <- anppAll[(anppAll$species %in% c("Achillea.millefolium","Asclepias.tuberosa",
                                          "Koeleria.cristata","Lupinus.perennis","Poa.pratensis",
                                          "Sorghastrum.nutans","Agropyron.repens","Andropogon.gerardi",
                                          "Bromus.inermis","Petalostemum.villosum","Amorpha.canescens",
                                          "Bouteloua.gracilis","Schizachyrium.scoparium", "Solidago.rigida",
                                          "Anemone.cylindrica", "Lespedeza.capitata", "Unsorted.Biomass",
                                          "Green.Biomass")),]

#get species anpp as columns
anpp <- dcast(anpp16, sampling_num + year + month + plot + ring + CO2_trt + N_trt + spp_count + group_count + experiment + monospecies + monogroup + water_trt + temp_trt ~ species2, value.var='anpp')
anpp[is.na(anpp)] <- 0
names(anpp)[names(anpp)=="Achillea.millefolium"] <- "ACMI"
names(anpp)[names(anpp)=="Asclepias.tuberosa"] <- "ASTU"
names(anpp)[names(anpp)=="Koeleria.cristata"] <- "KOCR"
names(anpp)[names(anpp)=="Lupinus.perennis"] <- "LUPE"
names(anpp)[names(anpp)=="Poa.pratensis"] <- "POPR"
names(anpp)[names(anpp)=="Sorghastrum.nutans"] <- "SONU"
names(anpp)[names(anpp)=="Agropyron.repens"] <- "AGRE"
names(anpp)[names(anpp)=="Andropogon.gerardi"] <- "ANGE"
names(anpp)[names(anpp)=="Bromus.inermis"] <- "BRIN"
names(anpp)[names(anpp)=="Petalostemum.villosum"] <- "PEVI"
names(anpp)[names(anpp)=="Amorpha.canescens"] <- "AMCA"
names(anpp)[names(anpp)=="Bouteloua.gracilis"] <- "BOGR"
names(anpp)[names(anpp)=="Schizachyrium.scoparium"] <- "SCSC"
names(anpp)[names(anpp)=="Solidago.rigida"] <- "SORI"
names(anpp)[names(anpp)=="Anemone.cylindrica"] <- "ANCY"
names(anpp)[names(anpp)=="Lespedeza.capitata"] <- "LECA"
names(anpp)[names(anpp)=="Unsorted.Biomass"] <- "Unsort"

#merge anpp and trt data
anppTrt <- merge(anpp, trt, all=T)
























