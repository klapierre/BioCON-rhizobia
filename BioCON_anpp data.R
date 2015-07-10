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

#get rid of 16 spp weeds in plots where they are labeled as spp
anppTrt$ACMI1 <- with(anppTrt, ifelse(Achillea.millefolium==1, ACMI, 0))
anppTrt$AGRE1 <- with(anppTrt, ifelse(Agropyron.repens==1, AGRE, 0))
anppTrt$AMCA1 <- with(anppTrt, ifelse(Amorpha.canescens==1, AMCA, 0))
anppTrt$ANGE1 <- with(anppTrt, ifelse(Andropogon.gerardi==1, ANGE, 0))
anppTrt$ANCY1 <- with(anppTrt, ifelse(Anemone.cylindrica==1, ANCY, 0))
anppTrt$ASTU1 <- with(anppTrt, ifelse(Asclepias.tuberosa==1, ASTU, 0))
anppTrt$BOGR1 <- with(anppTrt, ifelse(Bouteloua.gracilis==1, BOGR, 0))
anppTrt$BRIN1 <- with(anppTrt, ifelse(Bromus.inermis==1, BRIN, 0))
anppTrt$KOCR1 <- with(anppTrt, ifelse(Koeleria.cristata==1, KOCR, 0))
anppTrt$LECA1 <- with(anppTrt, ifelse(Lespedeza.capitata==1, LECA, 0))
anppTrt$LUPE1 <- with(anppTrt, ifelse(Lupinus.perennis==1, LUPE, 0))
anppTrt$PEVI1 <- with(anppTrt, ifelse(Petalostemum.villosum==1, PEVI, 0))
anppTrt$POPR1 <- with(anppTrt, ifelse(Poa.pratensis==1, POPR, 0))
anppTrt$SCSC1 <- with(anppTrt, ifelse(Schizachyrium.scoparium==1, SCSC, 0))
anppTrt$SORI1 <- with(anppTrt, ifelse(Solidago.rigida==1, SORI, 0))
anppTrt$SONU1 <- with(anppTrt, ifelse(Sorghastrum.nutans==1, SONU, 0))
anppNoTrt <- subset(anppTrt, select=-c(Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, CO2_trt, N_trt, spp_count, group_count, experiment, monospecies, monogroup, sampling_num, water_trt, temp_trt, C.3, C.4, Forb, Legume, legume_spp, legume_num, trt, spp_trt, ACMI, ASTU, KOCR, LUPE, POPR, SONU, AGRE, ANGE, BRIN, PEVI, AMCA, BOGR, SCSC, SORI, ANCY, LECA))
names(anppNoTrt)[names(anppNoTrt)=="ACMI1"] <- "Achillea.millefolium"
names(anppNoTrt)[names(anppNoTrt)=="ASTU1"] <- "Asclepias.tuberosa"
names(anppNoTrt)[names(anppNoTrt)=="KOCR1"] <- "Koeleria.cristata"
names(anppNoTrt)[names(anppNoTrt)=="LUPE1"] <- "Lupinus.perennis"
names(anppNoTrt)[names(anppNoTrt)=="POPR1"] <- "Poa.pratensis"
names(anppNoTrt)[names(anppNoTrt)=="SONU1"] <- "Sorghastrum.nutans"
names(anppNoTrt)[names(anppNoTrt)=="AGRE1"] <- "Agropyron.repens"
names(anppNoTrt)[names(anppNoTrt)=="ANGE1"] <- "Andropogon.gerardi"
names(anppNoTrt)[names(anppNoTrt)=="BRIN1"] <- "Bromus.inermis"
names(anppNoTrt)[names(anppNoTrt)=="PEVI1"] <- "Petalostemum.villosum"
names(anppNoTrt)[names(anppNoTrt)=="AMCA1"] <- "Amorpha.canescens"
names(anppNoTrt)[names(anppNoTrt)=="BOGR1"] <- "Bouteloua.gracilis"
names(anppNoTrt)[names(anppNoTrt)=="SCSC1"] <- "Schizachyrium.scoparium"
names(anppNoTrt)[names(anppNoTrt)=="SORI1"] <- "Solidago.rigida"
names(anppNoTrt)[names(anppNoTrt)=="ANCY1"] <- "Anemone.cylindrica"
names(anppNoTrt)[names(anppNoTrt)=="LECA1"] <- "Lespedeza.capitata"
names(anppNoTrt)[names(anppNoTrt)=="Unsort"] <- "Unsorted.Biomass"
anppTrue16 <- melt(anppNoTrt, id.vars=c('plot', 'ring', 'year', 'month'), variable.name='species', value.name='anpp')

#merge again with trt data
anppTrue16Trt <- merge(anppTrue16, trt)

#fixing problem with "unsorted" and "green" biomass in monocultures not being labeled as the monospecies; Peter okayed this assumption Apr 28 2015
anppTrue16Trt$species <- as.character(anppTrue16Trt$species)
anppTrue16Trt$monospecies <- as.character(anppTrue16Trt$monospecies)
anppTrue16TrtNonzero <- subset(anppTrue16Trt, subset=(anpp>0))
anppTrue16TrtNonzero$fix <- ifelse(anppTrue16TrtNonzero$species=='Unsorted.Biomass', 1, 0)
anppTrue16TrtNonzero$fix2 <- anppTrue16TrtNonzero$spp_count+anppTrue16TrtNonzero$fix
anppTrue16TrtNonzero$species <- as.factor(with(anppTrue16TrtNonzero, ifelse(fix2==2, monospecies, species)))

#take max of spring and fall biomass values for each species in each plot and year
anppMax <- ddply(anppTrue16TrtNonzero, c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup', 'species', 'Achillea.millefolium', 'Agropyron.repens', 'Amorpha.canescens', 'Andropogon.gerardi', 'Anemone.cylindrica', 'Asclepias.tuberosa', 'Bouteloua.gracilis', 'Bromus.inermis', 'Koeleria.cristata', 'Lespedeza.capitata', 'Lupinus.perennis', 'Petalostemum.villosum', 'Poa.pratensis', 'Schizachyrium.scoparium', 'Solidago.rigida', 'Sorghastrum.nutans', 'C.3', 'C.4', 'Forb', 'Legume', 'legume_num', 'legume_spp', 'trt', 'spp_trt'), summarise, anpp=max(anpp))
















