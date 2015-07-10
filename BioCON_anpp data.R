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

#remove temp and water manipulation data, bareground plots
anppSub <- subset(anppTrt, subset=(temp_trt!='HTamb' & temp_trt!='HTelv' & water_trt!='H2Oamb' & water_trt!='H2Oneg' & spp_count!=0))

#get rid of 16 spp weeds in plots where they are labeled as spp
anppSub$ACMI1 <- with(anppSub, ifelse(Achillea.millefolium==1, ACMI, 0))
anppSub$AGRE1 <- with(anppSub, ifelse(Agropyron.repens==1, AGRE, 0))
anppSub$AMCA1 <- with(anppSub, ifelse(Amorpha.canescens==1, AMCA, 0))
anppSub$ANGE1 <- with(anppSub, ifelse(Andropogon.gerardi==1, ANGE, 0))
anppSub$ANCY1 <- with(anppSub, ifelse(Anemone.cylindrica==1, ANCY, 0))
anppSub$ASTU1 <- with(anppSub, ifelse(Asclepias.tuberosa==1, ASTU, 0))
anppSub$BOGR1 <- with(anppSub, ifelse(Bouteloua.gracilis==1, BOGR, 0))
anppSub$BRIN1 <- with(anppSub, ifelse(Bromus.inermis==1, BRIN, 0))
anppSub$KOCR1 <- with(anppSub, ifelse(Koeleria.cristata==1, KOCR, 0))
anppSub$LECA1 <- with(anppSub, ifelse(Lespedeza.capitata==1, LECA, 0))
anppSub$LUPE1 <- with(anppSub, ifelse(Lupinus.perennis==1, LUPE, 0))
anppSub$PEVI1 <- with(anppSub, ifelse(Petalostemum.villosum==1, PEVI, 0))
anppSub$POPR1 <- with(anppSub, ifelse(Poa.pratensis==1, POPR, 0))
anppSub$SCSC1 <- with(anppSub, ifelse(Schizachyrium.scoparium==1, SCSC, 0))
anppSub$SORI1 <- with(anppSub, ifelse(Solidago.rigida==1, SORI, 0))
anppSub$SONU1 <- with(anppSub, ifelse(Sorghastrum.nutans==1, SONU, 0))
anppNoTrt <- subset(anppSub, select=-c(Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, CO2_trt, N_trt, spp_count, group_count, experiment, monospecies, monogroup, sampling_num, water_trt, temp_trt, C.3, C.4, Forb, Legume, legume_spp, legume_num, leg_num_spp, trt, spp_trt, ACMI, ASTU, KOCR, LUPE, POPR, SONU, AGRE, ANGE, BRIN, PEVI, AMCA, BOGR, SCSC, SORI, ANCY, LECA))
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
anppMax <- ddply(anppTrue16TrtNonzero, c('year', 'plot', 'ring', 'CO2_trt', 'N_trt', 'spp_count', 'group_count', 'experiment', 'monospecies', 'monogroup', 'species', 'Achillea.millefolium', 'Agropyron.repens', 'Amorpha.canescens', 'Andropogon.gerardi', 'Anemone.cylindrica', 'Asclepias.tuberosa', 'Bouteloua.gracilis', 'Bromus.inermis', 'Koeleria.cristata', 'Lespedeza.capitata', 'Lupinus.perennis', 'Petalostemum.villosum', 'Poa.pratensis', 'Schizachyrium.scoparium', 'Solidago.rigida', 'Sorghastrum.nutans', 'C.3', 'C.4', 'Forb', 'Legume', 'legume_num', 'legume_spp', 'leg_num_spp', 'trt', 'spp_trt'), summarise, anpp=max(anpp))

#subset out monocultures
anppMono <- subset(anppMax, subset=(spp_count==1 & anpp!=0))
#subset out polycultures
anppPoly <- subset(anppMax, subset=(spp_count!=1 & anpp!=0))

#calculate average biomass of each species in monoculture by year, CO2 and N trt
anppMonoAvg <- aggregate(anppMono$anpp, by=list(year=anppMono$year, CO2_trt=anppMono$CO2_trt, N_trt=anppMono$N_trt, species=anppMono$species, trt=anppMono$trt), FUN=mean)
names(anppMonoAvg)[names(anppMonoAvg)=="x"] <- "mono_anpp"

#make species biomass as columns
anppMonoAvgWide <- dcast(anppMonoAvg, year + CO2_trt + N_trt + trt ~ species, value.var='mono_anpp')

names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Achillea.millefolium"] <- "ACMImono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Agropyron.repens"] <- "AGREmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Amorpha.canescens"] <- "AMCAmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Andropogon.gerardi"] <- "ANGEmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Anemone.cylindrica"] <- "ANCYmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Asclepias.tuberosa"] <- "ASTUmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Bouteloua.gracilis"] <- "BOGRmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Bromus.inermis"] <- "BRINmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Koeleria.cristata"] <- "KOCRmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Lespedeza.capitata"] <- "LECAmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Lupinus.perennis"] <- "LUPEmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Petalostemum.villosum"] <- "PEVImono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Poa.pratensis"] <- "POPRmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Schizachyrium.scoparium"] <- "SCSCmono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Solidago.rigida"] <- "SORImono"
names(anppMonoAvgWide)[names(anppMonoAvgWide)=="Sorghastrum.nutans"] <- "SONUmono"

#calculate total plot biomass
anppSum <- with(anppPoly, aggregate(anpp, by=list(year=year, plot=plot, ring=ring, CO2_trt=CO2_trt, N_trt=N_trt, spp_count=spp_count, group_count=group_count, experiment=experiment, monogroup=monogroup, Achillea.millefolium=Achillea.millefolium, Agropyron.repens=Agropyron.repens, Amorpha.canescens=Amorpha.canescens, Andropogon.gerardi=Andropogon.gerardi, Anemone.cylindrica=Anemone.cylindrica, Asclepias.tuberosa=Asclepias.tuberosa, Bouteloua.gracilis=Bouteloua.gracilis, Bromus.inermis=Bromus.inermis, Koeleria.cristata=Koeleria.cristata, Lespedeza.capitata=Lespedeza.capitata, Lupinus.perennis=Lupinus.perennis, Petalostemum.villosum=Petalostemum.villosum, Poa.pratensis=Poa.pratensis, Schizachyrium.scoparium=Schizachyrium.scoparium, Solidago.rigida=Solidago.rigida, Sorghastrum.nutans=Sorghastrum.nutans, C.3=C.3, C.4=C.4, Forb=Forb, Legume=Legume, legume_spp=legume_spp, legume_num=legume_num, leg_num_spp=leg_num_spp, trt=trt, spp_trt=spp_trt), FUN=sum))
names(anppSum)[names(anppSum)=="x"] <- "total_biomass"

#merge monoculture averages with polyculture data
anppRY <- merge(anppSum, anppMonoAvgWide, all=T)

#calculate expected yield from monoculture biomass
anppRY$ACMImono2 <- with(anppRY, ifelse(Achillea.millefolium==1, ACMImono, 0))
anppRY$AGREmono2 <- with(anppRY, ifelse(Agropyron.repens==1, AGREmono, 0))
anppRY$AMCAmono2 <- with(anppRY, ifelse(Amorpha.canescens==1, AMCAmono, 0))
anppRY$ANGEmono2 <- with(anppRY, ifelse(Andropogon.gerardi==1, ANGEmono, 0))
anppRY$ANCYmono2 <- with(anppRY, ifelse(Anemone.cylindrica==1, ANCYmono, 0))
anppRY$ANCYmono2[is.na(anppRY$ANCYmono2)] <- 0
anppRY$ASTUmono2 <- with(anppRY, ifelse(Asclepias.tuberosa==1, ASTUmono, 0))
anppRY$BOGRmono2 <- with(anppRY, ifelse(Bouteloua.gracilis==1, BOGRmono, 0))
anppRY$BRINmono2 <- with(anppRY, ifelse(Bromus.inermis==1, BRINmono, 0))
anppRY$KOCRmono2 <- with(anppRY, ifelse(Koeleria.cristata==1, KOCRmono, 0))
anppRY$LECAmono2 <- with(anppRY, ifelse(Lespedeza.capitata==1, LECAmono, 0))
anppRY$LUPEmono2 <- with(anppRY, ifelse(Lupinus.perennis==1, LUPEmono, 0))
anppRY$PEVImono2 <- with(anppRY, ifelse(Petalostemum.villosum==1, PEVImono, 0))
anppRY$POPRmono2 <- with(anppRY, ifelse(Poa.pratensis==1, POPRmono, 0))
anppRY$SCSCmono2 <- with(anppRY, ifelse(Schizachyrium.scoparium==1, SCSCmono, 0))
anppRY$SORImono2 <- with(anppRY, ifelse(Solidago.rigida==1, SORImono, 0))
anppRY$SONUmono2 <- with(anppRY, ifelse(Sorghastrum.nutans==1, SONUmono, 0))
anppRY$exp_yield <- with(anppRY, ifelse(spp_count==4, (ACMImono2+AGREmono2+AMCAmono2+ANGEmono2+ANCYmono2+ASTUmono2+BOGRmono2+BRINmono2+KOCRmono2+LECAmono2+LUPEmono2+PEVImono2+POPRmono2+SCSCmono2+SORImono2+SONUmono2)/4, 
                                        ifelse(spp_count==9, (ACMImono2+AGREmono2+AMCAmono2+ANGEmono2+ANCYmono2+ASTUmono2+BOGRmono2+BRINmono2+KOCRmono2+LECAmono2+LUPEmono2+PEVImono2+POPRmono2+SCSCmono2+SORImono2+SONUmono2)/9, (ACMImono2+AGREmono2+AMCAmono2+ANGEmono2+ANCYmono2+ASTUmono2+BOGRmono2+BRINmono2+KOCRmono2+LECAmono2+LUPEmono2+PEVImono2+POPRmono2+SCSCmono2+SORImono2+SONUmono2)/16)))

#calculate relative yield total
anppRY$RYT <- with(anppRY, total_biomass/exp_yield)

#check RYT for normality within plots
shapiro.test(anppRY$RYT)
qqnorm(anppRY$RYT)

#make normal
anppRY$logRYT <- log10(anppRY$RYT)
shapiro.test(anppRY$logRYT)
qqnorm(anppRY$logRYT)

# anppRY$sqrtRYT <- sqrt(anppRY$RYT)
# shapiro.test(anppRY$sqrtRYT)
# qqnorm(anppRY$sqrtRYT)

##############################################
##############################################

#RY of legumes in monoculture vs polyculture, with and without other legumes
anppPolyRY <- merge(anppPoly, anppMonoAvg, all=T)
anppPolyRY$RY <- with(anppPolyRY, ifelse(spp_count==4, (anpp/(mono_anpp/4)),
                                         ifelse(spp_count==9, (anpp/(mono_anpp/9)), (anpp/(mono_anpp/16)))))

anppPolyRY <- anppPolyRY[complete.cases(anppPolyRY),]

#make column identifying if a species is a legume or non-legume
anppPolyRY$spp_type <- with(anppPolyRY, ifelse(species=='Amorpha.canescens', 'legume',
                                               ifelse(species=='Lespedeza.capitata', 'legume',
                                                      ifelse(species=='Lupinus.perennis', 'legume',
                                                             ifelse(species=='Petalostemum.villosum', 'legume', 'non-legume')))))

anppPolyRY$other_legs <- with(anppPolyRY, as.factor(ifelse(leg_num_spp=='0 legumes', '0 legumes',
                                                           ifelse(leg_num_spp=='2 legumes', '2 legumes',
                                                                  ifelse(leg_num_spp=='4 legumes', '4 legumes', 'single legume')))))

#normalize data

qqnorm(anppPolyRY$RY)
anppPolyRY$logRY <- log10(anppPolyRY$RY)
qqnorm(anppPolyRY$logRY)

##############################################
##############################################

#scale up RY to spp_num*RY
anppPolyRY$RYscaled <- with(anppPolyRY, spp_count*RY)
anppPolyRY$RYscaled_log <- log10(anppPolyRY$RYscaled)

#normality for anpp
qqnorm(anppPolyRY$anpp)
anppPolyRY$log_anpp <- log10(anppPolyRY$anpp)
qqnorm(anppPolyRY$log_anpp)

##############################################
##############################################

#difference between 1 and 4 legume plots
legume4 <- subset(anppPolyRY, spp_count==4 & other_legs=='4 legumes')
legume4mean <- with(legume4, aggregate(RY, list(species=species, year=year, CO2_trt=CO2_trt, N_trt=N_trt, trt=trt), mean))
names(legume4mean)[names(legume4mean)=='x'] <- 'RY4'
legume4mean$logRY4 <- log10(legume4mean$RY4)

legume1 <- subset(anppPolyRY, spp_count==4 & other_legs=='single legume')
legume1mean <- with(legume1, aggregate(RY, list(species=species, year=year, CO2_trt=CO2_trt, N_trt=N_trt, trt=trt), mean))
names(legume1mean)[names(legume1mean)=='x'] <- 'RY1'
legume1mean$logRY1 <- log10(legume1mean$RY1)

legume4to1 <- merge(legume4mean, legume1mean, all=T)
legume4to1complete <- legume4to1[complete.cases(legume4to1[,c(7,9)]),]
legume4to1complete$diff <- legume4to1complete$RY4-legume4to1complete$RY1
legume4to1complete$diff_log <- legume4to1complete$logRY4-legume4to1complete$logRY1

shapiro.test(legume4to1complete$diff)
qqnorm(legume4to1complete$diff)

shapiro.test(legume4to1complete$diff_log)
qqnorm(legume4to1complete$diff_log)

##############################################
##############################################

#get anpp hedges d
anppTrtMeans <- ddply(subset(anppPolyRY, spp_count==4), c('year', 'spp_type', 'trt'), summarise,
                      anpp_mean=mean(anpp),
                      anpp_sd=sd(anpp),
                      anpp_n=length(anpp))

anppRRmean <- dcast(anppTrtMeans, year + spp_type ~ trt, value.var='anpp_mean')
anppRRsd <- dcast(anppTrtMeans, year + spp_type ~ trt, value.var='anpp_sd')
names(anppRRsd)[names(anppRRsd)=='Camb_Namb'] <- 'ctl_sd'
names(anppRRsd)[names(anppRRsd)=='Camb_Nenrich'] <- 'N_sd'
names(anppRRsd)[names(anppRRsd)=='Cenrich_Namb'] <- 'CO2_sd'
names(anppRRsd)[names(anppRRsd)=='Cenrich_Nenrich'] <- 'CO2_N_sd'
anppRRn <- dcast(anppTrtMeans, year + spp_type ~ trt, value.var='anpp_n')
names(anppRRn)[names(anppRRn)=='Camb_Namb'] <- 'ctl_n'
names(anppRRn)[names(anppRRn)=='Camb_Nenrich'] <- 'N_n'
names(anppRRn)[names(anppRRn)=='Cenrich_Namb'] <- 'CO2_n'
names(anppRRn)[names(anppRRn)=='Cenrich_Nenrich'] <- 'CO2_N_n'
anppRRmeansd <- merge(anppRRmean, anppRRsd)
anppRR <- merge(anppRRmeansd, anppRRn)
anppRR$hedgesd_CO2 <- with(anppRR, ((Cenrich_Namb-Camb_Namb)/sqrt((CO2_sd*(CO2_n-1)+ctl_sd*(ctl_n-1))/(CO2_n+ctl_n-2)))*(1-(3/(4*(CO2_n+ctl_n-2)-1))))
anppRR$hedgesd_N <- with(anppRR, ((Camb_Nenrich-Camb_Namb)/sqrt((N_sd*(N_n-1)+ctl_sd*(ctl_n-1))/(N_n+ctl_n-2)))*(1-(3/(4*(N_n+ctl_n-2)-1))))
anppRR$hedgesd_CO2_N <- with(anppRR, ((Cenrich_Nenrich-Camb_Namb)/sqrt((CO2_N_sd*(CO2_N_n-1)+ctl_sd*(ctl_n-1))/(CO2_N_n+ctl_n-2)))*(1-(3/(4*(CO2_N_n+ctl_n-2)-1))))
anppRR$hedgesd_var_CO2 <- with(anppRR, ((CO2_n+ctl_n)/(CO2_n*ctl_n))+((hedgesd_CO2^2)/(2*(CO2_n+ctl_n))))
anppRR$hedgesd_var_N <- with(anppRR, ((N_n+ctl_n)/(N_n*ctl_n))+((hedgesd_N^2)/(2*(N_n+ctl_n))))
anppRR$hedgesd_var_CO2_N <- with(anppRR, ((CO2_N_n+ctl_n)/(CO2_N_n*ctl_n))+((hedgesd_CO2_N^2)/(2*(CO2_N_n+ctl_n))))
anppRR$hedgesd_ci_CO2 <- 1.96*anppRR$hedgesd_var_CO2
anppRR$hedgesd_ci_N <- 1.96*anppRR$hedgesd_var_N
anppRR$hedgesd_ci_CO2_N <- 1.96*anppRR$hedgesd_var_CO2_N

##############################################
##############################################

#get anpp hedges d by legume spp
anppTrtMeansSpp <- ddply(subset(anppPolyRY, spp_count==4 & spp_type=='legume'), c('year', 'species', 'trt'), summarise,
                         anpp_mean=mean(anpp),
                         anpp_sd=sd(anpp),
                         anpp_n=length(anpp))

anppRRmeanSpp <- dcast(anppTrtMeansSpp, year + species ~ trt, value.var='anpp_mean')
anppRRsdSpp <- dcast(anppTrtMeansSpp, year + species ~ trt, value.var='anpp_sd')
names(anppRRsdSpp)[names(anppRRsdSpp)=='Camb_Namb'] <- 'ctl_sd'
names(anppRRsdSpp)[names(anppRRsdSpp)=='Camb_Nenrich'] <- 'N_sd'
names(anppRRsdSpp)[names(anppRRsdSpp)=='Cenrich_Namb'] <- 'CO2_sd'
names(anppRRsdSpp)[names(anppRRsdSpp)=='Cenrich_Nenrich'] <- 'CO2_N_sd'
anppRRnSpp <- dcast(anppTrtMeansSpp, year + species ~ trt, value.var='anpp_n')
names(anppRRnSpp)[names(anppRRnSpp)=='Camb_Namb'] <- 'ctl_n'
names(anppRRnSpp)[names(anppRRnSpp)=='Camb_Nenrich'] <- 'N_n'
names(anppRRnSpp)[names(anppRRnSpp)=='Cenrich_Namb'] <- 'CO2_n'
names(anppRRnSpp)[names(anppRRnSpp)=='Cenrich_Nenrich'] <- 'CO2_N_n'
anppRRmeansdSpp <- merge(anppRRmeanSpp, anppRRsdSpp)
anppRRSpp <- merge(anppRRmeansdSpp, anppRRnSpp)
anppRRSpp$hedgesd_CO2 <- with(anppRRSpp, ((Cenrich_Namb-Camb_Namb)/sqrt((CO2_sd*(CO2_n-1)+ctl_sd*(ctl_n-1))/(CO2_n+ctl_n-2)))*(1-(3/(4*(CO2_n+ctl_n-2)-1))))
anppRRSpp$hedgesd_N <- with(anppRRSpp, ((Camb_Nenrich-Camb_Namb)/sqrt((N_sd*(N_n-1)+ctl_sd*(ctl_n-1))/(N_n+ctl_n-2)))*(1-(3/(4*(N_n+ctl_n-2)-1))))
anppRRSpp$hedgesd_CO2_N <- with(anppRRSpp, ((Cenrich_Nenrich-Camb_Namb)/sqrt((CO2_N_sd*(CO2_N_n-1)+ctl_sd*(ctl_n-1))/(CO2_N_n+ctl_n-2)))*(1-(3/(4*(CO2_N_n+ctl_n-2)-1))))
anppRRSpp$hedgesd_var_CO2 <- with(anppRRSpp, ((CO2_n+ctl_n)/(CO2_n*ctl_n))+((hedgesd_CO2^2)/(2*(CO2_n+ctl_n))))
anppRRSpp$hedgesd_var_N <- with(anppRRSpp, ((N_n+ctl_n)/(N_n*ctl_n))+((hedgesd_N^2)/(2*(N_n+ctl_n))))
anppRRSpp$hedgesd_var_CO2_N <- with(anppRRSpp, ((CO2_N_n+ctl_n)/(CO2_N_n*ctl_n))+((hedgesd_CO2_N^2)/(2*(CO2_N_n+ctl_n))))
anppRRSpp$hedgesd_ci_CO2 <- 1.96*anppRRSpp$hedgesd_var_CO2
anppRRSpp$hedgesd_ci_N <- 1.96*anppRRSpp$hedgesd_var_N
anppRRSpp$hedgesd_ci_CO2_N <- 1.96*anppRRSpp$hedgesd_var_CO2_N

anppHedgesDaSpp <- anppRRSpp[,-c(3:14,18:23)]
anppHedgesDSppd <- melt(anppHedgesDaSpp, id.vars=c('year', 'species'), variable.name='trt', value.name='hedges_d')
anppHedgesDbSpp <- anppRRSpp[,-c(3:20)]
anppHedgesDSppci <- melt(anppHedgesDbSpp, id.vars=c('year', 'species'), variable.name='trt', value.name='hedges_d_ci')
anppHedgesDSppci$trt2 <- with(anppHedgesDSppci, ifelse(trt=='hedgesd_ci_CO2', 'hedgesd_CO2',
                                                       ifelse(trt=='hedgesd_ci_N', 'hedgesd_N', 'hedgesd_CO2_N')))
anppHedgesDSppci$trt <- NULL
names(anppHedgesDSppci)[names(anppHedgesDSppci)=='trt2'] <- 'trt'
anppHedgesDSpp <- merge(anppHedgesDSppd, anppHedgesDSppci)


##############################################
##############################################

#clean up workspace
rm(list=c('anpp', 'anpp16', 'anppAll', 'anppInitial', 'anppInitialYear', 'anppLast', 'anppLastLong', 'anppLastLongComplete', 'anppMax', 'anppMono', 'anppMonoAvgWide', 'anppMonoAvg', 'anppPoly',  'anppNoTrt', 'anppSum', 'anppTrt', 'anppTrue16', 'anppTrue16Trt', 'anppTrue16TrtNonzero', 'trt', 'anppSub', 'legume1', 'legume1mean', 'legume4', 'legume4mean', 'legume4to1', 'anppTrtMeans', 'anppRRmean', 'anppRRmeansd', 'anppRRn', 'anppRRsd', 'anppHedgesDaSpp', 'anppHedgesDbSpp', 'anppHedgesDSppci', 'anppHedgesDSppd', 'anppRRmeansdSpp', 'anppRRmeanSpp', 'anppRRnSpp', 'anppRRSpp', 'anppTrtMeansSpp'))
