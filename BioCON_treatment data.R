setwd('C:\\Users\\Kim\\Dropbox\\NSF BioCON rhizobia\\data\\BioCON data')

#treatment data
trt <- read.csv('e141_treatments.csv')

#generate a column with binary code for the species treatment for each plot
trt$spp_trt <- with(trt, as.factor(paste(Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, sep='_')))

#get total number of legumes in each plot
trt$legume_num <- trt$Amorpha.canescens + trt$Lespedeza.capitata + trt$Lupinus.perennis + trt$Petalostemum.villosum

#make a column with legume ID for plots with one legume
trt$legume_spp <- as.factor(paste(trt$Amorpha.canescens, trt$Lespedeza.capitata, trt$Lupinus.perennis, trt$Petalostemum.villosum, sep='_'))

#make a column with legume spp for one legume plots, "two", "three" or "all 4" for multiple legume plots
trt$leg_num_spp_initial <- with(trt, as.character(ifelse(Amorpha.canescens==1, 'AMCA',
                                              ifelse(Lespedeza.capitata==1, 'LECA',
                                                     ifelse(Lupinus.perennis==1, 'LUPE',
                                                            ifelse(Petalostemum.villosum==1, 'PEVI', '0 legumes'))))))
trt$leg_num_spp <- with(trt, as.factor(ifelse(legume_num==2, '2 legumes',
                                       ifelse(legume_num==3, '3 legumes',
                                       ifelse(legume_num==4, '4 legumes', leg_num_spp_initial)))))
trt <- trt[,-33]

#making a combined CO2 and N treatment column
trt$trt <- paste(trt$CO2_trt, trt$N_trt, sep='_')
