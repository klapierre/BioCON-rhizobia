################################################################################
##  1_BioCON_treatment data.R: Gathering BioCON plot treatment information.
##
##  Author: Kimberly Komatsu
################################################################################

library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

#### Read in treatment data ####
trt <- read.csv('e141_treatments.csv') %>% 
  #get total number of legumes in each plot
  mutate(legume_num=Amorpha.canescens + Lespedeza.capitata + Lupinus.perennis + Petalostemum.villosum) %>% 
  #make a column with legume spp for one legume plots, "two", "three" or "all 4" for multiple legume plots
  mutate(leg_num_spp=as.character(ifelse(legume_num==1 & Amorpha.canescens==1, 'AMCA',                                                 ifelse(legume_num==1 & Lespedeza.capitata==1, 'LECA',
                                  ifelse(legume_num==1 & Lupinus.perennis==1, 'LUPE',
                                  ifelse(legume_num==1 & Petalostemum.villosum==1, 'PEVI',
                                  ifelse(legume_num==2, '2 legumes',
                                  ifelse(legume_num==3, '3 legumes',
                                  ifelse(legume_num==4, '4 legumes', '0 legumes'))))))))) %>% 
  #making a combined CO2 and N treatment column
  mutate(trt=paste(CO2_trt, N_trt, sep='_')) %>% 
  #generate column with list of focal species in each plot
  mutate(Achillea.millefolium=ifelse(Achillea.millefolium==1, "Achillea millefolium", 0),
         Agropyron.repens=ifelse(Agropyron.repens==1, "Agropyron repens", 0),
         Amorpha.canescens=ifelse(Amorpha.canescens==1, "Amorpha canescens", 0),
         Andropogon.gerardi=ifelse(Andropogon.gerardi==1, "Andropogon gerardi", 0),
         Anemone.cylindrica=ifelse(Anemone.cylindrica==1, "Anemone cylindrica", 0),
         Asclepias.tuberosa=ifelse(Asclepias.tuberosa==1, "Asclepias tuberosa", 0),
         Bouteloua.gracilis=ifelse(Bouteloua.gracilis==1, "Bouteloua gracilis", 0),
         Bromus.inermis=ifelse(Bromus.inermis==1, "Bromus inermis", 0),
         Koeleria.cristata=ifelse(Koeleria.cristata==1, "Koeleria cristata", 0),
         Lespedeza.capitata=ifelse(Lespedeza.capitata==1, "Lespedeza capitata", 0),
         Lupinus.perennis=ifelse(Lupinus.perennis==1, "Lupinus perennis", 0),
         Petalostemum.villosum=ifelse(Petalostemum.villosum==1, "Petalostemum villosum", 0),
         Poa.pratensis=ifelse(Poa.pratensis==1, "Poa pratensis", 0),
         Schizachyrium.scoparium=ifelse(Schizachyrium.scoparium==1, "Schizachyrium scoparium", 0),
         Solidago.rigida=ifelse(Solidago.rigida==1, "Solidago rigida", 0),
         Sorghastrum.nutans=ifelse(Sorghastrum.nutans==1, "Sorghastrum nutans", 0)) %>% 
  mutate(spp_trt=paste(Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, sep=', ')) %>% 
  #make a column with legume ID for plots with one legume
  mutate(legume_spp=as.factor(paste(Amorpha.canescens, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, sep='_')))