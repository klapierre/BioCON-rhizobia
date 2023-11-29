################################################################################
##  1_BioCON_treatment data.R: Gathering BioCON plot treatment information.
##
##  Author: Kimberly Komatsu
################################################################################

library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

#### Read in treatment data ####
trt <- read.csv('e141_treatments.csv') %>% 
  #generate column with binary code for species treatment in each plot
  mutate(spp_trt=paste(Achillea.millefolium, Agropyron.repens, Amorpha.canescens, Andropogon.gerardi, Anemone.cylindrica, Asclepias.tuberosa, Bouteloua.gracilis, Bromus.inermis, Koeleria.cristata, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, Poa.pratensis, Schizachyrium.scoparium, Solidago.rigida, Sorghastrum.nutans, sep='_')) %>% 
  #get total number of legumes in each plot
  mutate(legume_num=(Amorpha.canescens + Lespedeza.capitata + Lupinus.perennis + Petalostemum.villosum)) %>% 
  #make a column with legume ID for plots with one legume
  mutate(legume_spp=as.factor(paste(Amorpha.canescens, Lespedeza.capitata, Lupinus.perennis, Petalostemum.villosum, sep='_'))) %>% 
  #make a column with legume spp for one legume plots, "two", "three" or "all 4" for multiple legume plots
  mutate(leg_num_spp=as.character(ifelse(legume_num==1 & Amorpha.canescens==1, 'AMCA',                                                 ifelse(legume_num==1 & Lespedeza.capitata==1, 'LECA',
                                  ifelse(legume_num==1 & Lupinus.perennis==1, 'LUPE',
                                  ifelse(legume_num==1 & Petalostemum.villosum==1, 'PEVI',
                                  ifelse(legume_num==2, '2 legumes',
                                  ifelse(legume_num==3, '3 legumes',
                                  ifelse(legume_num==4, '4 legumes', '0 legumes'))))))))) %>% 
  #making a combined CO2 and N treatment column
  mutate(trt=paste(trt$CO2_trt, trt$N_trt, sep='_'))