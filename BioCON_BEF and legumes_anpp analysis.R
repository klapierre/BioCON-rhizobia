library(plyr)
library(ggplot2)
library(reshape2)
library(grid)
library(nlme) 
library(lsmeans)
library(lavaan)

source('C:\\Users\\Kim\\Desktop\\general R code\\general-R-code\\bar graph summary stats.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\general R code\\general-R-code\\ggplot_theme set.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_treatment data.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_climate data.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_anpp data.r', chdir=T)
source('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data\\BioCON-rhizobia\\BioCON_cover data.r', chdir=T)

setwd('C:\\Users\\Kim\\Desktop\\BioCON rhizobia\\BioCON data')

#merge ANPP and climate data
anppClimate <- merge(anppRY, climate)

#mixed model
#IMPORTANT: 4 spp polycultures
#IMPORTANT: using growing year climate variables
RYTmixedCategorical <- lme(logRYT ~ gy_precip_cm*leg_num_spp, random=~1|plot, data=subset(anppClimate, spp_count==4 & trt=='Camb_Namb'))
summary(RYTmixedCategorical)
anova(RYTmixedCategorical)
lsmeans(RYTmixedCategorical, cld~leg_num_spp)

anppClimate$order <- factor(anppClimate$leg_num_spp, levels=c('0 legumes', 'AMCA', 'LECA', 'LUPE', 'PEVI', '2 legumes', '3 legumes', '4 legumes'))

ggplot(subset(anppClimate, spp_count==4 & trt=='Camb_Namb'), aes(x=gy_precip_cm, y=logRYT, colour=order)) +
  geom_point() +
  geom_smooth(method=lm) +
  ylab('log Relative Yield Total') +
  xlab('Total Annual Precipitation (cm)')


