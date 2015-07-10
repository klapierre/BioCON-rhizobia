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




