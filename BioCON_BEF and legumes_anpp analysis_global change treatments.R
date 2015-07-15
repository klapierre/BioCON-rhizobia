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
anppClimateInitial <- merge(anppRY, climate)
polyClimate <- merge(anppPolyRY, climate)
RRclimate <- merge(anppRR, climate)
anppMonoClimateInitial <- merge(anppMono, climate)

#remove outliers
anppClimate <- subset(anppClimateInitial, RYT<8.5)
anppMonoClimate <- subset(anppMonoClimateInitial, anpp<1100)

#mixed model
#do legumes or non-legumes respond more to CO2 and N?
#this is just species RY in 4 spp polyculture, not considering what type of polyculture they were in (BIC and AIC were lower when year and/or gy_precip were left out of model)
trtType4sppMixed <- lme(RYscaled_log ~ CO2_trt*N_trt*spp_type, random=~1|plot, data=subset(polyClimate, spp_count==4))
summary(trtType4sppMixed)
anova(trtType4sppMixed)

color <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")
ggplot(barGraphStats(data=subset(polyClimate, spp_count==4), variable='RYscaled_log', byFactorNames=c('spp_type', 'trt')), aes(x=spp_type, y=mean, fill=trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=color) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('log Relative Yield') +
  theme(axis.title.x=element_blank())


#mixed model
#do different legume species respond more to CO2 and N?
#this is just species biomass in monoculture
trtLegSppMixed <- lme(anpp ~ CO2_trt*N_trt*species, random=~1|plot, data=subset(anppMonoClimate, legume_num==1))
summary(trtLegSppMixed)
anova(trtLegSppMixed)

color <- c("#E69F00", "#009E73", "#0072B2", "#CC79A7")
ggplot(barGraphStats(data=subset(anppMonoClimate, legume_num==1), variable='anpp', byFactorNames=c('species', 'trt')), aes(x=species, y=mean, fill=trt)) +
  geom_bar(stat='identity', position=position_dodge()) +
  scale_fill_manual(values=color) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('ANPP (g m-2)') +
  theme(axis.title.x=element_blank())


#hedges' d of legumes vs non-legumes in 4 spp polycultures
#does hedgesd vary with gy precip?
summary(lm(hedgesd_CO2 ~ gy_precip_cm*spp_type, data=RRclimate)) #yes (marginal)
summary(lm(hedgesd_N ~ gy_precip_cm*spp_type, data=RRclimate)) #yes
summary(lm(hedgesd_CO2_N ~ gy_precip_cm*spp_type, data=RRclimate)) #yes

#mixed model with hedges d
anppHedgesDa <- RRclimate[,-c(5:16,20:25)]
anppHedgesD <- melt(anppHedgesDa, id.vars=c('year', 'gy_precip_cm', 'gy_temp', 'spp_type'), variable.name='trt', value.name='hedges_d')
hedgesdModel <- lm(hedges_d ~ year*trt*spp_type , data=anppHedgesD)
summary(hedgesdModel)
anova(hedgesdModel)

CO2plot <- ggplot(anppRR, aes(x=year, y=hedgesd_CO2, shape=spp_type)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, aes(linetype=spp_type), colour='black') +
  geom_errorbar(aes(ymin=hedgesd_CO2-hedgesd_ci_CO2, ymax=hedgesd_CO2+hedgesd_ci_CO2, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position=c(0.8,0.95), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-2.5,12), xlim=c(1997.5,2005.5)) +
  scale_y_continuous(breaks=seq(-2,12,2)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  annotate('text', x=1997.75, y=11.5, label='(a) CO2', size=7, hjust=0)
Nplot <- ggplot(anppRR, aes(x=year, y=hedgesd_N, shape=spp_type)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, aes(linetype=spp_type), colour='black') +
  geom_errorbar(aes(ymin=hedgesd_N-hedgesd_ci_N, ymax=hedgesd_N+hedgesd_ci_N, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-2.5,12), xlim=c(1997.5,2005.5)) +
  scale_y_continuous(breaks=seq(-2,12,2)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  annotate('text', x=1997.75, y=11.5, label='(b) N', size=7, hjust=0)
CO2Nplot <- ggplot(anppRR, aes(x=year, y=hedgesd_CO2_N, shape=spp_type)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, aes(linetype=spp_type), colour='black') +
  geom_errorbar(aes(ymin=hedgesd_CO2_N-hedgesd_ci_CO2_N, ymax=hedgesd_CO2_N+hedgesd_ci_CO2_N, width=0.2)) +
  ylab('Hedges d') +
  coord_cartesian(ylim=c(-2.5,12), xlim=c(1997.5,2005.5)) +
  scale_y_continuous(breaks=seq(-2,12,2)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45))+
  annotate('text', x=1997.75, y=11.5, label='(c) CO2 + N', size=7, hjust=0)

#put hedges d panels together
pushViewport(viewport(layout=grid.layout(1,3))) 
print(CO2plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(Nplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(CO2Nplot, vp=viewport(layout.pos.row=1, layout.pos.col=3))


#legume species' hedges d by specis
#mixed model with hedges d
hedgesdModelSpp <- lm(hedges_d ~ year*trt*species , data=anppHedgesDSpp)
summary(hedgesdModelSpp)
anova(hedgesdModelSpp)

amcaplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Amorpha.canescens'), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position=c(0.8,0.85), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  scale_colour_discrete(labels=c('eCO2', 'eN', 'eCO2 + eN')) +
  scale_shape_discrete(labels=c('eCO2', 'eN', 'eCO2 + eN')) +
  annotate('text', x=1997.75, y=18, label='(a) A. canescens', size=7, hjust=0)
lecaplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Lespedeza.capitata'), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  annotate('text', x=1997.75, y=18, label='(b) Le. capitata', size=7, hjust=0)
lupeplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Lupinus.perennis'), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  annotate('text', x=1997.75, y=18, label='(c) Lu. perennis', size=7, hjust=0)
peviplotSpp <- ggplot(subset(anppHedgesDSpp, species=='Petalostemum.villosum' & year<2005), aes(x=year, y=hedges_d, shape=trt, colour=trt)) +
  geom_point(size=4) +
  geom_smooth(method=lm, size=1, se=F) +
  #   geom_errorbar(aes(ymin=hedges_d-hedges_d_ci, ymax=hedges_d+hedges_d_ci, width=0.2)) +
  ylab('Hedges d') +
  theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  coord_cartesian(ylim=c(-15,20), xlim=c(1997.5,2005.5)) +
  scale_x_continuous(breaks=seq(1998,2005,1)) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  annotate('text', x=1997.75, y=18, label='(d) P. villosum', size=7, hjust=0)

#put hedges d panels together
pushViewport(viewport(layout=grid.layout(2,2))) 
print(amcaplotSpp, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(lecaplotSpp, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(lupeplotSpp, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(peviplotSpp, vp=viewport(layout.pos.row=2, layout.pos.col=2))





















